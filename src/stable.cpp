#include "stable.hpp"
#include "stdio.h"
#include <string.h>
#include "matchingPoint.hpp"
#include "preprocess.hpp"
#include "MotionFilter.hpp"
#include "motionCompensate.hpp"
#include "cuda_runtime_api.h"
#include "cuda_runtime.h"

#include <unistd.h>
#include <time.h>

//int numStableObj = 0;
CStability* pStableObj = NULL;


VideoWriter g_WriteFile;

#define DEBUGTIME 	100
/*
*待添加
*1.稳像区域接口
*
*2.显示模式
*	0x00：画面显示模式（默认）
*	0x01：画中画演示模式
*	0x02：Side-By-Side演示模式
*
*3.防抖接口：
*	0x00： 关闭稳像校正
*	0x01： X/Y方向抖动校正
*	0x02： 旋转抖动校正
*	0x03： Zoom抖动校正
*/
extern "C" void tran_gray_cuda(unsigned char *src,unsigned char* dst,int nWidth,int nHeight);
extern "C" void Sobel_cuda(unsigned char *src, unsigned char *dst, int width, int height);

static unsigned int matime =0 ;
void Create_stable(void)
{
	#if 0
		if(numStableObj >= 4)
		{
			printf("overflow\n");
			return ;
		}
	#endif
	pStableObj = new CStability();
	if(NULL == pStableObj)
	{
		printf("error in new CStability\n");
		return ;
	}
	//numStableObj++;
	pStableObj->init();
	return;
}

void run_stable(Mat src,int nWidth,int nheight,uchar mode,unsigned int edge_h,unsigned int edge_v,affine_param* apout)
{
	// 1920   : edge_h  	320 pixel		edge_v	180	pixel	   // distance to the edge in pixel
	// 720 	: edge_h 		32   pixel		edge_v 	32	pixel
	pStableObj->RunStabilize(src,nWidth,nheight,mode,edge_h,edge_v,apout);
}

void destroy_stable(void)
{
	if(pStableObj != NULL)
		pStableObj->CloseStabilize(pStableObj->tss);
}	

void inline AnalysisMeResult( CStability * mcs)
{
    CStability * cs = mcs;
    stb_t* s = cs->tss;

    if (cs->MeErr)
    {
        if (cs->MeErr_cif == 0)
        {
            //用cif参数代替D1参数
            s->cur_af.cos = s->ap_cif.cos;
            s->cur_af.dx = s->ap_cif.dx * 2.0;
            s->cur_af.dy = s->ap_cif.dy * 2.0;
            s->cur_af.sin = s->ap_cif.sin;
            s->cur_af.theta = s->ap_cif.theta;
            s->cur_af.scale = s->ap_cif.scale;
            cs->MeErr = 0;
            s->stat[1] ++;
        }
        else if (cs->MeErr_qcif == 0)
        {
            //用qcif参数代替D1参数
            s->cur_af.cos = s->ap_qcif.cos;
            s->cur_af.dx = s->ap_qcif.dx * 4.0;
            s->cur_af.dy = s->ap_qcif.dy * 4.0;
            s->cur_af.sin = s->ap_qcif.sin;
            s->cur_af.theta = s->ap_qcif.theta;
            s->cur_af.scale = s->ap_qcif.scale;
            cs->MeErr = 0;
            s->stat[2] ++;
        }
        else if (cs->MeErr == 1)
        {
            cs->MeErr = 0;
            s->stat[3] ++;
        }
        else if (cs->MeErr_cif == 1)
        {
            //用cif参数代替D1参数
            s->cur_af.cos = s->ap_cif.cos;
            s->cur_af.dx = s->ap_cif.dx * 2.0;
            s->cur_af.dy = s->ap_cif.dy * 2.0;
            s->cur_af.sin = s->ap_cif.sin;
            s->cur_af.theta = s->ap_cif.theta;
            s->cur_af.scale = s->ap_cif.scale;
            s->stat[4] ++;
            cs->MeErr = 0;

        }
        else if (cs->MeErr_qcif == 1)
        {
            //用qcif参数代替D1参数
            s->cur_af.cos = s->ap_qcif.cos;
            s->cur_af.dx = s->ap_qcif.dx * 4.0;
            s->cur_af.dy = s->ap_qcif.dy * 4.0;
            s->cur_af.sin = s->ap_qcif.sin;
            s->cur_af.theta = s->ap_qcif.theta;
            s->cur_af.scale = s->ap_qcif.scale;
            s->stat[5] ++;
            cs->MeErr = 0;

        }
        else
        {
            cs->MeErr = 1;
            s->stat[7] ++;
        }
    }
    else
    {
        s->stat[0] ++;
    }

    if (!(cs->MeErr))
    {
        /*
        画面静止的时候，也可能有误差，
        导致画面静而不止，此处去掉微小误差，
        确保画面静止时不动*/
        if ((fabs(s->cur_af.dx) < 0.6)
                && (fabs(s->cur_af.dy) < 0.6)
                && (fabs(s->cur_af.theta) < 0.002))
        {
            s->cur_af.dx = 0.0;
            s->cur_af.dy = 0.0;
            s->cur_af.cos = 1.0;
            s->cur_af.sin = 0.0;
            s->cur_af.scale = 1.0;
            s->cur_af.theta = 0.0;
        }
    }
}
/********************************************/

CStability::CStability()
{
	m_modify = NULL;
	tss = NULL;
	pThis = this;
	
}

CStability::~CStability()
{
	
}

void CStability::init()
{
	int imgw = MAX_WIDTH;
	int imgh  = MAX_HEIGHT;

	if(tss == NULL)
		tss = new stb_t();
	stb_t* s = tss;
	allocspace();

	s->cur_af.cos = 1.0;
	s->cur_af.sin = 0.0;
	s->cur_af.dx = 0.0;
	s->cur_af.dy = 0.0;
	s->cur_af.theta = 0.0;
	s->cur_af.scale = 1.0;
	
	s->last_af.cos = 1.0;
	s->last_af.sin = 0.0;
	s->last_af.dx = 0.0;
	s->last_af.dy = 0.0;
	s->last_af.theta = 0.0;
	s->last_af.scale = 1.0;

	s->bak_af.cos = 1.0;
	s->bak_af.sin = 0.0;
	s->bak_af.dx = 0.0;
	s->bak_af.dy = 0.0;
	s->bak_af.theta = 0.0;
	s->bak_af.scale = 1.0;

	s->ap_cif.cos = 1.0;
	s->ap_cif.sin = 0.0;
	s->ap_cif.dx = 0.0;
	s->ap_cif.dy = 0.0;
	s->ap_cif.theta = 0.0;
	s->ap_cif.scale = 1.0;

	s->ap_qcif.cos = 1.0;
	s->ap_qcif.sin = 0.0;
	s->ap_qcif.dx = 0.0;
	s->ap_qcif.dy = 0.0;
	s->ap_qcif.theta = 0.0;
	s->ap_qcif.scale = 1.0;
	
	s->flt.xout = 0;
	s->flt.yout = 0;
	s->flt.sout = 1.0;
	s->flt.thetaout = 0.0;

	m_modify->cos = 1.0;
	m_modify->sin = 0.0;
	m_modify->dx = 0.0;
	m_modify->dy = 0.0;
	m_modify->theta = 0.0;
	m_modify->scale = 1.0;

	s->bReset = 0x01;

    	s->D1FpNum = 0;
    	s->CifFpNum = 0;
    	s->QcifFpNum = 0;

	s->i_width  = imgw;
	s->i_height = imgh;

	s->grid_w = (((imgw>>2)/22)>>2)<<2;
	s->grid_w = s->grid_w < 8 ? 8 : s->grid_w;
	s->grid_h = (((imgh>>2)/18)>>2)<<2;
	s->grid_h = s->grid_h < 8 ? 8 : s->grid_h;

	//s->fD1Cur->i_stride[0] = imgw;
	//s->fD1Cur->i_stride[1] = imgw>>1;
	//s->fD1Cur->i_stride[2] = imgw>>2;

	//s->fD1Out->i_stride[0] = imgw;
	//s->fD1Out->i_stride[1] = imgw>>1;
	//s->fD1Out->i_stride[2] = imgw>>2;

	memset(s->D1Fp,0,   (imgw>>2)*(imgh>>2)*sizeof(FPOINT));
	memset(s->CifFp,0,  (imgw>>2)*(imgh>>2)*sizeof(FPOINT));
	memset(s->QcifFp,0, (imgw>>2)*(imgh>>2)*sizeof(FPOINT));

	OpenStabilize(s);


	g_WriteFile.open("test.avi", CV_FOURCC('M','J','P','G'), 25.0, cvSize(720, 576));
			
}

void CStability::OpenStabilize(stb_t* s)
{
    //int i = 0;
    //int MP = 4;  //Kalman: number of measure vector dimensions
    //int DP = 8;  //Kalman: number of state   vector dimensions
    //int CP = 0;  //Kalman: number of control vector dimensions
   s->g_pKalman = kkalman.init();
    /*Kalman Filter Init*/
   //s->g_pKalman = kkalman.KalmanOpen(DP, MP, CP);
    if (s->g_pKalman == NULL)
    {
        return ;
    }

    kkalman.KalmanInitParam(s->g_pKalman, 0.0, 0.0, 0.0, 1.0, 0.0);	
    return ;
}

void CStability::allocspace()
{
	int imgw = MAX_WIDTH;
	int imgh  = MAX_HEIGHT;
	unsigned int datablock = imgw*imgh*1;

	stb_t* ms = tss;
	
	#if 0

	cudaError_t cudaStatus;

	ms->fD1Cur = new stb_frame_t();
	ms->fD1Out = new stb_frame_t();

	cudaStatus = cudaMalloc((void**)&ms->fD1Cur->buffer[0], datablock);
	cudaStatus = cudaMalloc((void**)&ms->fD1Out->buffer[0], datablock);

	cudaStatus = cudaMalloc((void**)&ms->fD1Ref, datablock);
	cudaStatus = cudaMalloc((void**)&ms->fD1CurSobel, datablock);
	cudaStatus = cudaMalloc((void**)&ms->fD1RefSobel, datablock);
	
	cudaStatus = cudaMalloc((void**)&ms->fCifCur, datablock>>2);
	cudaStatus = cudaMalloc((void**)&ms->fCifRef, datablock>>2);
	cudaStatus = cudaMalloc((void**)&ms->fCifCurSobel, datablock>>2);
	cudaStatus = cudaMalloc((void**)&ms->fCifRefSobel, datablock>>2);

	cudaStatus = cudaMalloc((void**)&ms->fQcifCur, datablock>>4);
	cudaStatus = cudaMalloc((void**)&ms->fQcifRef, datablock>>4);
	cudaStatus = cudaMalloc((void**)&ms->fQcifCurSobel, datablock>>4);
	cudaStatus = cudaMalloc((void**)&ms->fQcifRefSobel, datablock>>4);
	

	ms->D1Fp = new FPOINT[(MAX_WIDTH>>2)*(MAX_HEIGHT>>2)];
	ms->CifFp = new FPOINT[(MAX_WIDTH>>2)*(MAX_HEIGHT>>2)];
	ms->QcifFp = new FPOINT[(MAX_WIDTH>>2)*(MAX_HEIGHT>>2)];

	m_modify = new affine_param();
	
	#else

	ms->fD1Cur = new stb_frame_t();
	ms->fD1Out = new stb_frame_t();
	
	ms->fD1Cur->buffer[0] = new unsigned char[datablock];
	ms->fD1Out->buffer[0] = new unsigned char[datablock];
	
	ms->fD1Ref = new unsigned char[datablock];
	ms->fD1CurSobel = new unsigned char[datablock];
	ms->fD1RefSobel = new unsigned char[datablock];
	
	ms->fCifCur = new unsigned char[datablock>>2];
	ms->fCifRef = new unsigned char[datablock>>2];
	ms->fCifCurSobel = new unsigned char[datablock>>2];
	ms->fCifRefSobel = new unsigned char[datablock>>2];
	
	ms->fQcifCur = new unsigned char[datablock>>4];
	ms->fQcifRef = new unsigned char[datablock>>4];
	ms->fQcifCurSobel = new unsigned char[datablock>>4];
	ms->fQcifRefSobel = new unsigned char[datablock>>4];

	ms->D1Fp = new FPOINT[(MAX_WIDTH>>2)*(MAX_HEIGHT>>2)];
	ms->CifFp = new FPOINT[(MAX_WIDTH>>2)*(MAX_HEIGHT>>2)];
	ms->QcifFp = new FPOINT[(MAX_WIDTH>>2)*(MAX_HEIGHT>>2)];

	m_modify = new affine_param();
	#endif
	
	return ;
}

void CStability::destroy()
{
	printf("destroy\n");
}

void CStability::run()
{
	printf("run\n");
}

void CStability::CloseStabilize(stb_t *s)
{
	s->bReset = 1;
	//kkalman.KalmanClose(s->g_pKalman);
	//CloseStabilize(s);
}

void CStability::showPoints(unsigned char code)
{
	string pNum;
	unsigned int num[3] = {0};
	Point keypoints;
	stb_t* s = tss;
	h_mb_num = (s->i_width >> 2)/s->grid_w;
	v_mb_num = (s->i_height>>2)/s->grid_h;

	int height,width;
	unsigned char* dst1 = (unsigned char*)malloc(s->i_height*s->i_width);
	unsigned char* dst2 = (unsigned char*)malloc(s->i_height*s->i_width);
	unsigned char *src1,*src2;
	switch(code)
	{
		case 0 :
			src1 = s->fQcifCur;
			src2 = s->fQcifRef;
			height = s->i_height>>2;
			width  = s->i_width>>2;
			break;
		case 1:
			src1 = s->fCifCur;
			src2 = s->fCifRef;
			height = s->i_height>>1;
			width  = s->i_width>>1;
			break;
		case 2:
			src1 = s->fD1Cur->buffer[0];
			src2 = s->fD1Ref;
			height = s->i_height;
			width  = s->i_width;
			break;
		default :
			break;
	}

	Mat mat1 = Mat(height,width,CV_8UC1,dst1);
	memcpy(mat1.data,src1,height*width);
	Mat mat2 = Mat(height,width,CV_8UC1,dst2);
	memcpy(mat2.data,src2,height*width);
	num[0] = 0;
	for(int n = 0;n<h_mb_num*v_mb_num;n++)
	{	
		if(s->D1Fp[n].ftv)
		{	
			if(code == 0)
			{
				keypoints.x = s->QcifFp[n].x;
				keypoints.y = s->QcifFp[n].y;
				//printf("!!!cur  keypoints.x = %d \n",keypoints.x);
				//printf("!!!cur  keypoints.y = %d \n\n",keypoints.y);
				pNum = format("%d",num[0]);
				circle(mat1,keypoints,1,Scalar(255,255,255),0.5,8);
				//putText(mat1,pNum , Point(keypoints.x,keypoints.y), 0, 0.4, Scalar::all(255), 0.8, 8);
			
				keypoints.x = s->QcifFp[n].x+s->QcifFp[n].dx;
				keypoints.y = s->QcifFp[n].y+s->QcifFp[n].dy;
				//printf("!!!ref  QcifFp.dx = %d \n",s->QcifFp[n].dx);
				//printf("!!!ref  QcifFp.dy = %d \n\n",s->QcifFp[n].dy);
				//printf("!!!ref  keypoints.x = %d \n",keypoints.x);
				//printf("!!!ref  keypoints.y = %d \n\n",keypoints.y);
				circle(mat2,keypoints,1,Scalar(255,255,255),0.5,8);
				//putText(mat2,pNum , Point(keypoints.x,keypoints.y), 0, 0.4, Scalar::all(255), 0.8, 8);
				//printf("**************showPoints end*****************\n\n");		
			}
			else if(code == 1)
			{
				keypoints.x = s->CifFp[n].x;
				keypoints.y = s->CifFp[n].y;
				pNum = format("%d",n);
				circle(mat1,keypoints,1,Scalar(255,255,255),1,8);
				//putText(mat1, pNum, Point(keypoints.x,keypoints.y), 0, 0.4, Scalar::all(255), 1, 8);
				/******ref******/
				keypoints.x = s->CifFp[n].x+s->CifFp[n].dx;
				keypoints.y = s->CifFp[n].y+s->CifFp[n].dy;
				circle(mat2,keypoints,1,Scalar(255,255,255),1,8);
				//putText(mat2,pNum , Point(keypoints.x,keypoints.y), 0, 0.4, Scalar::all(255), 0.8, 8);
				/******end******/	
			}	
			else if(code == 2)
			{			
				keypoints.x = s->D1Fp[n].x;
				keypoints.y = s->D1Fp[n].y;
				pNum = format("%d",n);
				circle(mat1,keypoints,2,Scalar(255,255,255),1,8);
				//putText(mat1, pNum, Point(keypoints.x,keypoints.y), 0, 0.4, Scalar::all(255), 0.8, 8);
			
				/******ref******/
				keypoints.x = s->D1Fp[n].x+s->D1Fp[n].dx;
				keypoints.y = s->D1Fp[n].y+s->D1Fp[n].dy;
				circle(mat2,keypoints,2,Scalar(255,255,255),1,8);
				//putText(mat2,pNum , Point(keypoints.x,keypoints.y), 0, 0.4, Scalar::all(255), 0.8, 8);
				/******end******/		
			}	
			num[0]++;
		}
	}

	printf("s->CifFpNum = %d\n",s->CifFpNum);
	namedWindow("cur",CV_WINDOW_NORMAL);
	imshow("cur",mat1);
	namedWindow("ref",CV_WINDOW_NORMAL);
	imshow("ref",mat2);

	return ;
}

void extractUYVY2Gray(Mat src, Mat dst)
{
	int ImgHeight, ImgWidth,ImgStride;

	ImgWidth = src.cols;
	ImgHeight = src.rows;
	ImgStride = ImgWidth*2;
	uint8_t  *  pDst8_t;
	uint8_t *  pSrc8_t;

	pSrc8_t = (uint8_t*)(src.data);
	pDst8_t = (uint8_t*)(dst.data);
//#pragma UNROLL 4
	#pragma omp parallel for
	for(int y = 0; y < ImgHeight*ImgWidth; y++)
	{
		pDst8_t[y] = pSrc8_t[y*2+1];
	}
}



int CStability::RunStabilize(Mat src,int nWidth, int nHeight,uchar mode,unsigned int cedge_h,unsigned int cedge_v,affine_param* apout)
{
	//unsigned int tt0 = OSA_getCurTimeInMsec();
	int i;
	CStability *cs = pThis;
	stb_t* s = cs->tss;
	affine_param *ap_modify = cs->m_modify;
	static char pp = 0;
	int ttt,nnn = 0;
	cudaEvent_t	start, stop;
	static float pauseInput = 1000;
	//time12[0] = OSA_getCurTimeInMsec();
	/*   create the Mat obj again and again for attach to the new address */
	//mfout = Mat(s->i_height, s->i_width, CV_8UC1,s->fD1Out->buffer[0]);

	if(nWidth != 1920)
	{
		s->i_width = nWidth;
		s->i_height = nHeight;

		s->grid_w = (((nWidth>>2)/22)>>2)<<2;
		s->grid_w = s->grid_w < 8 ? 8 : s->grid_w;
		s->grid_h = (((nHeight>>2)/18)>>2)<<2;
		s->grid_h = s->grid_h < 8 ? 8 : s->grid_h;
	}
	
	edge_h = cedge_h / s->grid_h ;
	edge_v = cedge_v / s->grid_w ;
		
	//for (int i = 0; i < (MAX_WIDTH>> 2) * (MAX_HEIGHT>> 2); i++)
    	//{
       //	s->D1Fp[i].ftv = 0x00;
       //	s->CifFp[i].ftv = 0x00;
       // 	s->QcifFp[i].ftv = 0x00;
    	//}
	memset(s->D1Fp,0,   (MAX_WIDTH>>2)*(MAX_HEIGHT>>2)*sizeof(FPOINT));
	memset(s->CifFp,0,  (MAX_WIDTH>>2)*(MAX_HEIGHT>>2)*sizeof(FPOINT));
	memset(s->QcifFp,0, (MAX_WIDTH>>2)*(MAX_HEIGHT>>2)*sizeof(FPOINT));

	//unsigned int tt1 = OSA_getCurTimeInMsec();
	//printf("tt1 - tt0 = %u\n",tt1 - tt0);
	// record the video
	#if 0
	Mat temp = Mat(576,720,CV_8UC3);

	cudaMemcpy(temp.data, src.data, nWidth*nHeight*3, cudaMemcpyDeviceToHost);
	
	if(1 /*&& isfilt*/)
	{
		g_WriteFile << temp;
		imshow("vvvvvv",temp);
		waitKey(1);
	}
	return 0;
	#endif

	

#if 0
	cudaError_t cudaStatus;

	unsigned char *dfcur = src.data;
	unsigned char *src_gray = NULL;
	unsigned char *dfcif = NULL;
	unsigned char *dfcur_sobel = NULL;
	unsigned char *dfcif_sobel= NULL;
	unsigned char *dfqcif = NULL;
	unsigned char *dfqcif_sobel= NULL;
	cudaStatus = cudaMalloc((void**)&src_gray,sizeof(unsigned char)*nWidth*nHeight);

	tran_gray_cuda(src.data,src_gray, nWidth,nHeight);

	cudaStatus = cudaMalloc((void**)&dfcur_sobel,s->i_height*s->i_width);
	cudaStatus = cudaMalloc((void**)&dfcif, s->i_height>>1*s->i_width>>1);
	cudaStatus = cudaMalloc((void**)&dfcif_sobel, s->i_height>>1*s->i_width>>1);
	cudaStatus = cudaMalloc((void**)&dfqcif, s->i_height>>2*s->i_width>>2);
	cudaStatus = cudaMalloc((void**)&dfqcif_sobel, s->i_height>>2*s->i_width>>2);



printf("33333333333333333333\n");
return 0;	
	Sobel_cuda(src.data, dfcur_sobel, s->i_width, s->i_height);
	Sobel_cuda(dfcif, dfcif_sobel, s->i_width>>1, s->i_height>>1);
	Sobel_cuda(dfqcif, dfcif_sobel, s->i_width>>2, s->i_height>>2);
	
printf("444444444444444444444\n");

	
	
#else

	mfout = Mat(s->i_height, s->i_width, CV_8UC1);

	mfcur = Mat(s->i_height, s->i_width, CV_8UC1, s->fD1Cur->buffer[0]);
	mfCifCur = Mat(s->i_height>>1, s->i_width>>1, CV_8UC1,s->fCifCur);
	mfQcifCur = Mat(s->i_height>>2, s->i_width>>2, CV_8UC1,s->fQcifCur);

	mfcur_sobel = Mat(s->i_height, s->i_width, CV_8UC1,s->fD1CurSobel);
	mfCifCur_sobel = Mat(s->i_height>>1,s->i_width>>1, CV_8UC1,s->fCifCurSobel);
	mfQCifCur_sobel = Mat(s->i_height>>2,s->i_width>>2, CV_8UC1,s->fQcifCurSobel);	

	mfCur_ref = Mat(s->i_height, s->i_width, CV_8UC1,s->fD1Ref);
	mfCifCur_ref = Mat(s->i_height>>1, s->i_width>>1, CV_8UC1,s->fCifRef);
	mfQCifCur_ref = Mat(s->i_height>>2,s->i_width>>2, CV_8UC1,s->fQcifRef);

	mfcur_sobel_ref = Mat(s->i_height, s->i_width, CV_8UC1,s->fD1RefSobel);
	mfCifCur_sobel_ref = Mat(s->i_height>>1, s->i_width>>1, CV_8UC1,s->fCifRefSobel);
	mfQCifCur_sobel_ref = Mat(s->i_height>>2, s->i_width>>2, CV_8UC1,s->fQcifRefSobel);
	
	//src.copyTo(mfcur);
	//extractUYVY2Gray(tmp,mfcur);
	//memcpy(mfcur.data,src.data,nWidth*nHeight);
	//time12[2] = OSA_getCurTimeInMsec();
	#if 0
	Mat tmp = Mat(nHeight,nWidth,CV_8UC3);
	Mat tmpget = Mat(nHeight,nWidth,CV_8UC3);
	cudaMemcpy(tmp.data, src.data, nWidth*nHeight*3, cudaMemcpyDeviceToHost);
	cvtColor(tmp,mfcur,CV_BGR2GRAY);
	#endif
	
	memcpy(mfcur.data,src.data,nWidth*nHeight);	
	//time12[3] = OSA_getCurTimeInMsec();

	/*	pre-process	*/
	//unsigned int tt0 = OSA_getCurTimeInMsec();
	preprocess(mfcur, mfCifCur, mfQcifCur,mfcur_sobel,mfCifCur_sobel,mfQCifCur_sobel,s);
	
	// 8 - 16 ms ------  11ms
#endif

	//time12[4] = OSA_getCurTimeInMsec();
   	//edge_v = 5;
   	//edge_h = 3;
	/*  init the param of kalman	 */


	if(s->bReset)
	{
		time12[0] = 0;
		time12[3] = 0;
		time12[4] = 0;
		s->bReset = 0;
		InitFilter(&(s->last_af),ap_modify,&(s->flt));
		kkalman.KalmanInitParam(s->g_pKalman, 0.0, 0.0, 0.0, 1.0, 0.0);
		for (i = 0; i < 8; i++)
		{
			s->stat[i] = 0;
		}
		MeErr_cif = 10;
		MeErr_qcif = 10;
		MeErr = 10;
		//	这里用Mat.data交换指针，那么之前申请的空间无意义了，先这样调试，待查
	    FRAME_EXCHANGE(s->fD1Ref,		s->fD1Cur->buffer[0]);	
	    FRAME_EXCHANGE(s->fD1RefSobel, 	s->fD1CurSobel);
	    FRAME_EXCHANGE(s->fCifRef,      	s->fCifCur);
	    FRAME_EXCHANGE(s->fCifRefSobel,	s->fCifCurSobel);
	    FRAME_EXCHANGE(s->fQcifRef,    	s->fQcifCur);
	    FRAME_EXCHANGE(s->fQcifRefSobel, 	s->fQcifCurSobel);

	}
	 else
	 {

	 	//time12[5] = OSA_getCurTimeInMsec();
		
		 /*  	Find the feature points	*/

		//unsigned int tt0 = OSA_getCurTimeInMsec();

		kfindFtp.findFtp(&(s->edgeTh), s->QcifFp, &(s->QcifFpNum),s->i_width >> 2, s->i_height >> 2, 
				s->i_width >> 1, s->i_width,s->grid_w, s->grid_h, edge_h, edge_v,
				s->CifFp, &(s->CifFpNum),s->D1Fp, &(s->D1FpNum),
				mfcur_sobel,mfCifCur_sobel,mfQCifCur_sobel);
		//unsigned int tt1 = OSA_getCurTimeInMsec();
		//printf("tt1 - tt0 = %u\n",tt1 - tt0);
		//   1 ms
		//time12[6] = OSA_getCurTimeInMsec();

		/*		此处待添加调试  FTP_SHOW		*/

		//unsigned int tt0 = OSA_getCurTimeInMsec();
		if (s->QcifFpNum < 3) //特征点太少
		{
			//time12[7] = 0;
			//time12[8] = 0;
			
			MeErr_qcif = 10;
			MeErr_cif = 10;
			MeErr = 10;
		}
		else		
		{
		//	time12[7] = OSA_getCurTimeInMsec();
			unsigned int ts1 = OSA_getCurTimeInMsec();
			RunMatchingPoint(s,&MeErr,&MeErr_cif,&MeErr_qcif);	
			unsigned int ts2 = OSA_getCurTimeInMsec();
 		//	time12[8] = OSA_getCurTimeInMsec();
 			printf("matching time = %d\n",ts2 - ts1);
			matime++;
		}

		//unsigned int tt1 = OSA_getCurTimeInMsec();
		//printf("tt1 - tt0 = %u\n",tt1 - tt0);
		
		//time12[9] = OSA_getCurTimeInMsec();
		AnalysisMeResult(cs);
		// time12[10] = OSA_getCurTimeInMsec();
		MotionFilter(cs);
		// time12[11] = OSA_getCurTimeInMsec();
		 
		memcpy(apout,cs->m_modify,sizeof(affine_param));
			
		//time12[12] = OSA_getCurTimeInMsec();
		//MotionProcess(cs,tmp,tmpget,mode);
		//MotionProcess(cs,src,dst,mode);
		//cudaMemcpy(dst.data, tmpget.data, s->i_height*s->i_width*3, cudaMemcpyHostToDevice);
		#if 0
		#if 1
		float elapsedTime;
		( (		cudaEventCreate	(	&start)	) );
		( (		cudaEventCreate	(	&stop)	) );
		( (		cudaEventRecord	(	start,	0)	) );
		#endif
		MotionProcess(cs, src,dst,mode);
		#if 1
		((	cudaEventRecord(	stop,	0	)));
		((	cudaEventSynchronize(	stop)	));
		((	cudaEventSynchronize(	start)	));
		(	cudaEventElapsedTime(	&elapsedTime,	start,	stop));	
		printf("Time :	%3.1f	ms \n", elapsedTime);
		((	cudaEventDestroy(	start	)));
		((	cudaEventDestroy(	stop	)));		
		#endif
		time12[13] = OSA_getCurTimeInMsec();
		#if 0
			Mat mattmp = Mat(576,720,CV_8UC3);
			cudaMemcpy(mattmp.data, dst.data, s->i_height*s->i_width*3, cudaMemcpyDeviceToHost);
			printf("imshow\n");
			imshow("mfoutmfout",mattmp);
			waitKey(0);
		#endif
		#endif
	    FRAME_EXCHANGE(s->fD1Ref,s->fD1Cur->buffer[0]);	
	    FRAME_EXCHANGE(s->fD1RefSobel, mfcur_sobel.data);
	    FRAME_EXCHANGE(s->fCifRef,      mfCifCur.data);
	    FRAME_EXCHANGE(s->fCifRefSobel, mfCifCur_sobel.data);
	    FRAME_EXCHANGE(s->fQcifRef,     mfQcifCur.data);
	    FRAME_EXCHANGE(s->fQcifRefSobel, mfQCifCur_sobel.data);	
	    //time12[14] = OSA_getCurTimeInMsec();

#if 0
	printf("cs->m_modify.dx = %f\n",cs->m_modify->dx);
	printf("cs->m_modify.dy = %f\n",cs->m_modify->dy);
	printf("cs->m_modify.cos = %f\n",cs->m_modify->cos);
	printf("cs->m_modify.sin = %f\n",cs->m_modify->sin);

	printf("get src time : %u\n",time12[3] - time12[2]);
	printf("preprocess time : %u\n",time12[4] - time12[3]);	
  	printf("match time : %u\n",time12[8] - time12[7]);	
	//printf("motioncpmpensate time : %u\n",time12[13] - time12[12]);	
	printf("ellllllllllllllll time : %u\n",time12[14] - time12[0]);
#endif
	 }
	
	//analytime();
	//unsigned int tt1 = OSA_getCurTimeInMsec();
	//printf("tt1 - tt0 = %u\n",tt1 - tt0);
	return 0;
}

void CStability::analytime()
{
	int i;
	unsigned int tmp;
	if(anytimenum == DEBUGTIME -1)
	{
		for(i = 0;i<14;i++)
			avr[i] = anytime[i]/DEBUGTIME;
		avr[7] = anytime[7]/matime;
		avr[15] = anytime[15]/DEBUGTIME;

		anytimenum = 0;

		//for(i = 0;i<11;i++)
		{
			i = 2;
			printf("\nget src min[%d] = %u\n",i,mintime[i]);
			printf("get src avr[%d] = %u\n",i,avr[i]);	
			printf("get src max[%d] = %u\n",i,maxtime[i]);
			
			i = 3;
			printf("\npreprocess min[%d] = %u\n",i,mintime[i]);
			printf("preprocess avr[%d] = %u\n",i,avr[i]);	
			printf("preprocess max[%d] = %u\n",i,maxtime[i]);

			//i = 5;
			//printf("findpoint min[%d] = %u\n",i,mintime[i]);
			//printf("findpoint avr[%d] = %u\n",i,avr[i]);	
			//printf("findpoint max[%d] = %u\n",i,maxtime[i]);	

			i = 7;
			printf("match min[%d] = %u\n",i,mintime[i]);
			printf("match avr[%d] = %u\n",i,avr[i]);	
			printf("match max[%d] = %u\n",i,maxtime[i]);	

			//i = 9;
			//printf("analy min[%d] = %u\n",i,mintime[i]);
			//printf("analy avr[%d] = %u\n",i,avr[i]);	
			//printf("analy max[%d] = %u\n",i,maxtime[i]);	
			
			i = 10;
			printf("motionfilter min[%d] = %u\n",i,mintime[i]);
			printf("motionfilter avr[%d] = %u\n",i,avr[i]);	
			printf("motionfilter max[%d] = %u\n",i,maxtime[i]);	

			//i = 12;
			//printf("motioncpmpensate min[%d] = %u\n",i,mintime[i]);
			//printf("motioncpmpensate avr[%d] = %u\n",i,avr[i]);	
			//printf("motioncpmpensate max[%d] = %u\n",i,maxtime[i]);	

			i=13;
			printf("allcost min[%d] = %u\n",i,mintime[i]);
			printf("allcost avr[%d] = %u\n",i,avr[i]);	
			printf("allcost max[%d] = %u\n",i,maxtime[i]);				
		}
			
				
	}
	else if(anytimenum > DEBUGTIME -1)
	{
		anytimenum = 0;	
	}
	
	if(anytimenum == 0)
	{
		memset(anytime,0,20*sizeof(unsigned long));
		memset(mintime,10000,20*sizeof(unsigned int));
		memset(maxtime,0,20*sizeof(unsigned int));
		matime = 0;
	}
		
	for(i = 0;i<=12;i++)
	{
		tmp = (time12[i+1] - time12[i]);
		anytime[i] += tmp;
		
		if(tmp < mintime[i])
			mintime[i] = tmp;		
		else if(tmp>maxtime[i])
			maxtime[i] = tmp;
	}	
		tmp = (time12[14] - time12[0]);
		anytime[13] += tmp;
		
		if(tmp < mintime[13])
			mintime[13] = tmp;		
		else if(tmp>maxtime[13])
			maxtime[13] = tmp;

	anytimenum ++ ;

	return ;	
}

