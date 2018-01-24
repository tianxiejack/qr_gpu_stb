#include "stable.hpp"
#include "stdio.h"
#include <string.h>
#include "matchingPoint.hpp"
#include "preprocess.hpp"
#include "MotionFilter.hpp"
#include "motionCompensate.hpp"

//int numStableObj = 0;
CStability* pStableObj = NULL;

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

void run_stable(Mat src,Mat dst,int nWidth,int nheight,uchar mode,unsigned int edge_h,unsigned int edge_v)
{
	// 1920   : edge_h  	320 pixel		edge_v	180	pixel	   // distance to the edge in pixel
	// 720 	: edge_h 		32   pixel		edge_v 	32	pixel
	pStableObj->RunStabilize(src,dst,nWidth,nheight,mode,edge_h,edge_v);
}

void destroy_stable(void)
{
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

	s->fD1Cur->i_stride[0] = imgw;
	s->fD1Cur->i_stride[1] = imgw>>1;
	s->fD1Cur->i_stride[2] = imgw>>2;

	s->fD1Out->i_stride[0] = imgw;
	s->fD1Out->i_stride[1] = imgw>>1;
	s->fD1Out->i_stride[2] = imgw>>2;

	memset(s->D1Fp,0,   (imgw>>2)*(imgh>>2)*sizeof(FPOINT));
	memset(s->CifFp,0,  (imgw>>2)*(imgh>>2)*sizeof(FPOINT));
	memset(s->QcifFp,0, (imgw>>2)*(imgh>>2)*sizeof(FPOINT));

	OpenStabilize(s);
}

void CStability::OpenStabilize(stb_t* s)
{
    int i = 0;
    int MP = 4;  //Kalman: number of measure vector dimensions
    int DP = 8;  //Kalman: number of state   vector dimensions
    int CP = 0;  //Kalman: number of control vector dimensions

    /*Kalman Filter Init*/
    s->g_pKalman = kkalman.KalmanOpen(DP, MP, CP);
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
	kkalman.KalmanClose(s->g_pKalman);
	CloseStabilize(s);
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


int CStability::RunStabilize(Mat src,Mat dst,int nWidth, int nHeight,uchar mode,unsigned int cedge_h,unsigned int cedge_v)
{
	int i;
	CStability *cs = pThis;
	stb_t* s = cs->tss;
	affine_param *ap_modify = cs->m_modify;
	static char pp = 0;
	int ttt,nnn = 0;
	/*   create the Mat obj again and again for attach to the new address */
	//mfout = Mat(s->i_height, s->i_width, CV_8UC1,s->fD1Out->buffer[0]);
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

	//dst = Mat(s->i_height, s->i_width, CV_8UC1,s->fD1Out->buffer[0]);
	src.copyTo(mfcur);

	if(nWidth == 1920)
	{
		edge_h = cedge_h / s->grid_w ;
		edge_v = cedge_v / s->grid_h ;
	}
	else
	{
		s->grid_w = (((nWidth>>2)/22)>>2)<<2;
		s->grid_w = s->grid_w < 8 ? 8 : s->grid_w;
		s->grid_h = (((nHeight>>2)/18)>>2)<<2;
		s->grid_h = s->grid_h < 8 ? 8 : s->grid_h;

		edge_h = cedge_h / s->grid_h ;
		edge_v = cedge_v / s->grid_w ;
	}
	
	
	for (int i = 0; i < (s->i_width>> 2) * (s->i_height>> 2); i++)
    	{
       	s->D1Fp[i].ftv = 0x00;
       	s->CifFp[i].ftv = 0x00;
        	s->QcifFp[i].ftv = 0x00;
    	}

	/*	pre-process	*/
	preprocess(mfcur, mfCifCur, mfQcifCur,mfcur_sobel,mfCifCur_sobel,mfQCifCur_sobel,s);

	/*  init the param of kalman	 */
	if(s->bReset)
	{
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
		 /*  	Find the feature points	*/
		kfindFtp.findFtp(&(s->edgeTh), s->QcifFp, &(s->QcifFpNum),s->i_width >> 2, s->i_height >> 2, 
				s->i_width >> 1, s->i_width,s->grid_w, s->grid_h, edge_h, edge_v,
				s->CifFp, &(s->CifFpNum),s->D1Fp, &(s->D1FpNum),
				mfcur_sobel,mfCifCur_sobel,mfQCifCur_sobel);

		/*		此处待添加调试  FTP_SHOW		*/
		if (s->QcifFpNum < 3) //特征点太少
		{
			MeErr_qcif = 10;
			MeErr_cif = 10;
			MeErr = 10;
		}
		else		
			RunMatchingPoint(s,&MeErr,&MeErr_cif,&MeErr_qcif);		

		AnalysisMeResult(cs);
		
		MotionFilter(cs);

		MotionProcess(cs,mfcur,dst,mode);
		
	    FRAME_EXCHANGE(s->fD1Ref,s->fD1Cur->buffer[0]);	
	    FRAME_EXCHANGE(s->fD1RefSobel, mfcur_sobel.data);
	    FRAME_EXCHANGE(s->fCifRef,      mfCifCur.data);
	    FRAME_EXCHANGE(s->fCifRefSobel, mfCifCur_sobel.data);
	    FRAME_EXCHANGE(s->fQcifRef,     mfQcifCur.data);
	    FRAME_EXCHANGE(s->fQcifRefSobel, mfQCifCur_sobel.data);		 
	 }
	return 0;
}
