
#include "matchingPoint.hpp"
#include "math.h"
#include "stdio.h"

//#include <emmintrin.h>
//#include <xmmintrin.h>
//#include <pmmintrin.h>
//#include <mmintrin.h>
#include <omp.h>
//#include <smmintrin.h>
#include <arm_neon.h>
#include "sse2neon.hpp"
#include "osa.h"

void IMG_mad_8x8(const unsigned char * refImg,const unsigned char * srcImg,int pitch,int sx, int sy,unsigned int * match)                                                              
{                                                              
   int i, j, x, y, matx, maty;                                
   unsigned int matpos, matval;                                   
   unsigned char *psrc,*pref;														  
   unsigned int acc ;
   unsigned short _acc[8];
   matval = ~0U;                                              
   matx = maty = 0;                           
   
	 for (y = 0; y < sy; y++)
	 {
	   for (x = 0; x < sx; x++)                                   
	   {
	   	  acc = 0;    
		  
		 #if 0
			  for(j = 0;j<8;j++)
			  {
			  	psrc = (unsigned char *)(srcImg + j*pitch);
				pref = (unsigned char *)(refImg + (j+y)*pitch + x);
				for(i =0;i<8;i++)
				{
					acc += abs(psrc[i] - pref[i]);
				}
			  }	
		#else
				#pragma UNROLL(8)
				for(j = 0;j<8;j++)
				{
					int srcPtr = j*pitch;
					int refPtr = (j+y)*pitch + x;

					#if 1
					__m128i _a128,_b128,rslt;

					_a128 = _mm_loadl_epi64((__m128i*)(srcImg + srcPtr));
					_b128 = _mm_loadl_epi64((__m128i*)(refImg + refPtr));
					rslt = _mm_sad_epu8(_a128,_b128);
					//_acc[j] = _mm_extract_epi32(rslt,0)
					// extern int _mm_extract_epi32(_m128i src,const int ndx);
					_acc[j] = vgetq_lane_u16((uint16x8_t)rslt,0);
					//_mm_empty();
					#else

					uint8x8_t _a,_b,rslt;

					_a = vld1_u8((const uint8_t*)(srcImg + srcPtr));
					_b = vld1_u8((const uint8_t*)(refImg + refPtr));
					rslt = vabd_u8(_a,_b);
					//uint8_t vget_lane_u8 (uint8x8_t __a, const int __b)\u037e
		
					_acc[j] += vget_lane_u8(rslt,0);
					_acc[j] += vget_lane_u8(rslt,1);
					_acc[j] += vget_lane_u8(rslt,2);
					_acc[j] += vget_lane_u8(rslt,3);
					_acc[j] += vget_lane_u8(rslt,4);
					_acc[j] += vget_lane_u8(rslt,5);
					_acc[j] += vget_lane_u8(rslt,6);
					_acc[j] += vget_lane_u8(rslt,7);				
					
					#endif
				}
				acc = (_acc[0] + _acc[1] + _acc[2] + _acc[3] + _acc[4] + _acc[5] + _acc[6] + _acc[7]);
		#endif

		 if (acc < matval)                                  
	 	 {           
		    matval = acc;                                  
		    matx   = x;                                    
		    maty   = y;
	   	 }    
	   }
	 }
	 															  
   matpos    = (0xffff0000 & (matx << 16)) | (0x0000ffff & maty);                           
   match[0] = matpos;                                         
   match[1] = matval;                                         
}                                                              


int brinv_1D(float a[], int n)
{
    int i, j, k, l, u, v;
    float d, p;
    int is[4];
    int js[4];
    //is=(int*)malloc(n*sizeof(int));
    //js=(int*)malloc(n*sizeof(int));
    for (k = 0; k <= n - 1; k++)
    {
        d = 0.0;
        for (i = k; i <= n - 1; i++)
        {
            for (j = k; j <= n - 1; j++)
            {
                l = i * n + j;
                p = fabs(a[l]);
                if (p > d)
                {
                    d = p;
                    is[k] = i;
                    js[k] = j;
                }
            }
        }
        if (d + 1.0 == 1.0)
        {
//          free (is);
//          free (js);
            return 0;
        }
        if (is[k] != k)
        {
            for (j = 0; j <= n - 1; j++)
            {
                u = k * n + j;
                v = is[k] * n + j;
                p = a[u];
                a[u] = a[v];
                a[v] = p;
            }
        }
        if (js[k] != k)
        {
            for (i = 0; i <= n - 1; i++)
            {
                u = i * n + k;
                v = i * n + js[k];
                p = a[u];
                a[u] = a[v];
                a[v] = p;
            }
        }
        l = k * n + k;
        a[l] = 1.0 / a[l];
        for (j = 0; j <= n - 1; j++)
        {
            if (j != k)
            {
                u = k * n + j;
                a[u] = a[u] * a[l];
            }
        }
        for (i = 0; i <= n - 1; i++)
        {
            if (i != k)
            {
                for (j = 0; j <= n - 1; j++)
                {
                    if (j != k)
                    {
                        u = i * n + j;
                        a[u] = a[u] - a[i * n + k] * a[k * n + j];
                    }
                }
            }
        }
        for (i = 0; i <= n - 1; i++)
        {
            if (i != k)
            {
                u = i * n + k;
                a[u] = -a[u] * a[l];
            }
        }
    }
    for (k = n - 1; k >= 0; k--)
    {
        if (js[k] != k)
        {
            for (j = 0; j <= n - 1; j++)
            {
                u = k * n + j;
                v = js[k] * n + j;
                p = a[u];
                a[u] = a[v];
                a[v] = p;
            }
        }
        if (is[k] != k)
        {
            for (i = 0; i <= n - 1; i++)
            {
                u = i * n + k;
                v = i * n + is[k];
                p = a[u];
                a[u] = a[v];
                a[v] = p;
            }
        }
    }
//  free( is);
//  free( js);
    return -1;
}


static void disove_AtA(float *pA, float *pQ, int m, int n)
{
    int x, y, i;
    float sum;

    for (y = 0; y < n; y++)
    {
        for (x = y; x < n; x++)
        {
            sum = 0.0;
            for (i = 0; i < m; i++)
                sum += pA[i * n + y] * pA[i * n + x];

            pQ[y * n + x] = sum;
            pQ[x * n + y] = sum;
        }
    }
}


static void disove_AtB(float *pA, float *pB, float *pAB, int m, int n)
{
    int y, i;
    float sum;

    for (y = 0; y < n; y++)
    {
        sum = 0;
        for (i = 0; i < m; i++)
            sum += pA[i * n + y] * pB[i];

        *pAB++ = sum;
    }
}


static void disove_rslt(float *pQ, float *pAB, float *pBB, int n)
{
    int y, i;
    float sum;

    for (y = 0; y < n; y++)
    {
        sum = 0;
        for (i = 0; i < n; i++)
            sum += pQ[y * n + i] * (float)pAB[i];

        *pBB++ = sum;
    }
}

/*
 * 运动估计
 * fp：指向特征点结构数组
 * fp_num：特征点个数
 * fimg：当前图像指针
 * pref：参考图像指针
 * i_width：图像宽度
 * i_height：图像高度
 * WinX：搜索窗口的一办宽度，必须8对齐。建议WinX最大96
，也即水平搜索窗口宽度最大不超过192个像素
 * WinY：搜索窗口的一半高度：建议WinY最大64，也即垂直搜索窗口最大宽度是128
 *
 *
 */
void GetAllMvRefine(FPOINT *fp,int fp_num,unsigned char *fimg,unsigned char *pref,int i_width, int i_height,int WinX, int WinY)
{
    int i ;
#pragma omp parallel for
    for (i = 0; i < fp_num; i++)
    { 
   	    int rx, ry,scaleXmem,scaleYmem;
   	    unsigned char *refbufpp,*cur_img;
  	    unsigned int bm[2];
	    cur_img = (unsigned char*)(fimg + (fp[i].y - 3) * i_width + (fp[i].x - 3));

	    rx = fp[i].x + fp[i].px - 3;
	    rx = rx < 0 ? 0 : rx;
	    rx = (rx + 4) & 0xfff8;
	    rx = rx > (i_width - 8)  ? i_width - 8  : rx;
	    ry = fp[i].y + fp[i].py - 3;
	    ry = ry < 0 ? 0 : ry;
	    ry = ry > (i_height - 8) ? i_height - 8 : ry;
	    rx = rx < WinX ? 0 : rx - WinX;
	    scaleXmem = (rx + WinX * 2 + 8) > i_width ? (i_width - rx) : (WinX * 2 + 8);
	    ry = ry < WinY ? 0: ry - WinY;
	    scaleYmem = (ry + WinY * 2 + 8) > i_height ? (i_height - ry) : (WinY * 2 + 8);

	    refbufpp = (unsigned char *)(pref + ry*i_width + rx);
		
        IMG_mad_8x8(refbufpp, cur_img,i_width, scaleXmem -8, scaleYmem -8, bm); 
	
	fp[i].dx = (bm[0] >> 16) + rx - fp[i].x + 3;
	fp[i].dy = (bm[0] & 0xffff) + ry - fp[i].y + 3;
    }
}

/*
    计算运动参数模型
    pAp：存放计算获得运动模型参数的结果
    fp：特征点指针，当然是指经过运动估计后的
    fpnum：特征点个数
    返回：
        0：模型参数正常
        1：输入的特征点个数太少
        2：模型参数误差比较大
*/
int GetMoveParam(affine_param *pAp,FPOINT *fp,int fpnum,unsigned char *pData)
{
    int j;
    float avgErr;

    //if(fpnum < 8) return 2;       //  输入的特征点个数太少，建议更新参考帧s
    for (j = 0; j < 5; j++)
    {
        avgErr = cal_affine_param(pAp, fp, fpnum);

        if (avgErr < 1.0)    
		break;
        fpnum = RefineFeaturePoint(fp, fpnum, avgErr);
        if (fpnum < 3)       
		break;
    }
    if (avgErr < 5.0)    return 0;  //  可靠地建立了模型，正常返回
    if (avgErr < 10.0)   return 1;  //  计算的模型参数误差比较大

    return 2;   //  计算的模型参数误差太大，更新参考帧

}


/*
    计算相似变换参数，返回拟合精度
    pAp：相似变换参数指针
    fp：特征点指针
    fpnum：特征点数目
    返回拟合精度
*/
float cal_affine_param(affine_param *pAp,FPOINT *fp,int fpnum)
{
    int i;
    const int nSize = fpnum;
    const int m = nSize * 2;
    const int n = 4;
    int idx = 0;
    float vX, vY, vD, vDist, avgErr;

    float *pMatrixA;
    float *pMatrixQ;
    float *pMatrixB;
    float *pfMatrixBB;
    float *pMatrixAB;
    void * pBuf = (void*)malloc(sizeof(float)*(m * n + n * n + m + n + n));
    pMatrixA = (float *)pBuf;
    pMatrixQ = pMatrixA + m * n;
    pMatrixB = pMatrixQ + n * n;
    pfMatrixBB = pMatrixB + m;
    pMatrixAB = pfMatrixBB + n;

    memset(pBuf, 0, (m * n + n * n + m + n + n)*sizeof(float));

    if (nSize < 8)
    {
        pAp->cos = pfMatrixBB[0] = 1.0;
        pAp->sin = pfMatrixBB[1] = 0.0;
        vX = 0;
        for (i = 0; i < nSize; i++)      vX += (float)(fp[i].dx);
        vY = 0;
        for (i = 0; i < nSize; i++)      vY += (float)(fp[i].dy);

        pfMatrixBB[2] = -vX / (float)(nSize);
        pfMatrixBB[3] = -vY / (float)(nSize);

    }
    else
    {
        for (i = 0; i < fpnum; i++)
        {
            const int offsetA = idx << 3;  //one idx prefer 8 elements
            const int offsetB = idx << 1;  //one idx prefer 2 elements

            pMatrixA[offsetA + 0] =  (float)(fp[i].x + fp[i].dx);
            pMatrixA[offsetA + 1] =  (float)(-(fp[i].y + fp[i].dy));
            pMatrixA[offsetA + 2] =  (float)1;
            pMatrixA[offsetA + 3] =  (float)0;

            pMatrixA[offsetA + 4] =  (float)(fp[i].y + fp[i].dy);
            pMatrixA[offsetA + 5] =  (float)(fp[i].x + fp[i].dx);
            pMatrixA[offsetA + 6] =  0;
            pMatrixA[offsetA + 7] =  1;

            pMatrixB[offsetB + 0] = (float)(fp[i].x);
            pMatrixB[offsetB + 1] = (float)(fp[i].y);

            ++idx;
        }

        disove_AtA(pMatrixA, pMatrixQ, m, 4);
        disove_AtB(pMatrixA, pMatrixB, pMatrixAB, m, 4);
        if (!brinv_1D(pMatrixQ, 4))
            return -1000000000.0;
        disove_rslt(pMatrixQ, pMatrixAB, pfMatrixBB, 4);
    }

    vD = 0.0;
    for (i = 0; i < fpnum; ++i)
    {
        //x' = a*x - b*y + c;
        //y' = b*x + a*y + d;
        vX = pfMatrixBB[0] * ((float)(fp[i].x + fp[i].dx)) - pfMatrixBB[1] * ((float)(fp[i].y + fp[i].dy)) + pfMatrixBB[2];
        vY = pfMatrixBB[1] * ((float)(fp[i].x + fp[i].dx)) + pfMatrixBB[0] * ((float)(fp[i].y + fp[i].dy)) + pfMatrixBB[3];
        vDist = (vX - ((float)(fp[i].x))) * (vX - ((float)(fp[i].x))) + (vY - ((float)(fp[i].y))) * (vY - ((float)(fp[i].y)));
        fp[i].dist = vDist;
        vD += vDist;
    }
    avgErr = vD / (float)fpnum;

    pAp->cos = pfMatrixBB[0];
    pAp->sin = pfMatrixBB[1];
    pAp->dx = pfMatrixBB[2];
    pAp->dy = pfMatrixBB[3];
    pAp->theta = atan(pfMatrixBB[1] / pfMatrixBB[0]);
    pAp->scale = sqrt(pfMatrixBB[1] * pfMatrixBB[1] + pfMatrixBB[0] * pfMatrixBB[0]);
	
    free(pBuf);
    return avgErr;
}


int RefineFeaturePoint(FPOINT *fp, int fpnum, float threshhold)
{
    int i, num;

    num = 0;
    for (i = 0; i < fpnum; i++)
    {
        if (fp[i].dist < threshhold)
        {
            memcpy(fp + num, fp + i, sizeof(FPOINT));
            num++;
        }
    }

    return num;
}


/*
*   按照QCIF的特征点和运动搜索，来获得CIF的特征点和运动预测, 以QCIF中的运动向量2
倍作为CIF的运动向量的预测
*   QcifFp：QCIF图像的特征点指针
*   QcifFpNum：QCIF图像的特征点个数
*   CifFp：CIF图像的特征点指针
*   返回CIF的特征点个数
*/
int SetCifFpPmvFromQcif(FPOINT *QcifFp, int QcifFpNum, FPOINT *CifFp)
{
    int n;

    for (n = 0; n < QcifFpNum; n++)
    {
        CifFp[n].px = 2 * QcifFp[n].dx;
        CifFp[n].py = 2 * QcifFp[n].dy;
    }

    return 0;
}


/*
*   按照QCIF的特征点和运动搜索，来获得CIF的特征点和运动预测, 以QCIF变换参数来预测CIF的运动向量的预测
*   QcifFp：QCIF图像的特征点指针
*   QcifFpNum：QCIF图像的特征点个数
*   CifFp：CIF图像的特征点指针
*   pAp：变换方程参数
*/
int SetCifFpPmvFromAp(FPOINT *QcifFp, int QcifFpNum, FPOINT *CifFp, affine_param *pAp)
{
    int n;

    for (n = 0; n < QcifFpNum; n++)
    {
        CifFp[n].px = ( cos(pAp->theta) * CifFp[n].x  + sin(pAp->theta) * CifFp[n].x - 2 * pAp->dx)  - CifFp[n].x;
        CifFp[n].py = ( -sin(pAp->theta) * CifFp[n].x + cos(pAp->theta) * CifFp[n].y - 2 * pAp->dy)  - CifFp[n].y;
    }

    return 0;
}


/*
*   按照CIF的特征点和运动搜索，来获得D1的特征点和运动预测,以CIF中的运动向量2倍作为D1的运动向量的预测
*   CifFp：CIF图像的特征点指针
*   CifFpNum：CIF图像的特征点个数
*   D1Fp：D1图像的特征点指针
*/
void SetD1FpPmvFromCif(FPOINT *CifFp, int CifFpNum, FPOINT *D1Fp)
{
    int n;

    for (n = 0; n < CifFpNum; n++)
    {

        D1Fp[n].px = 2 * CifFp[n].dx;
        D1Fp[n].py = 2 * CifFp[n].dy;
    }

}


/*
*   按照CIF的特征点和运动搜索，来获得D1的特征点和运动预测, 以CIF变换参数来预测D1的运动向量的预测
*   CifFp：CIF图像的特征点指针
*   CifFpNum：CIF图像的特征点个数
*   D1Fp：D1图像的特征点指针
*   pAp: 变换方程参数
*/
void SetD1FpPmvFromAp(FPOINT *CifFp, int CifFpNum, FPOINT *D1Fp, affine_param *pAp)
{
    int n;
    for (n = 0; n < CifFpNum; n++)
    {
        D1Fp[n].px = ( cos(pAp->theta) * D1Fp[n].x + sin(pAp->theta) * D1Fp[n].x - 2 * pAp->dx)  - D1Fp[n].x;
        D1Fp[n].py = (-sin(pAp->theta) * D1Fp[n].x + cos(pAp->theta) * D1Fp[n].y - 2 * pAp->dy)  - D1Fp[n].y;
    }
}


void RunMatchingPoint(stb_t* s,int* MeErr,int* MeErr_cif,int *MeErr_qcif)
{
	stb_t *ms = s;
	unsigned int t1 =  OSA_getCurTimeInMsec();
	 /*在QCIF图像上做运动搜索*/
	GetAllMvRefine(ms->QcifFp,ms->QcifFpNum,ms->fQcifCur,ms->fQcifRef,ms->i_width>>2,ms->i_height>>2,32,32);
	unsigned int t2 =  OSA_getCurTimeInMsec();
	*MeErr_qcif = GetMoveParam(&(ms->ap_qcif),ms->QcifFp,ms->QcifFpNum,ms->fQcifCur);
	unsigned int t3 =  OSA_getCurTimeInMsec();
    /*把QCIF的特征点映射到CIF上，从而获得CIF图像的特征点，并依据QCIF的运动搜索对CIF图像各特征点的运动做预计。*/
	if (*MeErr_qcif)
	{
		SetCifFpPmvFromQcif(ms->QcifFp, ms->QcifFpNum, ms->CifFp);
	}
	else
	{
		SetCifFpPmvFromAp(ms->QcifFp, ms->QcifFpNum, ms->CifFp, &ms->ap_qcif);
	}
	unsigned int t4 = OSA_getCurTimeInMsec();
	/*在CIF图像上做运动搜索*/
    GetAllMvRefine(ms->CifFp, ms->CifFpNum, ms->fCifCur, ms->fCifRef,ms->i_width >> 1, ms->i_height >> 1, 16, 16);
	unsigned int t5 =  OSA_getCurTimeInMsec();
    *MeErr_cif = GetMoveParam(&(ms->ap_cif), ms->CifFp, ms->CifFpNum, ms->fCifCur);
	unsigned int t6 =  OSA_getCurTimeInMsec();
    /*把CIF的特征点映射到D1上，从而获得D1图像的特征点，并依据CIF的运动搜索对D1图像各特征点的运动做预计*/
    if (*MeErr_cif)
    {
        SetD1FpPmvFromCif(ms->CifFp, ms->CifFpNum, ms->D1Fp);
    }
    else
    {
        SetD1FpPmvFromAp(ms->CifFp, ms->CifFpNum, ms->D1Fp, &ms->ap_cif);
    }
	unsigned int t7 =  OSA_getCurTimeInMsec();
    /* 在D1图像上做运动搜索*/
    GetAllMvRefine(ms->D1Fp, ms->D1FpNum, ms->fD1Cur->buffer[0],ms->fD1Ref, ms->i_width, ms->i_height, 16, 16);
	unsigned int t8 =  OSA_getCurTimeInMsec();
	
    /*计算运动模型参数，并剔除奇异点(产生奇异点的原因，主要是局部运动，运动估计错误等)*/
    *MeErr = GetMoveParam(&(ms->cur_af), ms->D1Fp, ms->D1FpNum,ms->fD1Cur->buffer[0]);
	unsigned int t9 =  OSA_getCurTimeInMsec();

	printf("\nt2-t1 time : %u\n",t2-t1);
	printf("t5-t4 time : %u\n",t5-t4);
	printf("t8-t7 time : %u\n",t8-t7);
	
    return;
}
