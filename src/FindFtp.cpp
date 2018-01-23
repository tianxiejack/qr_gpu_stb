
#include "FindFtp.hpp"
#include "stdio.h"

CFindFtp::CFindFtp()
{
}

CFindFtp::~CFindFtp()
{
}

unsigned int CFindFtp::FindMaxVal(unsigned char*  dat, int width, int height, int stride)
{
	int i, j;
	unsigned char max_v;
	unsigned int max_n, ix, iy;
	unsigned char *psrc;
	max_v = 0;
	max_n = 0;
	ix = iy = 0;
	for(j=0; j<height; j++){
		psrc = dat+j*stride;
		for(i=0; i<width; i++){
			if(psrc[i] > max_v){
				max_v = psrc[i];
				ix = i;
				iy = j;
			}
		}
	}
	 max_n = iy*width+ix;
	return (max_n<<8) | max_v;
}

void CFindFtp::showPoints(Mat qcif,Mat cif,Mat d1,FPOINT* fp)
{
	int n;
	string pNum;
	unsigned int num[3] = {0};
	Point keypoints;
	for(n = 0;n<h_mb_num*v_mb_num;n++)
	{
		if(qcif_fp[n].ftv)
		{
			keypoints.x = qcif_fp[n].x;
			keypoints.y = qcif_fp[n].y;
			pNum = format("%d",num[0]);
			circle(qcif,keypoints,1,Scalar(255,0,255),0.5,8);
			//putText(qcif,pNum , Point(keypoints.x,keypoints.y), 0, 0.4, Scalar::all(255), 0.8, 8);
			num[0]++;
			//printf(" num = %d\n",num[0]);
		}
		if(cif_fp[n].ftv)
		{
			keypoints.x = cif_fp[n].x;
			keypoints.y = cif_fp[n].y;
			pNum = format("%d",n);
			circle(cif,keypoints,2,Scalar(255,0,255),1,8);
			//putText(cif, pNum, Point(keypoints.x,keypoints.y), 0, 0.4, Scalar::all(255), 1, 8);
			num[1]++;
		}
		if(d1_fp[n].ftv)
		{
			keypoints.x = d1_fp[n].x;
			keypoints.y = d1_fp[n].y;
			pNum = format("%d",n);
			circle(d1,keypoints,4,Scalar(255,0,255),2,8);
			//putText(d1, pNum, Point(keypoints.x,keypoints.y), 0, 0.7, Scalar::all(255), 2, 8);
			num[2]++;
		}
	}
	printf("###QCif_num = %d **** Cif_num = %d *** d1_num = %d\n",num[0],num[1],num[2]);

	namedWindow("showPoint in qcif",CV_WINDOW_NORMAL);
	namedWindow("showPoint in cif");
	namedWindow("showPoint in d1");
	imshow("showPoint in qcif",qcif);
	imshow("showPoint in cif",cif);
	imshow("showPoint in d1",d1);

	//waitKey(10000);
	
	return ;
}




int CFindFtp::findFtpInQcif()
{
    int i,j,n;
    //printf("v_mb_num - edge_v = %d\n",v_mb_num - edge_v);
    //printf("h_mb_num - edge_h = %d\n",h_mb_num - edge_h);
    pp = 0;
   // pSobelBuf[pp] =  (unsigned char*)(qcif_sobel + edge_v * mb_h * qcif_w + edge_h * mb_w);
   pSobelBuf =  (unsigned char*)qcif_sobel;
   qcifnum = 0;
    for (j = edge_v; j < (v_mb_num - edge_v); j++)
	{
	        for ( i = edge_h; i < (h_mb_num - edge_h); i++)
	        {
	            maxv = FindMaxVal(pSobelBuf+j*mb_h*qcif_w+i*mb_w, mb_w, mb_h, qcif_w );
	            max_n = maxv >> 8;
	            maxv &= 0xff;
	            MaxSobel[maxv]++;
	            max_y = max_n / mb_w;
	            max_x = max_n % mb_w;
	            if (maxv > edge_thr)
	            {
	                qcif_fp[j * h_mb_num + i].ftv = maxv;
	                qcif_fp[j * h_mb_num + i].x   = i * mb_w + max_x;
	                qcif_fp[j * h_mb_num + i].y   = j * mb_h + max_y;
	                qcifnum++;
	            }
	            else
	            {
	                qcif_fp[j * h_mb_num + i].ftv = 0;
	            }
	        }
	    }
	if (qcifnum < 3)
	{
	    *edgeTh = edge_thr;
	    *qcif_num = 0;
	    *cif_num = 0;
	    *d1_num = 0;
	    return -1;
	}
}

int CFindFtp::findFtpInCif()
{
	int n;
	for (n = 0; n < h_mb_num * v_mb_num; n++)
	    {
	        if (qcif_fp[n].ftv)
	        {        
	            cif_fp[n].ftv = qcif_fp[n].ftv;
	            /*
	            *   因为特征点是在QCIF图像上找的，而要在CIF图像上画坐标，这期间经过了1
	            次下采样和一次sobel
	            *   CIF-〉QCIF：    x' = 2*x +1 + 0.5, y' = 2*y + 1 + 0.5
	            *   QCIF-〉SOBEL：  x' = x, y' = y + 1
	            *   所以，最终有：x' = 2*x + 1.5; y' = 2*y+3.5
	            *   所以有：
	            *       CifFp[n].x = 2 * QcifFp[n].x + 2;//1.5
	            *       CifFp[n].y = 2 * QcifFp[n].y + 4;//3.5
	            *
	            */
	            cif_fp[n].x = 2 * qcif_fp[n].x;//+ 2;
	            cif_fp[n].y = 2 * qcif_fp[n].y + 2;//+ 4;
	        }
	        else
	        {
	            cif_fp[n].ftv = 0;
	        }
	    }

	  cifnum = ArrFtp(cif_sobel, cif_w,edge_thr, cif_fp, h_mb_num, v_mb_num);
	    if (cifnum < 3)
	    {
	        *edgeTh = edge_thr;
	        *qcif_num = 0;
	        *cif_num = 0;
	        *d1_num = 0;
	        return -1;
	    }


}

/*
    * 3: 将cif的基本特征点赋予D1，作为D1的基本特征点，并将其纠正到D1边缘上
 */
int CFindFtp::findFtpInD1()
{
	int n;
	cifnum = 0;
    for (n = 0; n < h_mb_num * v_mb_num; n++)
    {
            if (cif_fp[n].ftv)
            {
            	  cifnum ++;
                /*
                *   因为特征点是在CIF图像上找的，而要在D1图像上画坐标，这期间经过了1次下采样和1
                次sobel
                *   D1-〉CIF：  x' = 2*x +1 + 0.5, y' = 2*y + 1 + 0.5
                *   CIF-〉SOBEL：   x' = x, y' = y + 1
                *   所以，最终有：x' = 2*x + 1.5; y' = 2*y+3.5
                *   所以有：
                *       D1Fp[n].x = 2 * CifFp[n].x + 1;//1.5
                *       D1Fp[n].y = 2 * CifFp[n].y + 4;//3.5
                */
                d1_fp[n].ftv = cif_fp[n].ftv;
                d1_fp[n].x = 2 * cif_fp[n].x;// + 2;
                d1_fp[n].y = 2 * cif_fp[n].y + 2;
            }
            else
            {
                d1_fp[n].ftv = 0;
            }
        }

    d1num = ArrFtp(d1_sobel, d1_w,edge_thr, d1_fp, h_mb_num, v_mb_num);
    if (d1num < 3)
    {
        *edgeTh = edge_thr;
        *qcif_num = 0;
        *cif_num = 0;
        *d1_num = 0;
        return -1;
    }

}


int CFindFtp::tidyFtp()
{
    int i,j,ft_num;
    int h_thr,v_thr;
    if (d1num > MAX_FT_NUM)
    {
            max_y = (h_mb_num * v_mb_num) >> 1;
            max_n = 0;
            for (i = 0; i < 256; i++)
            {
                max_n += MaxSobel[255 - i];
                if (max_n > max_y) 
			break;
            }
            maxv = 255 - i;
            maxv = maxv > 0x80 ? 0x80 : maxv;
            maxv = maxv < 0x18 ? 0x18 : maxv;
            *edgeTh = maxv; //获得了边缘的阈值，在0x18和0x80之间。
            if (maxv > 0x18)
            {
                ft_num = 0;
                for (j = edge_v; j < (v_mb_num - edge_v); j++)
                {
                    for ( i = edge_h; i < (h_mb_num - edge_h); i++)
                    {
                        if (d1_fp[j * h_mb_num + i].ftv < maxv)
                            d1_fp[j * h_mb_num + i].ftv = 0;
                        else
                            ft_num++;
                    }
                }
            }
            else
                ft_num = d1num;
            //剔出挨得比较近的点
            for (j = edge_v; j < (v_mb_num - edge_v); j++)
            {
                if (ft_num < MAX_FT_NUM)    
			break;
                if (d1_fp[j * h_mb_num + edge_h].ftv)
                {
                    d1_fp[j * h_mb_num + edge_h].ftv = 0;
                    ft_num--;
                }

                if (ft_num < MAX_FT_NUM)    
			break;

                if (d1_fp[j * h_mb_num + (h_mb_num - edge_h - 1)].ftv)
                {
                    d1_fp[j * h_mb_num + (h_mb_num - edge_h - 1)].ftv = 0;
                    ft_num--;
                }
            }

            h_thr = 2;
            v_thr = 2;
            for (;;)
            {
                if (ft_num < MAX_FT_NUM)     break;

                for (j = edge_v; j < (v_mb_num - edge_v); j++)
                {
                    if (ft_num < MAX_FT_NUM)     break;
                    for ( i = edge_h; i < (h_mb_num - edge_h); i++)
                    {
                        if (i % h_thr)
                        {
                            if (d1_fp[j * h_mb_num + i].ftv)
                            {
                                d1_fp[j * h_mb_num + i].ftv = 0;
                                ft_num--;
                                if (ft_num < MAX_FT_NUM)     break;
                            }
                        }
                    }
                }
                h_thr <<=  1;

                for (j = edge_v; j < (v_mb_num - edge_v); j++)
                {
                    if (ft_num < MAX_FT_NUM)     break;
                    if ((j % v_thr) == 0) continue;
                    for ( i = edge_h; i < (h_mb_num - edge_h); i++)
                    {
                        if (d1_fp[j * h_mb_num + i].ftv)
                        {
                            d1_fp[j * h_mb_num + i].ftv = 0;
                            ft_num--;
                        }
                        if (ft_num < MAX_FT_NUM)     break;
                    }
                }
                v_thr <<= 1;
            }
            d1num = ft_num;
        }
	
    for (j = 0; j < v_mb_num; j++)
      {
          for ( i = 0; i < h_mb_num; i++)
          {
              cif_fp[j * h_mb_num + i].ftv = d1_fp[j * h_mb_num + i].ftv;
              qcif_fp[j * h_mb_num + i].ftv = d1_fp[j * h_mb_num + i].ftv;
          }
      }
      cifnum = d1num;
      qcifnum  = d1num;
      if (d1num < 3)
      {
          *edgeTh = edge_thr;
          *qcif_num = 0;
          *cif_num = 0;
          *d1_num = 0;
          return -1;
      }
}

int CFindFtp::arrangeFtp()
{
    int n = 0,i,j;	
    for (j = edge_v; j < (v_mb_num - edge_v); j++)
    {
        if (n == d1num)  break;
        for ( i = edge_h; i < (h_mb_num - edge_h); i++)
        {
            if (d1_fp[j * h_mb_num + i].ftv)
            {
                d1_fp[n].ftv = d1_fp[j * h_mb_num + i].ftv;		
                d1_fp[n].x = d1_fp[j * h_mb_num + i].x;
                d1_fp[n].y = d1_fp[j * h_mb_num + i].y;
		  d1_fp[j * h_mb_num + i].ftv = 0;

		//insert about the cif		
                cif_fp[n].ftv = cif_fp[j * h_mb_num + i].ftv;
                cif_fp[n].x = cif_fp[j * h_mb_num + i].x;
                cif_fp[n].y = cif_fp[j * h_mb_num + i].y;

		//insert about the qcif
                qcif_fp[n].ftv = qcif_fp[j * h_mb_num + i].ftv;
                qcif_fp[n].x = qcif_fp[j * h_mb_num + i].x;
                qcif_fp[n].y = qcif_fp[j * h_mb_num + i].y;
				
                n++;
                if (n == d1num)  break;
            }
	    else
	    {
		  cif_fp[j * h_mb_num + i].ftv = 0;
		  qcif_fp[j * h_mb_num + i].ftv = 0;
	    }
        }  
    }
	
    *edgeTh = edge_thr;
    *qcif_num = qcifnum;
    *cif_num = cifnum;
    *d1_num = d1num;
     return -1;
}

/*
*
将从低一级变换过来的特征点坐标放到边缘上，并剔出一些重复的点，重复的原因是由于?
匦掠成湓斐傻?
*   pp0,pp1：片内乒乓缓冲
*   pedge：边缘检测图像
*   i_width：图像宽度
*   edgeTh：边缘阈值
*   fp：特征点结构
*   h_bm_num：水平栅格的个数
*   v_bm_num：垂直栅格的个数
*   返回特征点个数
*/
int CFindFtp::ArrFtp(unsigned char *pedge,
    int i_width,unsigned char edgeTh,FPOINT *fp,
    int h_mb_num, int v_mb_num)
{
	int idx;
	int x,y,dx,dy;
	unsigned int maxval;
	int fpnum = 0;

	pSobelBuf = (unsigned char*)pedge;

	for (y = 0; y < v_mb_num; y++)
	{
	        for (x = 0; x < h_mb_num; x++)
	        {
	            if (fp[y * h_mb_num + x].ftv == 0)
	            	continue;

	            maxval =	 FindMaxVal(pSobelBuf+(fp[y* h_mb_num + x].y-4)*i_width + (fp[y* h_mb_num + x].x-6) , 12,  8, i_width);
	            max_n = maxval >> 8;
	            maxval &= 0xff;
	            if (maxval > edgeTh)
	            {
	                dx = max_n % 12;
	                dy = max_n / 12;
	                fp[y * h_mb_num + x].ftv = maxval;
	                fp[y * h_mb_num + x].x = (fp[y * h_mb_num + x].x-6) + dx;
	                fp[y * h_mb_num + x].y = fp[y * h_mb_num + x].y-4 + dy;
	                fpnum ++;
	            }
	            else
	            {
	                fp[y * h_mb_num + x].ftv = 0;
	            }
	        }
	    }

    if (fpnum == 0)  
		return 0;
    b_add = fpnum;
    for (x = 0; x < (h_mb_num * v_mb_num - 1); x++)
    {
        if (fp[x].ftv == 0)  continue;
        for (y = x + 1; y < h_mb_num * v_mb_num; y++)
        {
            if ((fp[x].x == fp[y].x) && (fp[x].y == fp[y].y) && (fp[y].ftv != 0))
            {
                fp[y].ftv = 0;
                b_add --;
            }
        }
    }
    return b_add;
}

/**
    * @Fun: FindFtp
    * @brief:  在QCIF,CIF,D1上寻找特征点：
    * @ 将输入的sobel图像划分成一个个mb_w(宽) * mb_h(高)的栅格,
    * @ 在每个栅格里寻求最大sobel值的点作为特征点,
    * @ 如果栅格中没有边缘(即sobel值小于阈值),则认为这个栅格里没有特征点.
    * @ 程序自动剔除挨得比较近的特征点,特征点的坐标和特征值(sobel值)填写到fp结构。
    * @
    * @param:edgeTh：用来返回判断为边缘的阈值，阈值最小0x18，最大0x80
    * @param:qcif_fp：存放QCIF中找到的特征点。
    * @                         fp一定要能够容纳可能的最多的特征点个数
    * @                         （即其大小不能小于(qcif_w/mb_w) * (qcif_h/sbm_h)个FPOINT）
    * @ param:qcif_sobel：QCIF的sobel边缘检测图像
    * @ param:qcif_num：返回QCIF图像的特征点个数
    * @ param:qcif_w：QCIF图像的宽度
    * @ param:qcif_h：QCIF图像的高度
    * @ param:cif_w：CIF图像的宽度
    * @ param:d1_w：D1图像的宽度
    * @ param:mb_w：栅格的宽度。要求w/mb_w能整除
    * @ param:mb_h：栅格的高度。要求h/mb_w能整除
    * @ param:edge_h：在左右边缘保留这个宽度不进行特征点的检测，单位是栅格的个数
    * @ param:edge_v：在上下边缘保留这个宽度不进行特征点的检测，单位是栅格的个数
    * @ param:cif_fp：存放CIF中找到的特征点。
    * @                 fp 一定要能够容纳可能的最多的特征点个数
    * @                 （即其大小不能小于(qcif_w/mb_w) * (qcif_h/bm_h)个FPOINT）
    * @ param:cif_sobel：CIF的sobel边缘检测图像
    * @ param:cif_num：返回CIF图像的特征点个数
    * @ param:d1_fp：存放D1中找到的特征点。
    * @             fp 一定要能够容纳可能的最多的特征点个数
    * @             （即其大小不能小于(qcif_w/mb_w) * (qcif_h/bm_h)个FPOINT）
    * @ param:d1_sobel：D1的sobel边缘检测图像
    * @ param:d1_num：返回D1图像的特征点个数
*/
void CFindFtp::findFtp(unsigned char *cedgeTh,FPOINT *cqcif_fp,
	    int *cqcif_num,int cqcif_w, int cqcif_h,int ccif_w,  int cd1_w,int cmb_w,   int cmb_h,int cedge_h, int cedge_v,
	    FPOINT *ccif_fp,int *ccif_num,FPOINT *cd1_fp,  int *cd1_num,
	    Mat fCurSobel,Mat CifSobel,Mat QCifSobel)
{
    int i, j, n;
    int next_i, next_j;
    int ErrNo = 0;

    unsigned int h_thr, v_thr;
    unsigned int  b_add; 
    qcifnum = 0 ;
    cifnum = 0;
    d1num = 0;
   
    MaxSobel = (int *)QCifSobel.data;
    edgeTh = cedgeTh;
    qcif_fp = cqcif_fp;
    qcif_sobel = QCifSobel.data;
    qcif_num = cqcif_num;
    qcif_w = cqcif_w;
    qcif_h = cqcif_h;
    cif_w = ccif_w;
    d1_w = cd1_w;
    mb_w = cmb_w;
    mb_h =cmb_h;
    edge_h = cedge_h;
    edge_v = cedge_v;
    cif_fp =ccif_fp;
    cif_sobel = CifSobel.data;
    cif_num = ccif_num;
    d1_fp = cd1_fp;
    d1_sobel = fCurSobel.data;
    d1_num = cd1_num;

    edge_thr = 0x18;    //边缘的最小阈值
    h_mb_num = qcif_w / mb_w;
    v_mb_num = qcif_h / mb_h;

    ErrNo = findFtpInQcif();
    if(ErrNo == -1)
    	return ;

    ErrNo = findFtpInCif();
    if(ErrNo == -1)
    	return ;

    ErrNo = findFtpInD1();
    if(ErrNo == -1)
    	return ;


    ErrNo = tidyFtp();
    if(ErrNo == -1)
    {
    	return ;
    }

    ErrNo = arrangeFtp();

//showPoints(QCifSobel,CifSobel,fCurSobel,NULL);
    return ;
}

