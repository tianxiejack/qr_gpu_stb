
#include "preprocess.hpp"

#include "opencv2/core/core.hpp"
#include "opencv2/core/types_c.h"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/features2d/features2d.hpp"
#include "stdio.h"
#include "osa.h"

using namespace cv;

#define DEBUGTIME	500

unsigned int timepoint[5];
static unsigned int anytime[5] = {0};
static unsigned int anytimenum = 0;
static unsigned int mintime[5];
static unsigned int maxtime[5] = {0};
static unsigned int avr[5];


void IMG_sobel(const unsigned char* in, unsigned char* out,short cols, short rows)
{
		/* ------------------------------------------------ */
              /*  Intermediate values.                            */
              /* ------------------------------------------------ */
              int H;    /* Horizontal mask result                 */
              int V;    /* Vertical mask result                   */
              int O;    /* Sum of horizontal and vertical masks   */
              int i;    /* Input pixel offset                     */
              int o;    /* Output pixel offset.                   */
              int xy;   /* Loop counter.                          */

              /* ------------------------------------------------ */
              /*  Input values.                                   */
              /* ------------------------------------------------ */
              int i00, i01, i02;
              int i10,      i12;
              int i20, i21, i22;

              /* ------------------------------------------------ */
              /*  Step through the entire image.  We step         */
              /*  through 'rows - 2' rows in the output image,    */
              /*  since those are the only rows that are fully    */
              /*  defined for our filter.                         */
              /* ------------------------------------------------ */
              for (xy = 0, i = cols + 1, o = 1;
                   xy < cols*(rows-2) - 2;
                   xy++, i++, o++)
              {

                  /* -------------------------------------------- */
                  /*  Read necessary data to process an input     */
                  /*  pixel.  The following instructions are      */
                  /*  written to reflect the position of the      */
                  /*  input pixels in reference to the pixel      */
                  /*  being processed, which would correspond     */
                  /*  to the blank space left in the middle.      */
                  /* -------------------------------------------- */
                  i00=in[i-cols-1]; i01=in[i-cols]; i02=in[i-cols+1];
                  i10=in[i     -1];                 i12=in[i     +1];
                  i20=in[i+cols-1]; i21=in[i+cols]; i22=in[i+cols+1];

                  /* -------------------------------------------- */
                  /*  Apply the horizontal mask.                  */
                  /* -------------------------------------------- */
                  H = -i00 - 2*i01 -   i02 +   i20 + 2*i21 + i22;

                  /* -------------------------------------------- */
                  /*  Apply the vertical mask.                    */
                  /* -------------------------------------------- */
                  V = -i00 +   i02 - 2*i10 + 2*i12 -   i20 + i22;

                  O = abs(H) + abs(V);

                  /* -------------------------------------------- */
                 /*  If the result is over 255 (largest valid    */
                  /*  pixel value), saturate (clamp) to 255.      */
                  /* -------------------------------------------- */
                  if (O > 255) O = 255;

                  /* -------------------------------------------- */
                  /*  Store the result.                           */
                  /* -------------------------------------------- */
                  out[o] = O;
              }
 }



//void pyramid(unsigned char*in_data,unsigned char*out_data,stb_t* s)
void pyramid(Mat in,Mat out,stb_t* s)
{
	unsigned char* in_data = in.data;
	unsigned char* out_data = out.data;

	/******	gaussPyr    ******/
	/*  	高斯下采样待补充      */

}

void pyramid_cv(Mat mfcur,Mat mfCifCur,Mat mQcifCur,stb_t* s)
{	
	cv::pyrDown(mfcur, mfCifCur, Size(s->i_width>>1, s->i_height>>1));
	cv::pyrDown(mfCifCur, mQcifCur, Size(s->i_width>>2, s->i_height>>2));
	
	return ;
}

void preprocess(Mat fcur,Mat cifCur,Mat QcifCur,Mat fCurSobel,Mat cifCurSobel,Mat QcifCurSobel,stb_t* s)
{
	unsigned int tm1,tm2,tm3,tmp;
	static unsigned int tt1 = 0,tt2 = 0;
	/*	 pyramid		*/
	timepoint[0] = OSA_getCurTimeInMsec();
	pyramid_cv(fcur,cifCur,QcifCur,s);
	timepoint[1] = OSA_getCurTimeInMsec();
	/*		sobel 	*/
	IMG_sobel(fcur.data, fCurSobel.data,fcur.cols, fcur.rows);
	IMG_sobel(cifCur.data, cifCurSobel.data,cifCur.cols, cifCur.rows);
	IMG_sobel(QcifCur.data, QcifCurSobel.data,QcifCur.cols, QcifCur.rows);
	timepoint[2] = OSA_getCurTimeInMsec();

	analytime();
	return ;
}


void analytime()
{
	int i;
	unsigned int tmp;

	if(anytimenum == DEBUGTIME)
	{
		for(i = 0;i<2;i++)
			avr[i] = anytime[i]/DEBUGTIME;
		
		anytimenum = 0;

		for(i = 0;i<2;i++)
		{
			printf(" min[%d] = %u\n",i,mintime[i]);
			printf(" avr[%d] = %u\n",i, avr[i]);	
			printf(" max[%d] = %u\n",i,maxtime[i]);			
		}			
	}
	else if(anytimenum >DEBUGTIME)
	{
		anytimenum = 0;	
	}
	else if(anytimenum == 0)
	{
		memset(anytime,0,sizeof(unsigned int));
		memset(mintime,10000,5*sizeof(unsigned int));
		memset(maxtime,0,5*sizeof(unsigned int));
	}
	else
	{
	
		for(i = 0;i<2;i++)
		{
			tmp = (timepoint[i+1] - timepoint[i]);
			anytime[i] += tmp;
			
			if(tmp < mintime[i])
				mintime[i] = tmp;		
			else if(tmp>maxtime[i])
				maxtime[i] = tmp;
		}	
	}
	anytimenum ++ ;

	return ;	
}

	
