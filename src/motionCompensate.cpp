 #include "stable.hpp"
#include "motionCompensate.hpp"

#include "cuda_runtime_api.h"
#include <device_launch_parameters.h>


extern "C" void RotImgProgress_cuda(unsigned char *src, unsigned char *dst, 
								float cos, float sin, float dx, float dy,
								int src_width, int src_height);



#if 1
void RotImg(unsigned char *forg,unsigned char *frot,int i_width,int i_height,float s_cos,float s_sin,float dx,float dy)
{
    float a, b, c, d;
    int m, n, p, q;
    int x, y, r_x, r_y, t_y_x, t_y_y;
    unsigned char *pdst;
    a = s_cos;
    b = s_sin;
    c = dx;
    d = dy;

#if 0
    m = (int)(a * 1024.0); //a / (a*a + b*b);
    n = (int)(b * 1024.0); //b / (a*a + b*b);
    p = (int)(-c * 1024.0) + 512; //-(a*c + b*d) / (a*a + b*b);
    q = (int)(-d * 1024.0) + 512; //-(a*d - b*c) / (a*a + b*b);
#else
    m = (int)((a / (a * a + b * b)) * 1024.0); //a / (a*a + b*b);
    n = (int)((b / (a * a + b * b)) * 1024.0); //b / (a*a + b*b);
    p = (int)(-((a * c + b * d) / (a * a + b * b)) * 1024.0) + 512; //-(a*c + b*d) / (a*a + b*b);
    q = (int)(-((a * d - b * c) / (a * a + b * b)) * 1024.0) + 512; //-(a*d - b*c) / (a*a + b*b);
#endif


    for (r_y = 0; r_y < i_height; r_y++)
    {
        t_y_x = n * r_y + p;
        t_y_y = m * r_y + q;
        pdst = frot + r_y * i_width;
        for (r_x = 0; r_x < i_width; r_x++)
        {
            x =  m * r_x + t_y_x;
            y = -n * r_x + t_y_y;
            x = x >> 10;
            y = y >> 10;
			
            if ((x < 0) || (x >= i_width) || (y < 0) || (y >= i_height)) 
            {
            		pdst[r_x] = 0;
            		continue;
            }
			
            pdst[r_x] = forg[y * i_width + x];
        }
    }	
}

#else

void RotImg(unsigned char *forg,unsigned char *frot,int i_width,int i_height,float s_cos,float s_sin,float dx,float dy)
{
    float a, b, c, d;
    int m, n, p, q;
    int x, y, r_x, r_y, t_y_x, t_y_y;
    unsigned char *pdst;
    a = s_cos;
    b = s_sin;
    c = dx;
    d = dy;

#if 0
    m = (int)(a * 1024.0); //a / (a*a + b*b);
    n = (int)(b * 1024.0); //b / (a*a + b*b);
    p = (int)(-c * 1024.0) + 512; //-(a*c + b*d) / (a*a + b*b);
    q = (int)(-d * 1024.0) + 512; //-(a*d - b*c) / (a*a + b*b);
#else
    m = (int)((a / (a * a + b * b)) * 1024.0); //a / (a*a + b*b);
    n = (int)((b / (a * a + b * b)) * 1024.0); //b / (a*a + b*b);
    p = (int)(-((a * c + b * d) / (a * a + b * b)) * 1024.0) + 512; //-(a*c + b*d) / (a*a + b*b);
    q = (int)(-((a * d - b * c) / (a * a + b * b)) * 1024.0) + 512; //-(a*d - b*c) / (a*a + b*b);
#endif


    for (r_y = 0; r_y < i_height; r_y++)
    {
        t_y_x = n * r_y + p;
        t_y_y = m * r_y + q;
        pdst = frot + r_y * i_width*3;
	#pragma omp parallel for
        for (r_x = 0; r_x < i_width; r_x++)
        {
            x =  m * r_x + t_y_x;
            y = -n * r_x + t_y_y;
            x = x >> 10;
            y = y >> 10;
			
            if ((x < 0) || (x >= i_width) || (y < 0) || (y >= i_height)) 
            {
            		pdst[r_x*3] = 0;
	   	  	pdst[r_x*3+1] = 0;
	     		pdst[r_x*3+2] = 0;
            		continue;
            }	
            pdst[r_x*3] = forg[y * i_width*3 + x*3];
	     pdst[r_x*3+1] = forg[y * i_width*3 + x*3+1];
	     pdst[r_x*3+2] = forg[y * i_width*3 + x*3+2];
        }
    }	


	return ;
}

#endif



void MotionProcess(CStability * mcs,Mat src,Mat dst,uchar mode)
{
	CStability* cs = mcs;
	stb_t* s = cs->tss;
	affine_param* ap = cs->m_modify;
	float cos = ap->cos;
	float sin = ap->sin;
	float dx = ap->dx;
	float dy = ap->dy;

	switch(mode)
	{
   	   case COMPENSATE_NORMAL:
   		   	   //do nothing;
   		   	   break;
   	   case COMPENSATE_X:
   		   	   dy = 0;
   		   	   cos = 1.0;
   		   	   sin = 0.0;
   		   	   break;
   	   case COMPENSATE_Y:
   		   	   dx = 0;
   		   	   cos = 1.0;
   		   	   sin = 0.0;
   		   	   break;
   	   case COMPENSATE_ROT:
   		   	   dx = 0.0;
   		   	   dy = 0.0;
   		   	   break;
   	   default:
   		   	   break;
	}
	//cos = 1.0;sin = 0.0;dx = 100.0;dy = 50.0;
	//RotImg(src.data,dst.data,s->i_width,s->i_height,cos,sin,dx,dy);
	
	RotImgProgress_cuda(src.data, dst.data, cos, sin, dx, dy,s->i_width, s->i_height);
}

void ImgProgress(unsigned char* src,unsigned char* dst,int nWidth,int nheight,affine_param* ap,unsigned char mode)
{
return ;
	float cos = ap->cos;
	float sin = ap->sin;
	float dx = ap->dx;
	float dy = ap->dy;

	switch(mode)
	{
   	   case COMPENSATE_NORMAL:
   		   	   //do nothing;
   		   	   break;
   	   case COMPENSATE_X:
   		   	   dy = 0;
   		   	   cos = 1.0;
   		   	   sin = 0.0;
   		   	   break;
   	   case COMPENSATE_Y:
   		   	   dx = 0;
   		   	   cos = 1.0;
   		   	   sin = 0.0;
   		   	   break;
   	   case COMPENSATE_ROT:
   		   	   dx = 0.0;
   		   	   dy = 0.0;
   		   	   break;
   	   default:
   		   	   break;
	}

	//intf("dx = %f,dy= %f\n",dx,dy);

	//unsigned char* tmp1 = (unsigned char*)malloc(720*576);
	//unsigned char* tmp2 = (unsigned char*)malloc(720*576);
	
	//cudaMemcpy(src, tmp1, 720*576, cudaMemcpyDeviceToHost);
	//RotImg(tmp1,tmp2,720,576,cos,sin,dx,dy);
	//cudaMemcpy(tmp2, dst, 720*576, cudaMemcpyHostToDevice);

	//free(tmp1);
	//free(tmp2);
	
	//RotImgProgress_cuda(src, dst, cos, sin, dx, dy,nWidth, nheight);
}



void Rotate_Tradition2(int x, int y, affine_param *ap, float *px1, float *py1)
{
	float m, n, p, q;
	float x0, y0;
	float a, b , c, d;

	a = ap->cos;
	b = ap->sin;
	c = ap->dx;
	d = ap->dy;

#if 0
	m = a;//a / (a*a + b*b);
	n = b ; //b / (a*a + b*b);
	p = -c; //-(a*c + b*d) / (a*a + b*b);
	q = -d; //-(a*d - b*c) / (a*a + b*b);
#else
	m = ((a / (a * a + b * b)) ); //a / (a*a + b*b);
	n = ((b / (a * a + b * b))); //b / (a*a + b*b);
	p = (-((a * c + b * d) / (a * a + b * b)) ); //-(a*c + b*d) / (a*a + b*b);
	q = (-((a * d - b * c) / (a * a + b * b)) ); //-(a*d - b*c) / (a*a + b*b);
#endif

// 	x0 = n * (float)y + p;
// 	y0 = m * (float)y + q;
// 
// 	*px1 = x0 + m * (double)x;
// 	*py1 = y0 - n * (double)x;
	
	*px1 = a*x - b*y + c;
	*py1 = b*x + a*y + d;
}

int getColor(float x, float y, unsigned char *c, int w)
{ 
	int a = (int)y;
	int b = (int)x;
	double dy = y-a;
	double dx = x-b;

	return (int)(0.5 + (c[a*w+b]*(1-dy)+c[(a+1)*w+b]*dy)*(1-dx)+(c[a*w+b+1]*(1-dy)+c[(a+1)*w+b+1]*dy)*dx);
}

