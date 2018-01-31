#include "cuda_runtime_api.h"
#include "device_launch_parameters.h"
#include "cuda.h"
//#include "esp_def.hpp"
#include "osa.h"

#define THREAD_NUM	32
#define BLOCK_NUM	32


__global__ static void cuda_tran_gray(unsigned char *src,unsigned char* dst,int nWidth,int nHeight)
{
	int ImgHeight, ImgWidth;
	int i;
	ImgWidth = nWidth;
	ImgHeight = nHeight;
	uint8_t *  pDst8_t;
	uint8_t *  pSrc8_t;

	pSrc8_t = (uint8_t*)(src);
	pDst8_t = (uint8_t*)(dst);

	//const int x = blockDim.x * blockIdx.x + threadIdx.x;
	const int tid = threadIdx.x;
	const int bid = blockIdx.x;

	//for(int x = 0; x < ImgHeight*ImgWidth; x++)
	for(i = bid*THREAD_NUM + tid;i<ImgHeight*ImgWidth;i += BLOCK_NUM * THREAD_NUM)
	{
		pDst8_t[i] = pSrc8_t[i*3+1];
	}
	
}

#if 0
__global__ static void YUVDownSample(unsigned char *yuv_buff,unsigned char *yuv2_buff,int PIC_W, int PIC_H) 
{
	int h,v;
	for(v=0; v<PIC_H/2; v++)
	{
		for(h=0; h<PIC_W/2; h++)
		{
			yuv2_buff[v*PIC_W/2+h] = yuv_buff[v*2*PIC_W+h*2];
		}
	}
}

#endif



__global__ void kernel_Sobel(unsigned char *src, unsigned char *dst, int width, int height)
{
	const int x = blockDim.x * blockIdx.x + threadIdx.x;
	const int y = blockDim.y * blockIdx.y + threadIdx.y;

	int Gx, Gy;
	int sobel;

	Gx = Gy = 0;

	//if(x>=1 && x<(width-1) && y>=1 && y<(height-1))
	{
		Gx = src[(y-1)*width + x-1]*(-1) + src[(y-1)*width + x]*(0) + src[(y-1)*width + x+1]*(1) +
		      src[y*width     + x-1]*(-2) + src[y*width     + x]*(0) + src[y*width     + x+1]*(2) +
		      src[(y+1)*width + x-1]*(-1) + src[(y+1)*width + x]*(0) + src[(y+1)*width + x+1]*(1);

/*		Gy = src[(y-1)*width + x-1]*(-1) + src[(y-1)*width + x]*(-2) + src[(y-1)*width + x+1]*(-1) +
		      src[y*width     + x-1]*(0)  + src[y*width     + x]*(0)  + src[y*width     + x+1]*(0) +
	  	      src[(y+1)*width + x-1]*(1)  + src[(y+1)*width + x]*(2)  + src[(y+1)*width + x+1]*(1);
*/
		 Gy = src[(y-1)*width + x-1]*(1) + src[(y-1)*width + x]*(2) + src[(y-1)*width + x+1]*(1) +
		      src[y*width     + x-1]*(0)  + src[y*width     + x]*(0)  + src[y*width     + x+1]*(0) +
	  	      src[(y+1)*width + x-1]*(-1)  + src[(y+1)*width + x]*(-2)  + src[(y+1)*width + x+1]*(-1);

		sobel = (int)sqrt((float)(Gx * Gx + Gy * Gy));
		//sobel = sobel < 20 ? 0 : sobel;
		dst[y*width + x]  = sobel;
	}
	//else
	{
		//dst[y*width + x] = 0;
	}

}

__global__ void kernel_RotImgProgress_gray(unsigned char *src, unsigned char *dst,
								float cos, float sin, float dx, float dy,
								int src_width, int src_height,
								int dst_width, int dst_height)
{
	const int x = blockDim.x * blockIdx.x + threadIdx.x;
	const int y = blockDim.y * blockIdx.y + threadIdx.y;

	float a, b , c, d;
	int x1, y1;
	float x2, y2;
	int yi, xi;
	float detx, dety;
	int centerx, centery;

	a = cos;
	b = sin;
	c = dx;
	d = dy;

	//centre point is (centerx, centery)
	centerx = src_width/2;//0;//src_width/2;//
	centery = src_height/2;//0;//src_height/2;//

	x1 = x - centerx;
	y1 = y - centery;
#if 0
	x2 = a*x1 - b*y1 + c + centerx;
	y2 = b*x1 + a*y1 + d + centery;
#else
	x2 = a*x1 + b*y1 + c + centerx;
	y2 = -b*x1 + a*y1 + d + centery;
#endif
	if(x2>=0 && x2<dst_width && y2>=0 && y2<dst_height)
	{
		xi = (int)x2;
		yi = (int)y2;
		detx = x2-xi;
		dety = y2-yi;
		dst[y*dst_width + x] = (unsigned char)(0.5 + (src[yi*src_width + xi]*(1-dety) + src[(yi+1)*src_width + xi]*dety)*(1-detx) +
								(src[yi*src_width + xi + 1]*(1-dety) + src[(yi+1)*src_width + xi + 1]*dety)*detx);
	}
	else
	{
		dst[y*dst_width + x] = 0;
	}
}


__global__ void kernel_RotImgProgress_(unsigned char *src, unsigned char *dst, 
								float cos, float sin, float dx, float dy,
								int src_width, int src_height, 
								int dst_width, int dst_height)
{
	const int x = blockDim.x * blockIdx.x + threadIdx.x;
	const int y = blockDim.y * blockIdx.y + threadIdx.y;


	float a, b , c, d;
	float x2, y2;
	int yi, xi;
	float detx, dety;

	a = cos;
	b = sin;
	c = dx;
	d = dy;

	x2 = a*x - b*y + c;
	y2 = b*x + a*y + d;

	if(x2>=0 && x2<dst_width && y2>=0 && y2<dst_height)
	{
		xi = (int)x2;
		yi = (int)y2;
		detx = x2-xi;
		dety = y2-yi;
		dst[y*dst_width + x] = (unsigned char)(0.5 + (src[yi*dst_width + xi]*(1-dety) + src[(yi+1)*dst_width + xi]*dety)*(1-detx) +
								(src[yi*dst_width + xi + 1]*(1-dety) + src[(yi+1)*dst_width + xi + 1]*dety)*detx);
	}
	else
	{
		dst[y*dst_width + x] = 0;
	}
}



__global__ void kernel_RotImgProgress_uyvy(unsigned char *src, unsigned char *dst,
								float cos, float sin, float dx, float dy,
								int src_width, int src_height,
								int dst_width, int dst_height)
{
	    float a, b, c, d;
	    int m, n, p, q;
	    int x, y, r_x, r_y, t_y_x, t_y_y;
	    int i_width,i_height;
	    unsigned char *pdst;
	    a = cos;
	    b = sin;
	    c = dx;
	    d = dy;
	    i_width = src_width;
	    i_height = src_height;

	    m = (int)((a / (a * a + b * b)) * 1024.0);
	    n = (int)((b / (a * a + b * b)) * 1024.0);
	    p = (int)(-((a * c + b * d) / (a * a + b * b)) * 1024.0) + 512; 
	    q = (int)(-((a * d - b * c) / (a * a + b * b)) * 1024.0) + 512; 

	    const int tid = threadIdx.x;
	    const int bid = blockIdx.x;
		
	    for (r_y = 0; r_y < i_height; r_y ++)
	    {
	        t_y_x = n * r_y + p;
	        t_y_y = m * r_y + q;
	        pdst = dst + r_y * i_width*3;
			
	        for (r_x = bid*THREAD_NUM + tid; r_x < i_width*3; r_x += BLOCK_NUM*THREAD_NUM)
	        {
	            x =  m * r_x + t_y_x;
	            y = -n * r_x + t_y_y;
	            x = x >> 10;
	            y = y >> 10;
				
	            if ((x < 0) || (x >= i_width*3) || (y < 0) || (y >= i_height)) 
	            		pdst[r_x] = 0;
		     else
				pdst[r_x] = src[y * i_width*3 + x];
	        }
	    }	
}


extern "C" void RotImgProgress_gray_cuda(unsigned char *src, unsigned char *dst,
								float cos, float sin, float dx, float dy,
								int src_width, int src_height,
								int dst_width, int dst_height)
{
	dim3 block((src_width+31)/32, (src_height+31)/32);
	dim3 thread(32,32);

 	//cudaEvent_t start, stop;
 	//cudaEventCreate(&start);
 	//cudaEventCreate(&stop);
 	//cudaEventRecord(start, NULL);

	kernel_RotImgProgress_gray<<<block, thread>>>(src, dst,
							cos, sin, dx, dy,
							src_width, src_height,
							dst_width, dst_height);

   //cudaEventRecord(stop, NULL);
   //cudaEventSynchronize(stop);
   //float time = 0.0f;
   //cudaEventElapsedTime(&time, start, stop);

   //OSA_printf("%s:RotImgProgress_cuda time = %f ms \n", __func__, time);
}

extern "C" void RotImgProgress_uyvy_cuda(unsigned char *src, unsigned char *dst,
								float cos, float sin, float dx, float dy,
								int src_width, int src_height,
								int dst_width, int dst_height)
{
	kernel_RotImgProgress_uyvy<<<BLOCK_NUM, THREAD_NUM>>>(src, dst,
							cos, sin, dx, dy,
							src_width, src_height,
							dst_width, dst_height);
}

extern "C" void Sobel_cuda(unsigned char *src, unsigned char *dst, int width, int height)
{
	kernel_Sobel<<<BLOCK_NUM, THREAD_NUM>>>(src, dst, width, height);
}



extern "C" void tran_gray_cuda(unsigned char *src,unsigned char* dst,int nWidth,int nHeight)
{
	cuda_tran_gray<<<BLOCK_NUM,THREAD_NUM>>>(src,dst,nWidth,nHeight);
}



extern "C" void RotImgProgress_cuda(unsigned char *src, unsigned char *dst, 
								float cos, float sin, float dx, float dy,
								int src_width, int src_height, 
								int dst_width, int dst_height)
{

	dim3 block((src_width+31)/32, (src_height+31)/32);
	dim3 thread(32,32);
	kernel_RotImgProgress_<<<block, thread>>>(src, dst, 
							cos, sin, dx, dy,
							src_width, src_height, dst_width, dst_height);

}






