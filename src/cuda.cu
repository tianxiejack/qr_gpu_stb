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


__global__ void kernel_RotImgProgress_(unsigned char *src, unsigned char *dst,							
								int src_width, int src_height,
								int m, int n,int p, int q)
{
	    int x, y;
	    unsigned char *pdst;

	    
	    const int r_x = blockDim.x * blockIdx.x + threadIdx.x;
	    const int r_y = blockDim.y * blockIdx.y + threadIdx.y;
	    if ((r_x < 0) || (r_x >= src_width) || (r_y < 0) || (r_y >= src_height)) 
		return;

	    x = m * r_x + n * r_y + p;
	    y = -n * r_x + m * r_y + q; 
	    x = x >> 10;
	    y = y >> 10;

	    if ((x < 0) || (x >= src_width) || (y < 0) || (y >= src_height)) {
	    	pdst = dst + r_y * src_width*3 + r_x*3;
      		pdst[0] = 0x00;
      		pdst[1] = 0x00;
      		pdst[2] = 0x00;
      	    }else{
	      	pdst = dst + r_y * src_width*3 + r_x*3;
	      	pdst[0] = src[y * src_width*3 + x*3];
	      	pdst[1] = src[y * src_width*3 + x*3+1];
	      	pdst[2] = src[y * src_width*3 + x*3+2]; 
          }
	   
}


/**************extern*******************/


extern "C" void tran_gray_cuda(unsigned char *src,unsigned char* dst,int nWidth,int nHeight)
{
	cuda_tran_gray<<<BLOCK_NUM,THREAD_NUM>>>(src,dst,nWidth,nHeight);
}

extern "C" void Sobel_cuda(unsigned char *src, unsigned char *dst, int width, int height)
{
	kernel_Sobel<<<BLOCK_NUM, THREAD_NUM>>>(src, dst, width, height);
}

extern "C" void RotImgProgress_cuda(unsigned char *src, unsigned char *dst, 
								float cos, float sin, float dx, float dy,
								int src_width, int src_height)
{

	dim3 block((src_width+31)/32, (src_height+31)/32);
	dim3 thread(32,32);

	float a, b, c,d ;
	int m,n,p,q;

	a = cos;
       b = sin;
       c = dx;
       d = dy;

      m = (int)((a / (a * a + b * b)) * 1024.0);
      n = (int)((b / (a * a + b * b)) * 1024.0);
      p = (int)(-((a * c + b * d) / (a * a + b * b)) * 1024.0) + 512; 
      q = (int)(-((a * d - b * c) / (a * a + b * b)) * 1024.0) + 512; 

	kernel_RotImgProgress_<<<block, thread>>>(src, dst, 
							src_width, src_height,
							m, n, p, q);

}

