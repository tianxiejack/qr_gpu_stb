#include "cuda_runtime_api.h"
#include "device_launch_parameters.h"
#include "cuda.h"
//#include "esp_def.hpp"
#include "osa.h"


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

__global__ void kernel_RotImgProgress_uyvy(unsigned char *src, unsigned char *dst,
								float cos, float sin, float dx, float dy,
								int src_width, int src_height,
								int dst_width, int dst_height)
{

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
	dim3 block((src_width+31)/32, (src_height+31)/32);
	dim3 thread(32,32);

 	//cudaEvent_t start, stop;
 	//cudaEventCreate(&start);
 	//cudaEventCreate(&stop);
 	//cudaEventRecord(start, NULL);

	kernel_RotImgProgress_uyvy<<<block, thread>>>(src, dst,
							cos, sin, dx, dy,
							src_width, src_height,
							dst_width, dst_height);

   //cudaEventRecord(stop, NULL);
   //cudaEventSynchronize(stop);
   //float time = 0.0f;
   //cudaEventElapsedTime(&time, start, stop);

   //OSA_printf("%s:RotImgProgress_cuda time = %f ms \n", __func__, time);
}

extern "C" void Sobel_cuda(unsigned char *src, unsigned char *dst, int width, int height)
{
	dim3 block((width+31)/32, (height+31)/32);
	dim3 thread(32,32);

// 	cudaEvent_t start, stop;
// 	cudaEventCreate(&start);
// 	cudaEventCreate(&stop);
// 	cudaEventRecord(start, NULL);

	kernel_Sobel<<<block, thread>>>(src, dst, width, height);

//	cudaEventRecord(stop, NULL);
// 	cudaEventSynchronize(stop);
// 	float time = 0.0f;
// 	cudaEventElapsedTime(&time, start, stop);

// 	OSA_printf("%s:Sobel_cuda time = %f ms \n", __func__, time);

}
