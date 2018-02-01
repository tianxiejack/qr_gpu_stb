#ifndef BASE_DATA_HPP_
#define BASE_DATA_HPP_

//#define NULL 0
#include "opencv2/core/core.hpp"
#include "opencv2/core/types_c.h"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/features2d/features2d.hpp"

#include "kalmanPort.hpp"

using namespace cv;
#define CUDA_MEM			0
#define MAX_WIDTH 			 1920
#define MAX_HEIGHT 			 1080
/*
#define POIINTSTHRESH		 5000
#define POINTS_NUM			 100
#define CUR					 		 0
#define REF							 1
#define EDGETHRESH          	 35
#define MAXVALTHRESH       100
#define MINVALTHRESH       80
*/
enum{
	COMPENSATE_NORMAL,		// 全校正
	COMPENSATE_X,					// 	仅X方向校正
	COMPENSATE_Y,   					// 仅Y方向校正
	COMPENSATE_ROT,				//仅旋转校正
};

typedef struct
{
	int i_stride[3];
	unsigned char *buffer[3];
	unsigned char *a;
} stb_frame_t;

typedef struct
{
    float cos;                      //相似变换得旋转参数（s*cos）
    float sin;                      //相似变换的旋转参数（s*sin）
    float dx;                       //X方向的偏移
    float dy;                       //Y方向的偏移
    float scale;                    //缩放系数
    float theta;                    //旋转角度
} affine_param;

typedef struct
{
    int ftv;               		   //特征值
    int x;                          //X坐标
    int y;                          //Y坐标
    int dx;                       //X方向运动向量
    int dy;                       //Y方向的运动向量
    int px;                       //X方向的预测运动向量
    int py;                       //Y方向的预测运动向量
    float dist;                  //该点与模型参数方程解算出的对应点间的距离
} FPOINT;

typedef struct
{
    float xout;                     //x方向运动滤波输出。[0]是滤波输出。
    float yout;                     //y方向运动滤波输出。[0]是滤波输出。
    float thetaout;                 //角度的滤波输出，[0]是滤波输出。
    float sout;                     //缩放倍数的滤波输出，[0]是滤波输出。

    float iirthetaout[3];           //角度的滤波输出，[0]是滤波输出。
    float iirsout[3];               //缩放倍数的滤波输出，[0]是滤波输出。
    float iirthetain[3];            //角度的滤波输入，[0]是滤波输出。
    float iirsin[3];                //缩放倍数的滤波输入，[0]是滤波输出。
} FILTER;



typedef struct
{
    void *pInnerBuf;                //1920 x 48 内存，地址是128对齐，在片内申请的一片内存，用于中间数据缓冲，以提高运行速度
    void *pInnerKalman;             //4096+512 内存，地址是128对齐，在片内申请的一片内存，用于中间数据缓冲，以提高运行速度

    int i_width;                    //图像宽（单位为像素）
    int i_height;                   //图像高（单位为像素）
    int grid_w;                     //特征点分布的栅格尺寸宽
    int grid_h;                     //特征点分布的栅格尺寸高

   unsigned char edgeTh;           //作为边缘的特征点的阈值，是由FindFeatruePoint函数计算出来的，其值的范围0x18--0x80之间

    stb_frame_t *fD1Cur;            //当前帧
    stb_frame_t *fD1Out;            //经过稳像后输出的D1图像帧

    unsigned char *fD1Ref;          //D1参考图像
    unsigned char *fD1RefFilt;      //经过中值滤波的D1参考图像
    unsigned char *fD1RefSobel;     //sobel边缘检测的D1参考图像
    unsigned char *fCifRef;         //CIF参考图像
    unsigned char *fCifRefFilt;     //经过中值滤波的CIF参考图像
    unsigned char *fCifRefSobel;    //sobel边缘检测的CIF参考图像
    unsigned char *fQcifRef;        //QCIF参考图像
    unsigned char *fQcifRefFilt;    //经过中值滤波的CIF参考图像
    unsigned char *fQcifRefSobel;   //sobel边缘检测的CIF参考图像

    unsigned char *fD1CurFilt;      //经过中值滤波的D1当前图像
    unsigned char *fD1CurSobel;     //sobel边缘检测的D1当前图像
    unsigned char *fCifCur;         //CIF当前图像
    unsigned char *fCifCurFilt;     //经过中值滤波的CIF当前图像
    unsigned char *fCifCurSobel;    //sobel边缘检测的CIF当前图像
    unsigned char *fQcifCur;        //QCIF当前图像
    unsigned char *fQcifCurFilt;    //经过中值滤波的QCIF当前图像
    unsigned char *fQcifCurSobel;   //sobel边缘检测的QCIF当前图像

    FPOINT *D1Fp;                   //D1图像的特征点
    int D1FpNum;                    //D1图像的特征点个数
    FPOINT *CifFp;                  //CIF图像的特征点
    int CifFpNum;                   //CIF图像的特征点个数
    FPOINT *QcifFp;                 //QCIF图像的特征点
    int QcifFpNum;                  //QCIF图像的特征点个数


    affine_param cur_af;            //当前由运动估计获得的运动模型参数
    affine_param last_af;           //保存上次图像所作的运动补偿参数
    affine_param bak_af;            //当前由运动估计获得的运动模型参数
    affine_param ap_qcif;
    affine_param ap_cif;

    FILTER flt;                     //滤波器

    HKalman g_pKalman;              //卡尔曼句柄

    int stat[8];
    volatile int bReset;
    volatile int sideBySide;
    volatile int featurePoint;
    volatile int elasped_interval;
}stb_t;

#endif
