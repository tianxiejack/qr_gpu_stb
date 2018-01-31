#ifndef STABLE_HPP_
#define STABLE_HPP_

enum{
	COMPENSATE_NORMAL,		// \u5168\u6821\u6b63
	COMPENSATE_X,					// 	\u4ec5X\u65b9\u5411\u6821\u6b63
	COMPENSATE_Y,   					// \u4ec5Y\u65b9\u5411\u6821\u6b63
	COMPENSATE_ROT,				//\u4ec5\u65cb\u8f6c\u6821\u6b63
};


typedef struct
{
    float cos;                      //相似变换得旋转参数（s*cos）
    float sin;                      //相似变换的旋转参数（s*sin）
    float dx;                       //X方向的偏移
    float dy;                       //Y方向的偏移
    float scale;                    //缩放系数
    float theta;                    //旋转角度
} affine_param;


// 1920   : edge_h  	320 pixel		edge_v	180	pixel	   // distance to the edge in pixel
// 720 	: edge_h 		32   pixel		edge_v 	32	pixel
extern void Create_stable();
extern void destroy_stable();
extern void run_stable(cv::Mat src,cv::Mat dst,int nWidth,int nheight,unsigned char mode,unsigned int edge_h,unsigned int edge_v,affine_param* apout);
extern void MotionProcess(CStability * mcs,Mat src,Mat dst,uchar mode);

#endif

