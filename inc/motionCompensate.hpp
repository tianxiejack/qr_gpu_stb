#ifndef MOTIONCOMPENSATE_HPP_
#define MOTIONCOMPENSATE_HPP_

void RotImg(unsigned char *forg,unsigned char *frot,int i_width,int i_height,float s_cos,float s_sin,float dx,float dy);
void MotionProcess(CStability* mcs,Mat src,Mat dst,uchar mode);

void Rotate_Tradition2(int x, int y, affine_param *ap, float *px1, float *py1);
int getColor(float x, float y, unsigned char *c, int w);
void RotImgProgress1(unsigned char *src, unsigned char *dst, int width,  int height,affine_param *ap, int width_dst, int height_dst);


#endif
