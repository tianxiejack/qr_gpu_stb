#ifndef MOTIONCOMPENSATE_HPP_
#define MOTIONCOMPENSATE_HPP_

void RotImg(unsigned char *forg,unsigned char *frot,int i_width,int i_height,float s_cos,float s_sin,float dx,float dy);
void MotionProcess(CStability* mcs,Mat src,Mat dst,uchar mode);

void Rotate_Tradition2(int x, int y, affine_param *ap, float *px1, float *py1);
int getColor(float x, float y, unsigned char *c, int w);

void ImgProgress(unsigned char* src,unsigned char* dst,int nWidth,int nheight,affine_param* ap,unsigned char mode);


#endif
