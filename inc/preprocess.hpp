#ifndef PRE_PROCESS_HPP_
#define PRE_PROCESS_HPP_

#include "baseData.hpp"

void IMG_sobel(const unsigned char* in, unsigned char* out,short cols, short rows);
void pyramid(Mat in,Mat out,stb_t* s);
void pyramid_cv(Mat in,Mat out,stb_t* s);
void preprocess(Mat fcur,Mat cifCur,Mat QcifCur,Mat fCurSobel,Mat cifCurSobel,Mat QcifCurSobel,stb_t* s);
void analytime();


#endif
