#ifndef MATCHING_POINT_H_
#define MATCHING_POINT_H_

#include "baseData.hpp"

void IMG_mad_8x8(const unsigned char * refImg,const unsigned char * srcImg,int pitch, int sx, int sy,unsigned int * match);

int brinv_1D(float a[], int n);
static void disove_AtA(float *pA, float *pQ, int m, int n);
static void disove_AtB(float *pA, float *pB, float *pAB, int m, int n);
static void disove_rslt(float *pQ, float *pAB, float *pBB, int n);

void GetAllMvRefine(FPOINT *fp,int fp_num,unsigned char *fimg,unsigned char *pref,int i_width, int i_height,int WinX, int WinY);
int RefineFeaturePoint(FPOINT *fp, int fpnum, float threshhold);
float cal_affine_param(affine_param *pAp,FPOINT *fp,int fpnum);
int SetCifFpPmvFromQcif(FPOINT *QcifFp, int QcifFpNum, FPOINT *CifFp);
int SetCifFpPmvFromAp(FPOINT *QcifFp, int QcifFpNum, FPOINT *CifFp, affine_param *pAp);
void SetD1FpPmvFromCif(FPOINT *CifFp, int CifFpNum, FPOINT *D1Fp);
void SetD1FpPmvFromAp(FPOINT *CifFp, int CifFpNum, FPOINT *D1Fp, affine_param *pAp);

void RunMatchingPoint(stb_t* s,int *MeErr,int* MeErr_cif,int *MeErr_qcif);


#endif
