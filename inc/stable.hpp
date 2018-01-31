#ifndef STABLE_HPP_
#define STABLE_HPP_

#include "baseData.hpp"
#include "kalman.hpp"
#include "FindFtp.hpp"


class CStability
{
public:
	CStability();
	 ~CStability();

	void init();
	void destroy();
	void run();

	void allocspace();
	void OpenStabilize(stb_t *s);
	void CloseStabilize(stb_t *s);
	int RunStabilize(Mat src,Mat dst,int nWidth, int nHeight,uchar mode,unsigned int edge_h,unsigned int edge_v,affine_param* apout);
	void showPoints(unsigned char code);

	void analytime();
		

	int MeErr, MeErr_qcif, MeErr_cif;
	int h_mb_num,v_mb_num;
	stb_t* tss;
	CStability* pThis;
	affine_param* m_modify;

public:
	unsigned int time[20];
	 unsigned long anytime[20] = {0};
	 unsigned int anytimenum = 0;
	 unsigned int mintime[20];
	 unsigned int maxtime[20] = {0};
	 unsigned int avr[20];
	 
	unsigned char edge_v,edge_h;
	CKalman kkalman;
	CFindFtp kfindFtp;
	Mat mfout,mfcur,mfcur_sobel,mfCur_ref,mfcur_sobel_ref;
	Mat mfCifCur,mfCifCur_sobel,mfCifCur_ref,mfCifCur_sobel_ref;
	Mat mfQcifCur,mfQCifCur_sobel,mfQCifCur_ref,mfQCifCur_sobel_ref;
};

#endif
