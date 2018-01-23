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
	int RunStabilize(Mat src,Mat dst,int nWidth, int nHeight,uchar mode);
	void showPoints(unsigned char code);

	int MeErr, MeErr_qcif, MeErr_cif;
	int h_mb_num,v_mb_num;
	stb_t* tss;
	CStability* pThis;
	affine_param* m_modify;

public:
	CKalman kkalman;
	CFindFtp kfindFtp;
	Mat mfout,mfcur,mfcur_sobel,mfCur_ref,mfcur_sobel_ref;
	Mat mfCifCur,mfCifCur_sobel,mfCifCur_ref,mfCifCur_sobel_ref;
	Mat mfQcifCur,mfQCifCur_sobel,mfQCifCur_ref,mfQCifCur_sobel_ref;
};

#endif
