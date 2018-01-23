#ifndef FIND_FTP_HPP
#define FIND_FTP_HPP

#include "baseData.hpp"

#define MAX_FT_NUM 64

class CFindFtp
{
public:
	CFindFtp();
	~CFindFtp();

	void findFtp(unsigned char *cedgeTh,FPOINT *cqcif_fp,
	    int *cqcif_num,int cqcif_w, int cqcif_h,int ccif_w,  int cd1_w,int cmb_w,   int cmb_h,int cedge_h, int cedge_v,
	    FPOINT *ccif_fp, int *ccif_num,FPOINT *cd1_fp,int *cd1_num,Mat fCurSobel,Mat CifSobel,Mat QCifSobel);

	int findFtpInQcif();
	int findFtpInCif();
	int findFtpInD1();
	int tidyFtp();
	int arrangeFtp();
	void showPoints(Mat qcif,Mat cif,Mat d1,FPOINT* fp);
	
	int ArrFtp(unsigned char *pedge,int i_width,unsigned char edgeTh,FPOINT *fp,
	    int h_mb_num, int v_mb_num);

	unsigned int FindMaxVal(unsigned char*  dat, int width, int height, int stride);
	unsigned char* pBuffer;

private:
	unsigned int edge_thr;
	unsigned int qcifnum,cifnum,d1num;
	int h_mb_num,v_mb_num;
	unsigned int maxv;
	int max_n, max_x, max_y;
	int next_i,next_j;
	unsigned int b_add;
    	unsigned char* pSobelBuf;

protected:
	unsigned char *edgeTh;
	unsigned char *qcif_sobel;
	unsigned char *cif_sobel;
	unsigned char *d1_sobel;
	int *cif_num;
	int *d1_num;
	int *qcif_num;
	int qcif_h,qcif_w,cif_w,d1_w;
	int mb_w,mb_h,edge_h,edge_v;
	FPOINT * qcif_fp;
	FPOINT *cif_fp;
	FPOINT *d1_fp;

	FPOINT * qcif_bak_fp;
	FPOINT *cif_bak_fp;
	FPOINT *d1_bak_fp;
	
	int *MaxSobel;
	int pp;
};


#endif
