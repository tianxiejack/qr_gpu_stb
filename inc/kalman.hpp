#ifndef KALMAN_H_
#define KALMAN_H_

#include <math.h>
#include <osa_debug.h>
#include <osa.h>
#include "baseData.hpp"
#include "kalmanPort.hpp"

#define FRAME_EXCHANGE(pFrame0, pFrame1)\
    {\
        UInt8 * pFrame = pFrame0 ;       \
        pFrame0 = pFrame1;      \
        pFrame1 = pFrame;        \
    }

class CKalman
{
public:
	CKalman(void);
	~CKalman(void);

public:
	int creat();
	HKalman init();
	int destroy();
	int run();

	HKalman KalmanOpen(int D, int M, int C);
	void KalmanClose(HKalman hKalman);
	void KalmanInitParam(HKalman hKalman, float theta, float delta_x, float delta_y, float scale, float DeltaT);
	int Kalman(HKalman hKalman, float *measure, float *control);
	void KalmanPredict(HKalman hKalman,  float * control );
	void KalmanCorrect(HKalman hKalman,  float * measurement );

	void MatrixMultiply(float* A, float* B, int m, int p, int n, float* C);
	void MatrixAddition(float* A, float* B, int m, int n, float* C);
	void MatrixTranspose(float* A, int m, int n, float* C);
	void MatrixSubtraction(float* A, float* B, int m, int n, float* C);
	void MatrixCopy(float *A, float *B, int m, int n);
	int MatrixBrinv( float * A, int n);
	void atrixSubtraction(float* A, float* B, int m, int n, float* C);

	float* pSpace;
};

#endif
