#ifndef KALMAN_PORT_H_
#define KALMAN_PORT_H_

typedef struct _Kalman_t
{
	void* memId;
	int MP;  //number of measure vector dimensions  4
	int DP;  //number of state   vector dimensions  8
	int CP;  //number of control vector dimensions  0

	float  *state_pre;           //[DP * 1]8*1
	float* state_post;           //[DP * 1]8*1
	float* measure_param;        //[MP * 1]4*1
	float* transition_matrix;    //[DP * DP]8*8
	float* control_matrix;       //[DP * CP]8*0 if (CP > 0)
	float* measurement_matrix;   //[MP * DP]4*8
	float* process_noise_cov;    //[DP * DP]8*8
	float* measurement_noise_cov;//[MP * MP]4*4
	float* error_cov_pre;        //[DP * DP]8*8
	float* error_cov_post;       //[DP * DP]8*8
	float* gain;                 //[DP * MP]8*4  sum 356

	float* E_MUL_E;                 //[MP * MP]4*4
	float* E_Kt;                    //[1  * MP]1*4
	float* H_Pk_Ht;                 //[MP * MP]4*4
	float* EE_SUB_HPkHt;            //[MP * MP]4*4
	float*  K_T;                    //[MP * DP]4*8
	float* K_E_MUL_E_Kt;            //[DP * DP]8*8
	float*  Pk_SUB_A_Pk1_At;        //[DP * DP]8*8 sum 212
	float  dK;
	float   b, bK;                  // forgetting factor

	float   deltat; // 锟斤拷频锟斤拷锟斤拷时锟斤拷锟斤拷
	int     m_bInited;

	//temporary  variable
	float *B_Uk;         //[ DP * 1  ]8*1
	float *A_Pk;         //[ DP * DP ]8*8
	float *A_T;          //[ DP * DP ]8*8
	float *APA_T;        //[ DP * DP ]8*8

	float *H_T;         //[ DP * MP ]8*4
	float *Pk_Ht;       //[ DP * MP ]8*4
	float *Pk_Ht_R;     //[ MP * MP ]4*4
	float *Pk_Ht_R_Inv; //[ MP * MP ]4*4
	float *H_Xk;        //[ MP * 1  ]4*1
	float *Kk_H_Xk;     //[ DP * 1  ]8*1
	float *H_Pk;        //[ MP * DP ]4*8
	float *Kk_H_Pk;     //[ DP * DP ]8*8 sum 404
}mKalman;

typedef mKalman* HKalman;


#endif