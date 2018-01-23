
#include "kalman.hpp"

#define false  0
#define true   1

/** \brief Floor a integer value. */
#define Utils_floor(val, align)  (((val) / (align)) * (align))

/** \brief Align a integer value. */
#define Utils_align(val, align)  Utils_floor(((val) + (align)-1), (align))
#define NULL 0
#define MEM_ILLEGAL NULL



CKalman::CKalman()
{
}

CKalman::~CKalman()
{

}

void CKalman::MatrixSubtraction(float* A, float* B, int m, int n, float* C)
{
    int i, j;
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            C[n * i + j] = A[n * i + j] - B[n * i + j];
        }
    }
}

void CKalman::MatrixTranspose(float* A, int m, int n, float* C)
{
    int i, j;
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            C[m * j + i] = A[n * i + j];
        }
    }
}


int CKalman::MatrixBrinv( float * A, int n)
{
    int i, j, k, l, u, v;
    float d, p;
    int is[64];
    int js[64];

    for ( k = 0; k < n; k++ )
    {
        d = 0.0;
        for ( i = k; i < n; i++ )
        {
            for ( j = k; j < n; j++ )
            {
                l = i * n + j;
                p = fabs(A[l]);
                if ( p > d )
                {
                    d = p;
                    is[k] = i;
                    js[k] = j;
                }
            }
        }
		if ( d+1.0 == 1.0 ) /* 锟斤拷锟斤拷为锟斤拷锟斤拷锟斤拷 */
        {
            //free( is );
            //free( js );
            return ( 0 );
        }
        if ( is[k] != k )
        {
            for ( j = 0; j < n; j++ )
            {
                u = k * n + j;
                v = is[k] * n + j;
                p = A[u];
                A[u] = A[v];
                A[v] = p;
            }
        }
        if ( js[k] != k )
        {
            for ( i = 0; i < n; i++ )
            {
                u = i * n + k;
                v = i * n + js[k];
                p = A[u];
                A[u] = A[v];
                A[v] = p;
            }
        }
        l = k * n + k;
        A[l] = 1.0 / A[l];
        for ( j = 0; j < n; j++ )
        {
            if ( j != k )
            {
                u = k * n + j;
                A[u] = A[u] * A[l];
            }
        }
        for ( i = 0; i < n; i++ )
        {
            if ( i != k )
            {
                for ( j = 0; j < n; j++ )
                {
                    if ( j != k )
                    {
                        u = i * n + j;
                        A[u] = A[u] - A[i * n + k] * A[k * n + j];
                    }
                }
            }
        }
        for ( i = 0; i < n; i++ )
        {
            if ( i != k )
            {
                u = i * n + k;
                A[u] = -A[u] * A[l];
            }
        }
    }
    for ( k = n - 1; k >= 0; k-- )
    {
        if ( js[k] != k )
        {
            for ( j = 0; j <= n - 1; j++ )
            {
                u = k * n + j;
                v = js[k] * n + j;
                p = A[u];
                A[u] = A[v];
                A[v] = p;
            }
        }
        if ( is[k] != k )
        {
            for ( i = 0; i < n; i++ )
            {
                u = i * n + k;
                v = i * n + is[k];
                p = A[u];
                A[u] = A[v];
                A[v] = p;
            }
        }
    }
    //free( is );
    //free( js );

    return (1);
}

void CKalman::MatrixCopy(float *A, float *B, int m, int n)
{
    memcpy(A, B, sizeof(float)*m * n);
}

HKalman CKalman::KalmanOpen(int D, int M, int C)
{
    HKalman hKalman =  (HKalman)malloc(sizeof(mKalman));	
    if (hKalman == MEM_ILLEGAL)
        return NULL;
		
    if ( D <= 0 || M <= 0 || C < 0 )
    {
        OSA_printf("state and measurement vectors must have positive number of dimensions! \n");
        return NULL;
    }

    hKalman->deltat = 0.04;
    hKalman->b = 0.95;
    hKalman->bK = 1;
    hKalman->m_bInited = false;
    hKalman->DP = D;
    hKalman->MP = M;
    hKalman->CP = C;

    pSpace = (float*)malloc(sizeof(float)*1500);
    float* tmpp = pSpace;
    unsigned int datablock = 0;

    hKalman->state_pre = &tmpp[datablock];
    if (hKalman->state_pre == MEM_ILLEGAL)
        return NULL;

    datablock += hKalman->DP * 1 ;
    hKalman->state_post = &tmpp[datablock];
    if (hKalman->state_post == MEM_ILLEGAL)
        return NULL;

    datablock += hKalman->DP * 1 ;
    hKalman->measure_param = &tmpp[datablock];
    if (hKalman->measure_param == MEM_ILLEGAL)
        return NULL;

    datablock += hKalman->MP * 1;
    hKalman->transition_matrix = &tmpp[datablock];
    if (hKalman->transition_matrix == MEM_ILLEGAL)
        return NULL;

    datablock += hKalman->DP * hKalman->DP;
    hKalman->process_noise_cov = &tmpp[datablock];
    if (hKalman->process_noise_cov == MEM_ILLEGAL)
        return NULL;

    datablock += hKalman->DP *  hKalman->DP;
    hKalman->measurement_matrix = &tmpp[datablock];
    if (hKalman->measurement_matrix == MEM_ILLEGAL)
        return NULL;

    datablock += hKalman->MP *  hKalman->DP ;
    hKalman->measurement_noise_cov = &tmpp[datablock];
    if (hKalman->measurement_noise_cov == MEM_ILLEGAL)
        return NULL;

    datablock += hKalman->MP * hKalman->MP ;
    hKalman->error_cov_pre = &tmpp[datablock];
    if (hKalman->error_cov_pre == MEM_ILLEGAL)
        return NULL;

    datablock += hKalman->DP * hKalman->DP ;
    hKalman->error_cov_post = &tmpp[datablock];
    if (hKalman->error_cov_post == MEM_ILLEGAL)
        return NULL;

    datablock += hKalman->DP * hKalman->DP ;
    hKalman->gain = &tmpp[datablock];
    if (hKalman->gain == MEM_ILLEGAL)
        return NULL;

    datablock += hKalman->DP * hKalman->MP ;
    if ( hKalman->CP > 0 )
    {
        hKalman->control_matrix = &tmpp[datablock];
		datablock += hKalman->DP * hKalman->CP;
        if (hKalman->control_matrix == MEM_ILLEGAL)
            return NULL;
    }

    hKalman->B_Uk = &tmpp[datablock];
    if (hKalman->B_Uk == MEM_ILLEGAL)
        return NULL;

    datablock += hKalman->DP * 1 ;
    hKalman->A_Pk = &tmpp[datablock];
    if (hKalman->A_Pk == MEM_ILLEGAL)
        return NULL;

    datablock += hKalman->DP * hKalman->DP ;
    hKalman->A_T = &tmpp[datablock];
    if (hKalman->A_T == MEM_ILLEGAL)
        return NULL;

    datablock += hKalman->DP * hKalman->DP;
    hKalman->APA_T = &tmpp[datablock];
    if (hKalman->APA_T == MEM_ILLEGAL)
        return NULL;

	datablock +=  hKalman->DP *  hKalman->DP;
    hKalman->H_T = &tmpp[datablock];
    if (hKalman->H_T == MEM_ILLEGAL)
        return NULL;

	datablock += hKalman->DP * hKalman->MP ;
    hKalman->Pk_Ht = &tmpp[datablock];
    if (hKalman->Pk_Ht == MEM_ILLEGAL)
        return NULL;

	datablock += hKalman->DP * hKalman->MP ;
    hKalman->Pk_Ht_R = &tmpp[datablock];
    if (hKalman->Pk_Ht_R == MEM_ILLEGAL)
        return NULL;

	datablock += hKalman->MP * hKalman->MP ;
    hKalman->Pk_Ht_R_Inv = &tmpp[datablock];
    if (hKalman->Pk_Ht_R_Inv == MEM_ILLEGAL)
        return NULL;

	datablock += hKalman->MP * hKalman->MP ;
    hKalman->H_Xk = &tmpp[datablock];
    if (hKalman->H_Xk == MEM_ILLEGAL)
        return NULL;

	datablock += hKalman->MP * 1 ;
    hKalman->Kk_H_Xk = &tmpp[datablock];
    if (hKalman->Kk_H_Xk == MEM_ILLEGAL)
        return NULL;

	datablock += hKalman->DP * 1 ;
    hKalman->H_Pk = &tmpp[datablock];
    if (hKalman->H_Pk == MEM_ILLEGAL)
        return NULL;

	datablock += hKalman->MP * hKalman->DP;
    hKalman->Kk_H_Pk = &tmpp[datablock];
    if (hKalman->Kk_H_Pk == MEM_ILLEGAL)
        return NULL;

	datablock += hKalman->DP * hKalman->DP ;
    hKalman->E_MUL_E = &tmpp[datablock];
    if (hKalman->E_MUL_E == MEM_ILLEGAL)
        return NULL;

	datablock += hKalman->MP * hKalman->MP ;
    hKalman->E_Kt = &tmpp[datablock];
    if (hKalman->E_Kt == MEM_ILLEGAL)
        return NULL;

	datablock += hKalman->MP * 1 ;
    hKalman->H_Pk_Ht = &tmpp[datablock];
    if (hKalman->H_Pk_Ht == MEM_ILLEGAL)
        return NULL;

	datablock += hKalman->MP * hKalman->MP ;
    hKalman->K_T = &tmpp[datablock];
    if (hKalman->K_T == MEM_ILLEGAL)
        return NULL;

	datablock += hKalman->MP * hKalman->DP ;
    hKalman->K_E_MUL_E_Kt = &tmpp[datablock];
    if (hKalman->K_E_MUL_E_Kt == MEM_ILLEGAL)
        return NULL;

	datablock += hKalman->DP * hKalman->DP ;
    hKalman->Pk_SUB_A_Pk1_At = &tmpp[datablock];
    if (hKalman->Pk_SUB_A_Pk1_At == MEM_ILLEGAL)
        return NULL;

	datablock += hKalman->DP * hKalman->DP ;
    hKalman->EE_SUB_HPkHt = &tmpp[datablock];
    if (hKalman->EE_SUB_HPkHt == MEM_ILLEGAL)
        return NULL;
	datablock += hKalman->MP * hKalman->MP ;
    hKalman->m_bInited = true;

    return hKalman;
}

void CKalman::KalmanClose(HKalman hKalman)
{
    free(pSpace);
    if (hKalman == NULL)
        return ;
    if (hKalman->m_bInited == false)
        return;
}

void CKalman::KalmanInitParam(HKalman hKalman, float theta, float delta_x, float delta_y, float scale, float DeltaT)
{
    int x, y;
    if (hKalman == NULL)
        return ;
    if (!hKalman->m_bInited)
    {
        return;
    }
    hKalman->deltat = DeltaT;
    /* 锟斤拷碳锟斤拷锟斤拷锟斤拷锟叫拷锟斤拷锟斤拷锟斤拷 */
    for ( y = 0; y < hKalman->DP; y++ )
    {
        for ( x = 0; x < hKalman->DP; x++ )
        {
            hKalman->process_noise_cov[y * hKalman->DP + x] = 0.0;
        }
    }
    hKalman->process_noise_cov[1 * hKalman->DP + 1] = 0.00001; //100.0;  /* 协锟斤拷锟筋都为100.0锟斤拷锟斤拷锟洁互锟斤拷锟斤拷 theta */
    //hKalman->process_noise_cov[3*hKalman->DP+3] = 0.00001;//100.0;  /* 协锟斤拷锟筋都为100.0锟斤拷锟斤拷锟洁互锟斤拷锟斤拷  dx */
    //hKalman->process_noise_cov[5*hKalman->DP+5] = 0.00001;//100.0;  /* 协锟斤拷锟筋都为100.0锟斤拷锟斤拷锟洁互锟斤拷锟斤拷 dy  */
    hKalman->process_noise_cov[3 * hKalman->DP + 3] = 0.00001; /* 协锟斤拷锟筋都为100.0锟斤拷锟斤拷锟洁互锟斤拷锟斤拷  dx */
    hKalman->process_noise_cov[5 * hKalman->DP + 5] = 0.00001; /* 协锟斤拷锟筋都为100.0锟斤拷锟斤拷锟洁互锟斤拷锟斤拷 dy  */
    hKalman->process_noise_cov[7 * hKalman->DP + 7] = 0.00001; //100.0;  /* 协锟斤拷锟筋都为100.0锟斤拷锟斤拷锟洁互锟斤拷锟斤拷 scale */
//  hKalman->process_noise_cov[0*hKalman->DP+0] = 10;//100.0;  /* 协锟斤拷锟筋都为100.0锟斤拷锟斤拷锟洁互锟斤拷锟斤拷 */
//  hKalman->process_noise_cov[2*hKalman->DP+2] = 10;//100.0;  /* 协锟斤拷锟筋都为100.0锟斤拷锟斤拷锟洁互锟斤拷锟斤拷 */
//  hKalman->process_noise_cov[4*hKalman->DP+4] = 10;//100.0;  /* 协锟斤拷锟筋都为100.0锟斤拷锟斤拷锟洁互锟斤拷锟斤拷 */

    /* 锟桔诧拷锟斤拷锟斤拷协锟斤拷锟斤拷锟斤拷锟?*/
    for ( y = 0; y < hKalman->MP; y++ )
    {
        for ( x = 0; x < hKalman->MP; x++ )
        {
            hKalman->measurement_noise_cov[y * hKalman->MP + x] = 0.0;
        }
    }
    hKalman->measurement_noise_cov[0 * hKalman->MP + 0] = 0.25; //0.001;  /* 协锟斤拷锟筋都为0.001锟斤拷锟斤拷锟洁互锟斤拷锟斤拷 */
    hKalman->measurement_noise_cov[1 * hKalman->MP + 1] = 0.25; //0.001;  /* 协锟斤拷锟筋都为0.001锟斤拷锟斤拷锟洁互锟斤拷锟斤拷 */
    hKalman->measurement_noise_cov[2 * hKalman->MP + 2] = 0.25; //0.001;  /* 协锟斤拷锟筋都为0.001锟斤拷锟斤拷锟洁互锟斤拷锟斤拷 */
    hKalman->measurement_noise_cov[3 * hKalman->MP + 3] = 0.25; //0.001;  /* 协锟斤拷锟筋都为0.001锟斤拷锟斤拷锟洁互锟斤拷锟斤拷 */

    /* 状态锟斤拷锟斤拷锟斤拷锟叫拷锟斤拷锟?*/
    for ( y = 0; y < hKalman->DP; y++ )
    {
        for ( x = 0; x < hKalman->DP; x++ )
        {
            hKalman->error_cov_post[y * hKalman->DP + x] = 0.0;
        }
    }
    for ( y = 0; y < hKalman->DP; y++ )
    {
        hKalman->error_cov_post[y * hKalman->DP + y] = 1.0; /* 锟皆角筹拷始协锟斤拷锟筋都为1锟斤拷锟斤拷锟洁互锟斤拷锟斤拷 */
    }

    /* 状态转锟斤拷锟斤拷 */
    for ( y = 0; y < hKalman->DP; y++ )
    {
        for ( x = 0; x < hKalman->DP; x++ )
        {
            hKalman->transition_matrix[y * hKalman->DP + x] = 0.0;
        }
    }
    for ( y = 0; y < hKalman->DP; y++ )
    {
        hKalman->transition_matrix[y * hKalman->DP + y] = 1.0; /* 锟皆斤拷为1 */
    }
    hKalman->transition_matrix[0 * hKalman->DP + 1] = 1;
    hKalman->transition_matrix[2 * hKalman->DP + 3] = 1;
    hKalman->transition_matrix[4 * hKalman->DP + 5] = 1;
    hKalman->transition_matrix[6 * hKalman->DP + 7] = 1;

    /* 锟桔诧拷锟斤拷状态锟斤拷锟斤拷锟桔诧拷锟斤拷锟斤拷转锟狡撅拷锟斤拷 */
    for ( y = 0; y < hKalman->MP; y++ )
    {
        for ( x = 0; x < hKalman->DP; x++ )
        {
            hKalman->measurement_matrix[y * hKalman->DP + x] = 0.0;
        }
    }
    hKalman->measurement_matrix[0 * hKalman->DP + 0] = 1.0;
    hKalman->measurement_matrix[1 * hKalman->DP + 2] = 1.0;
    hKalman->measurement_matrix[2 * hKalman->DP + 4] = 1.0;
    hKalman->measurement_matrix[3 * hKalman->DP + 6] = 1.0;

    // 锟桔诧拷锟斤拷锟斤拷锟斤拷锟斤拷顺锟斤拷theta, x, y
    hKalman->measure_param[0] = (float)theta;
    hKalman->measure_param[1] = (float)delta_x;
    hKalman->measure_param[2] = (float)delta_y;
    hKalman->measure_param[3] = (float)scale;
    /* 锟斤拷始锟斤拷thelta, thelta_v, x, x_v, y, y_v锟斤拷状态锟斤拷 */
    // 状态锟斤拷锟斤拷锟斤拷锟斤拷顺锟斤拷thelta, thelta_v,  x, y, x_v, y_v
    hKalman->state_post[0] = theta;
    hKalman->state_post[1] = 0.0;
    hKalman->state_post[2] = delta_x;
    hKalman->state_post[3] = 0.0;
    hKalman->state_post[4] = delta_y;
    hKalman->state_post[5] = 0.0;
    hKalman->state_post[6] = scale;
    hKalman->state_post[7] = 0.0;
}


void CKalman::KalmanPredict(HKalman hKalman,  float * control )
{
    if (hKalman == NULL)
        return;
    if (!hKalman->m_bInited)
    {
        return;
    }
    /* update the state */
    /* x'(k) = A*x(k) */
    MatrixMultiply( hKalman->transition_matrix, hKalman->state_post, hKalman->DP , hKalman->DP , 1 , hKalman->state_pre );

    if ( control != NULL && hKalman->CP > 0 )
    {
        /* x'(k) = x'(k) + B*u(k) */
        MatrixMultiply( hKalman->control_matrix, control, hKalman->DP , hKalman->CP , 1 , hKalman->B_Uk);
        MatrixAddition( hKalman->state_pre, hKalman->B_Uk, hKalman->DP, 1, hKalman->state_pre);
    }

    /* update error covariance matrices */
    /* A_Pk = A*P(k) */
    MatrixMultiply( hKalman->transition_matrix, hKalman->error_cov_post, hKalman->DP, hKalman->DP, hKalman->DP, hKalman->A_Pk);

    /* P'(k) = A_Pk*At + Q */
    MatrixTranspose(hKalman->transition_matrix, hKalman->DP, hKalman->DP, hKalman->A_T);

    MatrixMultiply(hKalman->A_Pk, hKalman->A_T, hKalman->DP, hKalman->DP, hKalman->DP, hKalman->APA_T);

    MatrixAddition(hKalman->APA_T, hKalman->process_noise_cov, hKalman->DP, hKalman->DP, hKalman->error_cov_pre);
    return ;
}

void CKalman::KalmanCorrect(HKalman hKalman,  float * measurement )
{
    int i;
    if (hKalman == NULL)
        return;
    if (!hKalman->m_bInited)
    {
        return;
    }
    if ( measurement == NULL)
    {
        ;//AfxMessageBox("Measurement is Null!!!");
    }
    if (measurement != NULL)
    {
        for (i = 0; i < hKalman->MP; i++)
        {
            hKalman->measure_param[i] = measurement[i];
        }
    }
    /* H_T = Ht*/
    MatrixTranspose( hKalman->measurement_matrix , hKalman->MP , hKalman->DP , hKalman->H_T);
    /* Pk_Ht = P'(k) * H_T */
    MatrixMultiply( hKalman->error_cov_pre, hKalman->H_T, hKalman->DP , hKalman->DP , hKalman->MP , hKalman->Pk_Ht);

    /* Pk_Ht_R = H*Pk_Ht + R */
    MatrixMultiply( hKalman->measurement_matrix , hKalman->Pk_Ht , hKalman->MP , hKalman->DP , hKalman->MP , hKalman->Pk_Ht_R);
    MatrixAddition( hKalman->Pk_Ht_R , hKalman->measurement_noise_cov , hKalman->MP , hKalman->MP , hKalman->Pk_Ht_R);

    /* Pk_Ht_R_Inv = inv(Pk_Ht_R) */
#if 0
    MatrixInversion( hKalman->Pk_Ht_R , hKalman->MP, hKalman->Pk_Ht_R_Inv);
#else
    MatrixCopy(hKalman->Pk_Ht_R_Inv, hKalman->Pk_Ht_R, hKalman->MP, hKalman->MP);
    MatrixBrinv(hKalman->Pk_Ht_R_Inv, hKalman->MP);
#endif

    /* K(k) = Pk_Ht * Pk_Ht_R_Inv  */
    MatrixMultiply( hKalman->Pk_Ht , hKalman->Pk_Ht_R_Inv, hKalman->DP , hKalman->MP , hKalman->MP , hKalman->gain);

    //update state_post
    /* H_Xk = z(k) - H*x'(k) */
    MatrixMultiply( hKalman->measurement_matrix , hKalman->state_pre , hKalman->MP , hKalman->DP , 1, hKalman->H_Xk);
    MatrixSubtraction( hKalman->measure_param , hKalman->H_Xk , hKalman->MP , 1, hKalman->H_Xk);

    /* x(k) = x'(k) + K(k)*H_Xk */
    MatrixMultiply( hKalman->gain , hKalman->H_Xk, hKalman->DP , hKalman->MP, 1, hKalman->Kk_H_Xk );
    MatrixAddition( hKalman->state_pre , hKalman->Kk_H_Xk , hKalman->DP , 1 , hKalman->state_post);

    //update error_cov_post
    /* P(k) = P'(k) - K(k)*H* P'(k) */
    MatrixMultiply( hKalman->measurement_matrix , hKalman->error_cov_pre , hKalman->MP , hKalman->DP , hKalman->DP , hKalman->H_Pk );
    MatrixMultiply( hKalman->gain , hKalman->H_Pk , hKalman->DP , hKalman->MP, hKalman->DP , hKalman->Kk_H_Pk );
    MatrixSubtraction( hKalman->error_cov_pre , hKalman->Kk_H_Pk , hKalman->DP , hKalman->DP , hKalman->error_cov_post);

}

void CKalman::MatrixMultiply(float* A, float* B, int m, int p, int n, float* C)
{
    // A = input matrix (m x p)
    // B = input matrix (p x n)
    // m = number of rows in A
    // p = number of columns in A = number of rows in B
    // n = number of columns in B
    // C = output matrix = A*B (m x n)
    int i, j, k;

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            C[n * i + j] = 0;
            for (k = 0; k < p; k++)
            {
                C[n * i + j] = C[n * i + j] + A[p * i + k] * B[n * k + j];
            }
        }
    }
	
}

void CKalman::MatrixAddition(float* A, float* B, int m, int n, float* C)
{
    // A = input matrix (m x n)
    // B = input matrix (m x n)
    // m = number of rows in A = number of rows in B
    // n = number of columns in A = number of columns in B
    // C = output matrix = A+B (m x n)
    int i, j;
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            C[n * i + j] = A[n * i + j] + B[n * i + j];
        }
    }
}


int CKalman::Kalman(HKalman hKalman, float *measure, float *control)
{
    if (hKalman == NULL)
        return -1;

    KalmanPredict(hKalman, control);
	
    KalmanCorrect(hKalman, measure);
    return 0;
}

