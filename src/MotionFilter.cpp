
#include "MotionFilter.hpp"
#include "stable.hpp"


static float IIRCoef[8][6] =
{
    {
        1.768435e-002, -3.527182e-002, 1.768435e-002,
        1.0, -1.986156e+000, 9.862525e-001
    }, // IIR1
    {
        1.763447e-002, -3.488359e-002, 1.763447e-002,
        1.0, -1.972291e+000, 9.726761e-001
    }, // IIR2
    {
        1.767785e-002, -3.382677e-002, 1.767785e-002,
        1.0, -1.944431e+000, 9.459598e-001
    }, // IIR3
    {
        1.833034e-002, -3.058438e-002, 1.833034e-002,
        1.0, -1.887696e+000, 8.937721e-001
    }, // IIR4
    {
        2.207892e-002, -1.908483e-002, 2.207892e-002,
        1.0, -1.765662e+000, 7.907348e-001
    }, // IIR5
    {
        4.104093e-002, 2.736051e-002, 4.104093e-002,
        1.0, -1.483873e+000, 5.933153e-001
    }, // IIR6
    {
        0.4, 0.0,  0.0,
        1.0, -0.6, 0.0
    }, // IIR7
    {
        0.7, 0.0, 0.0,
        1.0, -0.3, 0.0
    } // IIR8
};

/*
*   初始化滤波器，以及与运动补偿有关的参数：
*   pAp_last：上次运动补偿的参数
*   pAp_adj：当前运动补偿参数
*/
void InitFilter(affine_param *pAp_last, affine_param *pAp_adj, FILTER *flt)
{
    pAp_last->cos = 1.0;
    pAp_last->sin = 0.0;
    pAp_last->dx = 0.0;
    pAp_last->dy = 0.0;
    pAp_last->theta = 0.0;
    pAp_last->scale = 1.0;

    pAp_adj->cos =  1.0;
    pAp_adj->sin = 0.0;
    pAp_adj->dx = 0.0;
    pAp_adj->dy = 0.0;;
    pAp_adj->theta = 0.0;
    pAp_adj->scale = 1.0;

    flt->iirthetaout[0] = flt->iirthetaout[1]   = flt->iirthetaout[2]   = 0.0;
    flt->iirthetain[0]  = flt->iirthetain[1]    = flt->iirthetain[2]    = 0.0;

    flt->iirsout[0] = flt->iirsout[1] = flt->iirsout[2] = 1.0;
    flt->iirsin[0]  = flt->iirsin[1] = flt->iirsin[2]   = 1.0;

    return ;
}


void ParamFilter(CStability* mcs, affine_param *pAp_cur)
{
    float m_pMeasure[4];
    int IIRNum;

    CStability * cs = mcs;

    FILTER* flt = &(cs->tss->flt);
    HKalman hKalman = cs->tss->g_pKalman;

    /* Kalman滤波：更新yk，进行滤波，得到新的xk估计（k时刻） */
    m_pMeasure[0] = (float)pAp_cur->theta;
    m_pMeasure[1] = (float)pAp_cur->dx;
    m_pMeasure[2] = (float)pAp_cur->dy;
    m_pMeasure[3] = (float)pAp_cur->scale;

    cs->kkalman.Kalman(hKalman, m_pMeasure, NULL);/* 进行Kalman滤波 */

    IIRNum = 4;

    flt->iirthetaout[2] = flt->iirthetaout[1];
    flt->iirthetaout[1] = flt->iirthetaout[0];
    flt->iirthetain[2]   = flt->iirthetain[1];
    flt->iirthetain[1]  = flt->iirthetain[0];
    flt->iirthetain[0]  = pAp_cur->theta;
    flt->iirthetaout[0] = IIRCoef[IIRNum][0]*flt->iirthetain[0]
                          + IIRCoef[IIRNum][1]*flt->iirthetain[1]
                          + IIRCoef[IIRNum][2]*flt->iirthetain[2]
                          - IIRCoef[IIRNum][4]*flt->iirthetaout[1]
                          - IIRCoef[IIRNum][5]*flt->iirthetaout[2];

    flt->iirsout[2] = flt->iirsout[1];
    flt->iirsout[1] = flt->iirsout[0];
    flt->iirsin[2]   = flt->iirsin[1];
    flt->iirsin[1]  = flt->iirsin[0];
    flt->iirsin[0]  = pAp_cur->scale;
    flt->iirsout[0] = IIRCoef[IIRNum][0]*flt->iirsin[0]
                      + IIRCoef[IIRNum][1]*flt->iirsin[1]
                      + IIRCoef[IIRNum][2]*flt->iirsin[2]
                      - IIRCoef[IIRNum][4]*flt->iirsout[1]
                      - IIRCoef[IIRNum][5]*flt->iirsout[2];


    flt->xout   = (float)( hKalman->state_post[2] - pAp_cur->dx);
    flt->yout   = (float)( hKalman->state_post[4] - pAp_cur->dy);


}

/*
    运动滤波和运动补偿
    pAp_cur：当前解算出的运动模型参数
    pAp_last：上次做的运动补偿参数
    pAp_adj：这里存放这个函数计算出来的运动补偿参数
    MeErr：指示pAp_cur是否有效，为0表示有效。一般它是GetMoveParam（）函数的返回值
    flt：指向滤波器，这里存放着滤波器的中间状态，运动模式、滤波器参数等等重要信息
*/
int MotionFilter(CStability * mcs)
{
    float c2o_x,c2o_y, c2o_theta, c2o_scale, adj_theta, adj_scale;
    affine_param tmpAP;

   CStability * cs = mcs;
   stb_t* ms = cs->tss;
   HKalman hKalman = ms->g_pKalman;
   affine_param *pAp_cur = &(ms->cur_af);
   affine_param *pAp_last = &(ms->last_af);
   affine_param *pAp_adj = cs->m_modify;
   int MeErr = cs->MeErr;
   FILTER* flt = &(ms->flt);

    if(!MeErr)//正确地解算出了运动模型参数
	{
        //ParamFilter(hKalman, pAp_cur, flt);

        /* 计算当前图像与上次输出之间的关系*/
        c2o_scale = pAp_cur->scale / pAp_last->scale;
        c2o_theta = pAp_cur->theta - pAp_last->theta;
        c2o_x = pAp_cur->dx + pAp_last->dx;
        c2o_y = pAp_cur->dy + pAp_last->dy;
        tmpAP.scale = pAp_cur->scale;
        tmpAP.theta = pAp_cur->theta;
        tmpAP.dx = c2o_x;
        tmpAP.dy = c2o_y;

        ParamFilter(cs, &tmpAP);

        /* 计算调整量 */
        adj_theta = flt->iirthetaout[0] - c2o_theta;
        adj_scale = flt->iirsout[0]/c2o_scale;
        //pAp_adj->cos = adj_scale * cos(adj_theta); //pAp_cur->cos;
        //pAp_adj->sin = adj_scale * sin(adj_theta); //-pAp_cur->sin;
        pAp_adj->cos = cos(adj_theta); //pAp_cur->cos;
        pAp_adj->sin = sin(adj_theta); //-pAp_cur->sin;
        pAp_adj->dx =  flt->xout;
        pAp_adj->dy =  flt->yout;

        pAp_last->cos = pAp_adj->cos;
        pAp_last->sin = pAp_adj->sin;
        pAp_last->dx += pAp_cur->dx;
        pAp_last->dy += pAp_cur->dy;
        pAp_last->theta =  adj_theta;
        pAp_last->scale = adj_scale;
        if(( abs(pAp_last->dx) > 1000000.0 ) || ( abs(pAp_last->dy) > 1000000.0 ))
        {
            InitFilter(pAp_last, pAp_adj, flt);
	     cs->kkalman.KalmanInitParam(hKalman, 0.0, 0.0, 0.0, 1.0, 0.0);
        }
    }
    else
    {
        InitFilter(pAp_last, pAp_adj, flt);
        cs->kkalman.KalmanInitParam(hKalman, 0.0, 0.0, 0.0, 1.0, 0.0);
    }
    return 0;
}
