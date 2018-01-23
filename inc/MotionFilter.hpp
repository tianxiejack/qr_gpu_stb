#ifndef MOTION_FILTER_HPP_
#define MOTION_FILTER_HPP_

#include "baseData.hpp"
#include "stable.hpp"

int MotionFilter(CStability * mcs);
void ParamFilter(HKalman hKalman, affine_param *pAp_cur, FILTER *flt);
void InitFilter(affine_param *pAp_last, affine_param *pAp_adj, FILTER *flt);

#endif
