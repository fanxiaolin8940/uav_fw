/*
 * DynModel.c
 *
 * Code generation for model "DynModel".
 *
 * Model version              : 1.341
 * Simulink Coder version : 8.8 (R2015a) 09-Feb-2015
 * C source code generated on : Thu Jul 21 22:30:55 2016
 *
 * Target selection: grt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64
 * Emulation hardware selection:
 *    Differs from embedded hardware (Generic->16-bit Embedded Processor)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "DynModel.h"
#include "DynModel_private.h"

/* Block signals (auto storage) */
B_DynModel_T DynModel_B;

/* Continuous states */
X_DynModel_T DynModel_X;

/* Block states (auto storage) */
DW_DynModel_T DynModel_DW;

/* External inputs (root inport signals with auto storage) */
ExtU_DynModel_T DynModel_U;

/* External outputs (root outports fed by signals with auto storage) */
ExtY_DynModel_T DynModel_Y;

/* Real-time model */
RT_MODEL_DynModel_T DynModel_M_;
RT_MODEL_DynModel_T *const DynModel_M = &DynModel_M_;

/*
 * This function updates continuous states using the ODE4 fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  time_T t = rtsiGetT(si);
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE4_IntgData *id = (ODE4_IntgData *)rtsiGetSolverData(si);
  real_T *y = id->y;
  real_T *f0 = id->f[0];
  real_T *f1 = id->f[1];
  real_T *f2 = id->f[2];
  real_T *f3 = id->f[3];
  real_T temp;
  int_T i;
  int_T nXc = 21;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);

  /* Save the state values at time t in y, we'll use x as ynew. */
  (void) memcpy(y, x,
                (uint_T)nXc*sizeof(real_T));

  /* Assumes that rtsiSetT and ModelOutputs are up-to-date */
  /* f0 = f(t,y) */
  rtsiSetdX(si, f0);
  DynModel_derivatives();

  /* f1 = f(t + (h/2), y + (h/2)*f0) */
  temp = 0.5 * h;
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (temp*f0[i]);
  }

  rtsiSetT(si, t + temp);
  rtsiSetdX(si, f1);
  DynModel_step();
  DynModel_derivatives();

  /* f2 = f(t + (h/2), y + (h/2)*f1) */
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (temp*f1[i]);
  }

  rtsiSetdX(si, f2);
  DynModel_step();
  DynModel_derivatives();

  /* f3 = f(t + h, y + h*f2) */
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (h*f2[i]);
  }

  rtsiSetT(si, tnew);
  rtsiSetdX(si, f3);
  DynModel_step();
  DynModel_derivatives();

  /* tnew = t + h
     ynew = y + (h/6)*(f0 + 2*f1 + 2*f2 + 2*f3) */
  temp = h / 6.0;
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + temp*(f0[i] + 2.0*f1[i] + 2.0*f2[i] + f3[i]);
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

real_T rt_urand_Upu32_Yd_f_pw(uint32_T *u)
{
  uint32_T lo;
  uint32_T hi;

  /* Uniform random number generator (random number between 0 and 1)

     #define IA      16807                      magic multiplier = 7^5
     #define IM      2147483647                 modulus = 2^31-1
     #define IQ      127773                     IM div IA
     #define IR      2836                       IM modulo IA
     #define S       4.656612875245797e-10      reciprocal of 2^31-1
     test = IA * (seed % IQ) - IR * (seed/IQ)
     seed = test < 0 ? (test + IM) : test
     return (seed*S)
   */
  lo = *u % 127773UL * 16807UL;
  hi = *u / 127773UL * 2836UL;
  if (lo < hi) {
    *u = 2147483647UL - (hi - lo);
  } else {
    *u = lo - hi;
  }

  return (real_T)*u * 4.6566128752457969E-10;
}

real_T rt_nrand_Upu32_Yd_f_pw(uint32_T *u)
{
  real_T y;
  real_T sr;
  real_T si;

  /* Normal (Gaussian) random number generator */
  do {
    sr = 2.0 * rt_urand_Upu32_Yd_f_pw(u) - 1.0;
    si = 2.0 * rt_urand_Upu32_Yd_f_pw(u) - 1.0;
    si = sr * sr + si * si;
  } while (si > 1.0);

  y = sqrt(-2.0 * log(si) / si) * sr;
  return y;
}

real_T rt_roundd(real_T u)
{
  real_T y;
  if (fabs(u) < 4.503599627370496E+15) {
    if (u >= 0.5) {
      y = floor(u + 0.5);
    } else if (u > -0.5) {
      y = 0.0;
    } else {
      y = ceil(u - 0.5);
    }
  } else {
    y = u;
  }

  return y;
}

real_T rt_modd(real_T u0, real_T u1)
{
  real_T y;
  real_T tmp;
  if (u1 == 0.0) {
    y = u0;
  } else {
    tmp = u0 / u1;
    if (u1 <= floor(u1)) {
      y = u0 - floor(tmp) * u1;
    } else if (fabs(tmp - rt_roundd(tmp)) <= DBL_EPSILON * fabs(tmp)) {
      y = 0.0;
    } else {
      y = (tmp - floor(tmp)) * u1;
    }
  }

  return y;
}

void rt_mrdivide_U1d1x3_U2d3x3_Yd1x3(const real_T u0[3], const real_T u1[9],
  real_T y[3])
{
  real_T A[9];
  int16_T r1;
  int16_T r2;
  int16_T r3;
  real_T maxval;
  real_T a21;
  int16_T rtemp;
  memcpy(&A[0L], &u1[0L], 9U * sizeof(real_T));
  r1 = 0;
  r2 = 1;
  r3 = 2;
  maxval = fabs(u1[0]);
  a21 = fabs(u1[1]);
  if (a21 > maxval) {
    maxval = a21;
    r1 = 1;
    r2 = 0;
  }

  if (fabs(u1[2]) > maxval) {
    r1 = 2;
    r2 = 1;
    r3 = 0;
  }

  A[r2] = u1[r2] / u1[r1];
  A[r3] /= A[r1];
  A[3 + r2] -= A[3 + r1] * A[r2];
  A[3 + r3] -= A[3 + r1] * A[r3];
  A[6 + r2] -= A[6 + r1] * A[r2];
  A[6 + r3] -= A[6 + r1] * A[r3];
  if (fabs(A[3 + r3]) > fabs(A[3 + r2])) {
    rtemp = r2 + 1;
    r2 = r3;
    r3 = rtemp - 1;
  }

  A[3 + r3] /= A[3 + r2];
  A[6 + r3] -= A[3 + r3] * A[6 + r2];
  y[r1] = u0[0] / A[r1];
  y[r2] = u0[1] - A[3 + r1] * y[r1];
  y[r3] = u0[2] - A[6 + r1] * y[r1];
  y[r2] /= A[3 + r2];
  y[r3] -= A[6 + r2] * y[r2];
  y[r3] /= A[6 + r3];
  y[r2] -= A[3 + r3] * y[r3];
  y[r1] -= y[r3] * A[r3];
  y[r1] -= y[r2] * A[r2];
}

/* Model step function */
void DynModel_step(void)
{
  /* local block i/o variables */
  real_T rtb_q0q1q2q3[4];
  real_T b_y;
  real_T c_y;
  real_T d_y;
  real_T rtb_Saturation;
  real_T rtb_IntegratorSecondOrder_o2_j;
  real_T rtb_IntegratorSecondOrder_o1_i;
  real_T rtb_Switch_i;
  real_T rtb_Saturation_i;
  real_T rtb_Saturation_gu;
  real_T rtb_Saturation1;
  real_T rtb_Switch_fk;
  real_T rtb_Sum_j;
  real_T rtb_Switch_d;
  real_T rtb_TrigonometricFunction;
  real_T rtb_Sum_e;
  real_T rtb_Add5[3];
  real_T rtb_ixk;
  real_T rtb_jxk;
  real_T rtb_Product_k4[3];
  real_T rtb_Saturation_kg[3];
  real_T rtb_Sum_p[3];
  real_T rtb_Sum1_p;
  real_T rtb_Add6_0[3];
  int16_T rtb_Compare_0;
  real_T rtb_Sum_h_idx_0;
  real_T rtb_Sum_h_idx_1;
  real_T rtb_Sum_h_idx_2;
  real_T rtb_Add2_idx_0;
  real_T rtb_Add2_idx_1;
  real_T rtb_Add2_idx_2;
  real_T rtb_Saturation_l_idx_2;
  real_T rtb_Saturation_l_idx_1;
  real_T rtb_Saturation_l_idx_0;
  real_T rtb_Saturation_g_idx_0;
  real_T rtb_Saturation_g_idx_1;
  real_T rtb_Saturation_g_idx_2;
  real_T rtb_ZeroOrderHold2_idx_0;
  real_T rtb_ZeroOrderHold2_idx_1;
  real_T rtb_ZeroOrderHold2_idx_2;
  real32_T rtb_DataTypeConversion12_idx_0;
  real32_T rtb_DataTypeConversion12_idx_1;
  real32_T rtb_DataTypeConversion12_idx_2;
  real32_T rtb_DataTypeConversion14_idx_0;
  real32_T rtb_DataTypeConversion14_idx_1;
  real32_T rtb_DataTypeConversion14_idx_2;
  real32_T rtb_DataTypeConversion15_idx_0;
  real32_T rtb_DataTypeConversion15_idx_1;
  real32_T rtb_DataTypeConversion15_idx_2;
  real_T u0;
  real_T u1;
  real_T u0_0;
  if (rtmIsMajorTimeStep(DynModel_M)) {
    /* set solver stop time */
    if (!(DynModel_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&DynModel_M->solverInfo,
                            ((DynModel_M->Timing.clockTickH0 + 1) *
        DynModel_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&DynModel_M->solverInfo,
                            ((DynModel_M->Timing.clockTick0 + 1) *
        DynModel_M->Timing.stepSize0 + DynModel_M->Timing.clockTickH0 *
        DynModel_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(DynModel_M)) {
    DynModel_M->Timing.t[0] = rtsiGetT(&DynModel_M->solverInfo);
  }

  /* Integrator: '<S4>/xe,ye,ze' */
  DynModel_B.xeyeze[0] = DynModel_X.xeyeze_CSTATE[0];
  DynModel_B.xeyeze[1] = DynModel_X.xeyeze_CSTATE[1];
  DynModel_B.xeyeze[2] = DynModel_X.xeyeze_CSTATE[2];

  /* Sum: '<S54>/Sum1' incorporates:
   *  UnaryMinus: '<S54>/Ze2height'
   */
  DynModel_B.Sum1 = -DynModel_B.xeyeze[2];

  /* Saturate: '<S56>/Limit  altitude  to troposhere' */
  if (DynModel_B.Sum1 > 11000.0) {
    rtb_Switch_i = 11000.0;
  } else if (DynModel_B.Sum1 < 0.0) {
    rtb_Switch_i = 0.0;
  } else {
    rtb_Switch_i = DynModel_B.Sum1;
  }

  /* Sum: '<S56>/Sum1' incorporates:
   *  Constant: '<S56>/Sea Level  Temperature'
   *  Gain: '<S56>/Lapse Rate'
   *  Saturate: '<S56>/Limit  altitude  to troposhere'
   */
  DynModel_B.Sum1_e = 288.15 - 0.0065 * rtb_Switch_i;
  if (rtmIsMajorTimeStep(DynModel_M)) {
    /* Sum: '<S63>/Sum2' incorporates:
     *  Constant: '<S63>/K2C'
     *  RandomNumber: '<S63>/Random Number'
     *  Sum: '<S63>/Add1'
     */
    rtb_Saturation = (273.15 + DynModel_B.Sum1_e) + DynModel_DW.NextOutput;

    /* Saturate: '<S63>/Saturation' */
    if (rtb_Saturation > 85.0) {
      rtb_Saturation = 85.0;
    } else {
      if (rtb_Saturation < 40.0) {
        rtb_Saturation = 40.0;
      }
    }

    /* End of Saturate: '<S63>/Saturation' */
  }

  /* Gain: '<S56>/1//T0' */
  rtb_IntegratorSecondOrder_o2_j = 0.00347041471455839 * DynModel_B.Sum1_e;

  /* Math: '<S56>/(T//T0)^(g//LR) ' */
  if (rtb_IntegratorSecondOrder_o2_j < 0.0) {
    rtb_IntegratorSecondOrder_o1_i = -pow(-rtb_IntegratorSecondOrder_o2_j,
      5.2558756014667134);
  } else {
    rtb_IntegratorSecondOrder_o1_i = pow(rtb_IntegratorSecondOrder_o2_j,
      5.2558756014667134);
  }

  /* End of Math: '<S56>/(T//T0)^(g//LR) ' */

  /* Saturate: '<S56>/Limit  altitude  to Stratosphere' incorporates:
   *  Constant: '<S56>/Altitude of Troposphere'
   *  Sum: '<S56>/Sum'
   */
  if (11000.0 - DynModel_B.Sum1 > 0.0) {
    rtb_Switch_i = 0.0;
  } else if (11000.0 - DynModel_B.Sum1 < -9000.0) {
    rtb_Switch_i = -9000.0;
  } else {
    rtb_Switch_i = 11000.0 - DynModel_B.Sum1;
  }

  /* Math: '<S56>/Stratosphere Model' incorporates:
   *  Gain: '<S56>/g//R'
   *  Product: '<S56>/Product1'
   *  Saturate: '<S56>/Limit  altitude  to Stratosphere'
   *
   * About '<S56>/Stratosphere Model':
   *  Operator: exp
   */
  rtb_Switch_i = exp(1.0 / DynModel_B.Sum1_e * (0.034163191409533639 *
    rtb_Switch_i));

  /* Product: '<S56>/Product2' incorporates:
   *  Gain: '<S56>/P0'
   */
  DynModel_B.Product2 = 101325.0 * rtb_IntegratorSecondOrder_o1_i * rtb_Switch_i;
  if (rtmIsMajorTimeStep(DynModel_M)) {
    /* Sum: '<S61>/Sum2' incorporates:
     *  Gain: '<S61>/Bar2mBar'
     *  Gain: '<S61>/Pa2Bar'
     *  RandomNumber: '<S61>/Random Number'
     */
    rtb_Saturation_i = 1.0E-5 * DynModel_B.Product2 * 1000.0 +
      DynModel_DW.NextOutput_a;

    /* Saturate: '<S61>/Saturation' */
    if (rtb_Saturation_i > 1200.0) {
      rtb_Saturation_i = 1200.0;
    } else {
      if (rtb_Saturation_i < 10.0) {
        rtb_Saturation_i = 10.0;
      }
    }

    /* End of Saturate: '<S61>/Saturation' */
  }

  /* Product: '<S56>/Product3' incorporates:
   *  Gain: '<S56>/rho0'
   *  Product: '<S56>/Product'
   */
  DynModel_B.Product3 = rtb_IntegratorSecondOrder_o1_i /
    rtb_IntegratorSecondOrder_o2_j * 1.225 * rtb_Switch_i;

  /* Integrator: '<S4>/ub,vb,wb' */
  DynModel_B.ubvbwb[0] = DynModel_X.ubvbwb_CSTATE[0];
  DynModel_B.ubvbwb[1] = DynModel_X.ubvbwb_CSTATE[1];
  DynModel_B.ubvbwb[2] = DynModel_X.ubvbwb_CSTATE[2];
  if (rtmIsMajorTimeStep(DynModel_M)) {
    /* Sum: '<S59>/Sum2' incorporates:
     *  DotProduct: '<S59>/Dot Product'
     *  Gain: '<S59>/Bar2mBar'
     *  Gain: '<S59>/Gain'
     *  Gain: '<S59>/Pa2Bar'
     *  Product: '<S59>/Product'
     *  RandomNumber: '<S59>/Random Number'
     */
    rtb_Saturation_gu = ((DynModel_B.ubvbwb[0] * DynModel_B.ubvbwb[0] +
                          DynModel_B.ubvbwb[1] * DynModel_B.ubvbwb[1]) +
                         DynModel_B.ubvbwb[2] * DynModel_B.ubvbwb[2]) *
      DynModel_B.Product3 * 0.5 * 1.0E-5 * 1000.0 + DynModel_DW.NextOutput_l;

    /* Saturate: '<S59>/Saturation' */
    if (rtb_Saturation_gu > 1000.0) {
      rtb_Saturation_gu = 1000.0;
    } else {
      if (rtb_Saturation_gu < 0.0) {
        rtb_Saturation_gu = 0.0;
      }
    }

    /* End of Saturate: '<S59>/Saturation' */

    /* Sum: '<S51>/Add1' incorporates:
     *  RandomNumber: '<S51>/Random Number'
     */
    rtb_Saturation1 = DynModel_B.Sum1 + DynModel_DW.NextOutput_n;

    /* Saturate: '<S51>/Saturation1' */
    if (!(rtb_Saturation1 >= 0.0)) {
      rtb_Saturation1 = 0.0;
    }

    /* End of Saturate: '<S51>/Saturation1' */

    /* Sum: '<S107>/Sum' incorporates:
     *  Gain: '<S111>/Unit Conversion'
     *  Product: '<S110>/rad lat'
     *  RandomNumber: '<S55>/Random Number'
     *  Sum: '<S55>/Add1'
     */
    rtb_IntegratorSecondOrder_o2_j = (DynModel_DW.NextOutput_o[0] +
      DynModel_B.xeyeze[0]) * 1.5708579706943943E-7 * 57.295779513082323 +
      43.718691;

    /* Switch: '<S115>/Switch' incorporates:
     *  Abs: '<S115>/Abs'
     *  Bias: '<S115>/Bias'
     *  Bias: '<S115>/Bias1'
     *  Constant: '<S115>/Constant2'
     *  Math: '<S115>/Math Function1'
     */
    if (fabs(rtb_IntegratorSecondOrder_o2_j) > 180.0) {
      rtb_IntegratorSecondOrder_o2_j = rt_modd(rtb_IntegratorSecondOrder_o2_j +
        180.0, 360.0) + -180.0;
    }

    /* End of Switch: '<S115>/Switch' */

    /* Abs: '<S112>/Abs1' */
    rtb_Switch_i = fabs(rtb_IntegratorSecondOrder_o2_j);

    /* Switch: '<S112>/Switch' incorporates:
     *  Bias: '<S112>/Bias'
     *  Bias: '<S112>/Bias1'
     *  Constant: '<S108>/Constant'
     *  Constant: '<S108>/Constant1'
     *  Constant: '<S114>/Constant'
     *  Gain: '<S112>/Gain'
     *  Product: '<S112>/Divide1'
     *  RelationalOperator: '<S114>/Compare'
     *  Signum: '<S112>/Sign1'
     *  Switch: '<S108>/Switch1'
     */
    if ((rtb_Switch_i > 90.0) > 0) {
      /* Signum: '<S112>/Sign1' */
      if (rtb_IntegratorSecondOrder_o2_j < 0.0) {
        rtb_IntegratorSecondOrder_o2_j = -1.0;
      } else {
        if (rtb_IntegratorSecondOrder_o2_j > 0.0) {
          rtb_IntegratorSecondOrder_o2_j = 1.0;
        }
      }

      rtb_Switch_fk = (-(rtb_Switch_i + -90.0) + 90.0) *
        rtb_IntegratorSecondOrder_o2_j;
      rtb_Compare_0 = 180;
    } else {
      rtb_Switch_fk = rtb_IntegratorSecondOrder_o2_j;
      rtb_Compare_0 = 0;
    }

    /* End of Switch: '<S112>/Switch' */

    /* Sum: '<S108>/Sum' incorporates:
     *  Gain: '<S111>/Unit Conversion'
     *  Product: '<S110>/rad long '
     *  RandomNumber: '<S55>/Random Number'
     *  Sum: '<S107>/Sum'
     *  Sum: '<S55>/Add1'
     */
    rtb_Sum_j = ((DynModel_DW.NextOutput_o[1] + DynModel_B.xeyeze[1]) *
                 2.1658460268129011E-7 * 57.295779513082323 +
                 DynModel_ConstB.Switch_d) + (real_T)rtb_Compare_0;

    /* Switch: '<S113>/Switch' incorporates:
     *  Abs: '<S113>/Abs'
     *  Bias: '<S113>/Bias'
     *  Bias: '<S113>/Bias1'
     *  Constant: '<S113>/Constant2'
     *  Math: '<S113>/Math Function1'
     */
    if (fabs(rtb_Sum_j) > 180.0) {
      rtb_Sum_j = rt_modd(rtb_Sum_j + 180.0, 360.0) + -180.0;
    }

    /* End of Switch: '<S113>/Switch' */

    /* Sum: '<S107>/Sum1' incorporates:
     *  RandomNumber: '<S55>/Random Number'
     *  Sum: '<S55>/Add1'
     *  UnaryMinus: '<S107>/Ze2height'
     */
    rtb_Sum1_p = -(DynModel_DW.NextOutput_o[2] + DynModel_B.xeyeze[2]);

    /* RandomNumber: '<S55>/Random Number1' */
    rtb_Add2_idx_0 = DynModel_DW.NextOutput_h[0];
    rtb_Add2_idx_1 = DynModel_DW.NextOutput_h[1];
    rtb_Add2_idx_2 = DynModel_DW.NextOutput_h[2];
  }

  /* Integrator: '<S8>/q0 q1 q2 q3' */
  if (DynModel_DW.q0q1q2q3_IWORK.IcNeedsLoading) {
    DynModel_X.q0q1q2q3_CSTATE[0] = DynModel_ConstB.q0;
    DynModel_X.q0q1q2q3_CSTATE[1] = DynModel_ConstB.q1;
    DynModel_X.q0q1q2q3_CSTATE[2] = DynModel_ConstB.q2;
    DynModel_X.q0q1q2q3_CSTATE[3] = DynModel_ConstB.q3;
  }

  rtb_q0q1q2q3[0] = DynModel_X.q0q1q2q3_CSTATE[0];
  rtb_q0q1q2q3[1] = DynModel_X.q0q1q2q3_CSTATE[1];
  rtb_q0q1q2q3[2] = DynModel_X.q0q1q2q3_CSTATE[2];
  rtb_q0q1q2q3[3] = DynModel_X.q0q1q2q3_CSTATE[3];

  /* Sqrt: '<S30>/sqrt' incorporates:
   *  Product: '<S31>/Product'
   *  Product: '<S31>/Product1'
   *  Product: '<S31>/Product2'
   *  Product: '<S31>/Product3'
   *  Sum: '<S31>/Sum'
   */
  rtb_Switch_i = sqrt(((rtb_q0q1q2q3[0] * rtb_q0q1q2q3[0] + rtb_q0q1q2q3[1] *
                        rtb_q0q1q2q3[1]) + rtb_q0q1q2q3[2] * rtb_q0q1q2q3[2]) +
                      rtb_q0q1q2q3[3] * rtb_q0q1q2q3[3]);

  /* Product: '<S29>/Product' */
  rtb_IntegratorSecondOrder_o1_i = rtb_q0q1q2q3[0] / rtb_Switch_i;

  /* Product: '<S29>/Product1' */
  rtb_IntegratorSecondOrder_o2_j = rtb_q0q1q2q3[1] / rtb_Switch_i;

  /* Product: '<S29>/Product2' */
  rtb_Switch_d = rtb_q0q1q2q3[2] / rtb_Switch_i;

  /* Product: '<S29>/Product3' */
  rtb_Switch_i = rtb_q0q1q2q3[3] / rtb_Switch_i;

  /* Sum: '<S19>/Sum' incorporates:
   *  Product: '<S19>/Product'
   *  Product: '<S19>/Product1'
   *  Product: '<S19>/Product2'
   *  Product: '<S19>/Product3'
   */
  DynModel_B.VectorConcatenate[0] = ((rtb_IntegratorSecondOrder_o1_i *
    rtb_IntegratorSecondOrder_o1_i + rtb_IntegratorSecondOrder_o2_j *
    rtb_IntegratorSecondOrder_o2_j) - rtb_Switch_d * rtb_Switch_d) -
    rtb_Switch_i * rtb_Switch_i;

  /* Gain: '<S22>/Gain' incorporates:
   *  Product: '<S22>/Product2'
   *  Product: '<S22>/Product3'
   *  Sum: '<S22>/Sum'
   */
  DynModel_B.VectorConcatenate[1] = (rtb_IntegratorSecondOrder_o2_j *
    rtb_Switch_d - rtb_Switch_i * rtb_IntegratorSecondOrder_o1_i) * 2.0;

  /* Gain: '<S25>/Gain' incorporates:
   *  Product: '<S25>/Product1'
   *  Product: '<S25>/Product2'
   *  Sum: '<S25>/Sum'
   */
  DynModel_B.VectorConcatenate[2] = (rtb_IntegratorSecondOrder_o1_i *
    rtb_Switch_d + rtb_IntegratorSecondOrder_o2_j * rtb_Switch_i) * 2.0;

  /* Gain: '<S20>/Gain' incorporates:
   *  Product: '<S20>/Product2'
   *  Product: '<S20>/Product3'
   *  Sum: '<S20>/Sum'
   */
  DynModel_B.VectorConcatenate[3] = (rtb_Switch_i *
    rtb_IntegratorSecondOrder_o1_i + rtb_IntegratorSecondOrder_o2_j *
    rtb_Switch_d) * 2.0;

  /* Sum: '<S23>/Sum' incorporates:
   *  Product: '<S23>/Product'
   *  Product: '<S23>/Product1'
   *  Product: '<S23>/Product2'
   *  Product: '<S23>/Product3'
   */
  DynModel_B.VectorConcatenate[4] = ((rtb_IntegratorSecondOrder_o1_i *
    rtb_IntegratorSecondOrder_o1_i - rtb_IntegratorSecondOrder_o2_j *
    rtb_IntegratorSecondOrder_o2_j) + rtb_Switch_d * rtb_Switch_d) -
    rtb_Switch_i * rtb_Switch_i;

  /* Gain: '<S26>/Gain' incorporates:
   *  Product: '<S26>/Product1'
   *  Product: '<S26>/Product2'
   *  Sum: '<S26>/Sum'
   */
  DynModel_B.VectorConcatenate[5] = (rtb_Switch_d * rtb_Switch_i -
    rtb_IntegratorSecondOrder_o1_i * rtb_IntegratorSecondOrder_o2_j) * 2.0;

  /* Gain: '<S21>/Gain' incorporates:
   *  Product: '<S21>/Product1'
   *  Product: '<S21>/Product2'
   *  Sum: '<S21>/Sum'
   */
  DynModel_B.VectorConcatenate[6] = (rtb_IntegratorSecondOrder_o2_j *
    rtb_Switch_i - rtb_IntegratorSecondOrder_o1_i * rtb_Switch_d) * 2.0;

  /* Gain: '<S24>/Gain' incorporates:
   *  Product: '<S24>/Product1'
   *  Product: '<S24>/Product2'
   *  Sum: '<S24>/Sum'
   */
  DynModel_B.VectorConcatenate[7] = (rtb_IntegratorSecondOrder_o1_i *
    rtb_IntegratorSecondOrder_o2_j + rtb_Switch_d * rtb_Switch_i) * 2.0;

  /* Sum: '<S27>/Sum' incorporates:
   *  Product: '<S27>/Product'
   *  Product: '<S27>/Product1'
   *  Product: '<S27>/Product2'
   *  Product: '<S27>/Product3'
   */
  DynModel_B.VectorConcatenate[8] = ((rtb_IntegratorSecondOrder_o1_i *
    rtb_IntegratorSecondOrder_o1_i - rtb_IntegratorSecondOrder_o2_j *
    rtb_IntegratorSecondOrder_o2_j) - rtb_Switch_d * rtb_Switch_d) +
    rtb_Switch_i * rtb_Switch_i;

  /* Product: '<S14>/Product' incorporates:
   *  Math: '<S4>/Transpose'
   */
  for (rtb_Compare_0 = 0; rtb_Compare_0 < 3; rtb_Compare_0++) {
    DynModel_B.Product[rtb_Compare_0] = 0.0;
    DynModel_B.Product[rtb_Compare_0] += DynModel_B.VectorConcatenate[3 *
      rtb_Compare_0] * DynModel_B.ubvbwb[0];
    DynModel_B.Product[rtb_Compare_0] += DynModel_B.VectorConcatenate[3 *
      rtb_Compare_0 + 1] * DynModel_B.ubvbwb[1];
    DynModel_B.Product[rtb_Compare_0] += DynModel_B.VectorConcatenate[3 *
      rtb_Compare_0 + 2] * DynModel_B.ubvbwb[2];
  }

  /* End of Product: '<S14>/Product' */
  if (rtmIsMajorTimeStep(DynModel_M)) {
    /* Sum: '<S55>/Add2' */
    rtb_Add2_idx_0 += DynModel_B.Product[0];
    rtb_Add2_idx_1 += DynModel_B.Product[1];
    rtb_Add2_idx_2 += DynModel_B.Product[2];

    /* Trigonometry: '<S55>/Trigonometric Function' */
    rtb_TrigonometricFunction = atan2(rtb_Add2_idx_1, rtb_Add2_idx_0);
  }

  /* Sum: '<S54>/Sum' incorporates:
   *  Gain: '<S92>/Unit Conversion'
   *  Product: '<S91>/rad lat'
   *  Product: '<S91>/x*cos'
   */
  rtb_IntegratorSecondOrder_o2_j = DynModel_B.xeyeze[0] * 1.5708579706943943E-7 *
    57.295779513082323 + 43.718691;

  /* Switch: '<S96>/Switch' incorporates:
   *  Abs: '<S96>/Abs'
   *  Bias: '<S96>/Bias'
   *  Bias: '<S96>/Bias1'
   *  Constant: '<S96>/Constant2'
   *  Math: '<S96>/Math Function1'
   */
  if (fabs(rtb_IntegratorSecondOrder_o2_j) > 180.0) {
    rtb_IntegratorSecondOrder_o2_j = rt_modd(rtb_IntegratorSecondOrder_o2_j +
      180.0, 360.0) + -180.0;
  }

  /* End of Switch: '<S96>/Switch' */

  /* Abs: '<S93>/Abs1' */
  rtb_Switch_i = fabs(rtb_IntegratorSecondOrder_o2_j);

  /* Switch: '<S93>/Switch' incorporates:
   *  Bias: '<S93>/Bias'
   *  Bias: '<S93>/Bias1'
   *  Constant: '<S89>/Constant'
   *  Constant: '<S89>/Constant1'
   *  Constant: '<S95>/Constant'
   *  Gain: '<S93>/Gain'
   *  Product: '<S93>/Divide1'
   *  RelationalOperator: '<S95>/Compare'
   *  Signum: '<S93>/Sign1'
   *  Switch: '<S89>/Switch1'
   */
  if ((rtb_Switch_i > 90.0) > 0) {
    /* Signum: '<S93>/Sign1' */
    if (rtb_IntegratorSecondOrder_o2_j < 0.0) {
      rtb_IntegratorSecondOrder_o2_j = -1.0;
    } else {
      if (rtb_IntegratorSecondOrder_o2_j > 0.0) {
        rtb_IntegratorSecondOrder_o2_j = 1.0;
      }
    }

    rtb_Switch_d = (-(rtb_Switch_i + -90.0) + 90.0) *
      rtb_IntegratorSecondOrder_o2_j;
    rtb_Compare_0 = 180;
  } else {
    rtb_Switch_d = rtb_IntegratorSecondOrder_o2_j;
    rtb_Compare_0 = 0;
  }

  /* End of Switch: '<S93>/Switch' */

  /* Sum: '<S89>/Sum' incorporates:
   *  Gain: '<S92>/Unit Conversion'
   *  Product: '<S91>/rad long '
   *  Product: '<S91>/y*cos'
   *  Sum: '<S54>/Sum'
   */
  rtb_Sum_e = (2.1658460268129011E-7 * DynModel_B.xeyeze[1] * 57.295779513082323
               + DynModel_ConstB.Switch_b) + (real_T)rtb_Compare_0;
  if (rtmIsMajorTimeStep(DynModel_M)) {
    /* Product: '<S58>/Product2' */
    rtb_Sum_h_idx_0 = 0.86602540378443871;
    rtb_Sum_h_idx_1 = 0.49999999999999994;
    rtb_Sum_h_idx_2 = 0.0;

    /* Sum: '<S58>/Sum2' incorporates:
     *  Product: '<S58>/Matrix Multiply2'
     *  RandomNumber: '<S58>/Random Number'
     *  Saturate: '<S58>/Saturation'
     */
    for (rtb_Compare_0 = 0; rtb_Compare_0 < 3; rtb_Compare_0++) {
      rtb_Sum_p[rtb_Compare_0] = (DynModel_B.VectorConcatenate[rtb_Compare_0 + 3]
        * 0.49999999999999994 + DynModel_B.VectorConcatenate[rtb_Compare_0] *
        0.86602540378443871) + DynModel_DW.NextOutput_am;
    }

    /* End of Sum: '<S58>/Sum2' */

    /* Saturate: '<S58>/Saturation' incorporates:
     *  Product: '<S58>/Matrix Multiply2'
     */
    if (rtb_Sum_p[0] > 2.0) {
      rtb_Saturation_g_idx_0 = 2.0;
    } else if (rtb_Sum_p[0] < -2.0) {
      rtb_Saturation_g_idx_0 = -2.0;
    } else {
      rtb_Saturation_g_idx_0 = rtb_Sum_p[0];
    }

    if (rtb_Sum_p[1] > 2.0) {
      rtb_Saturation_g_idx_1 = 2.0;
    } else if (rtb_Sum_p[1] < -2.0) {
      rtb_Saturation_g_idx_1 = -2.0;
    } else {
      rtb_Saturation_g_idx_1 = rtb_Sum_p[1];
    }

    if (rtb_Sum_p[2] > 2.0) {
      rtb_Saturation_g_idx_2 = 2.0;
    } else if (rtb_Sum_p[2] < -2.0) {
      rtb_Saturation_g_idx_2 = -2.0;
    } else {
      rtb_Saturation_g_idx_2 = rtb_Sum_p[2];
    }

    /* Gain: '<S52>/Output' incorporates:
     *  RandomNumber: '<S52>/White Noise'
     */
    DynModel_B.Output = 0.00019364916731037085 * DynModel_DW.NextOutput_lh;
  }

  /* Switch: '<S2>/Switch' incorporates:
   *  Gain: '<S2>/Gain4'
   *  Product: '<S2>/Matrix Multiply'
   */
  if (-DynModel_B.Sum1 >= 0.0) {
    /* Sum: '<S2>/Add8' incorporates:
     *  Gain: '<S2>/Gain3'
     */
    rtb_Switch_i = 5.0 * DynModel_B.Sum1 + -11.772;
    for (rtb_Compare_0 = 0; rtb_Compare_0 < 3; rtb_Compare_0++) {
      rtb_Add5[rtb_Compare_0] = DynModel_B.VectorConcatenate[rtb_Compare_0 + 6] *
        rtb_Switch_i;
    }
  } else {
    rtb_Add5[0] = 0.0;
    rtb_Add5[1] = 0.0;
    rtb_Add5[2] = 0.0;
  }

  /* End of Switch: '<S2>/Switch' */
  if (rtmIsMajorTimeStep(DynModel_M)) {
    /* Memory: '<S2>/Memory2' */
    DynModel_B.Memory2 = DynModel_DW.Memory2_PreviousInput;
  }

  /* Gain: '<S7>/RPM2RADS' incorporates:
   *  Gain: '<S7>/V2RPM'
   *  SecondOrderIntegrator: '<S47>/Integrator, Second-Order'
   *  SecondOrderIntegrator: '<S48>/Integrator, Second-Order'
   *  SecondOrderIntegrator: '<S49>/Integrator, Second-Order'
   *  SecondOrderIntegrator: '<S50>/Integrator, Second-Order'
   */
  DynModel_Y.Rotor_Speed[0] = 950.0 * DynModel_X.IntegratorSecondOrder_CSTATE[0]
    * 0.10471975511965977;
  DynModel_Y.Rotor_Speed[1] = 950.0 * DynModel_X.IntegratorSecondOrder_CSTATE_h
    [0] * 0.10471975511965977;
  DynModel_Y.Rotor_Speed[2] = 950.0 * DynModel_X.IntegratorSecondOrder_CSTATE_n
    [0] * 0.10471975511965977;
  DynModel_Y.Rotor_Speed[3] = 950.0 * DynModel_X.IntegratorSecondOrder_CSTATE_d
    [0] * 0.10471975511965977;

  /* MATLAB Function: '<S6>/multicopter' incorporates:
   *  Constant: '<S6>/h_ref10'
   *  Constant: '<S6>/h_ref12'
   *  Constant: '<S6>/h_ref8'
   */
  /* MATLAB Function 'DynModel/Dynamics/Subsystem/multicopter': '<S44>:1' */
  /* ===============================Parameters================================= */
  /*  */
  /*  a = [ax;ay;az];     Vector with the cross-sectional areas */
  /*  M = [M];            Frame Mass */
  /*  l = [l];            Lenght of the Quadcopter arm */
  /*  Kthr = [Kthr];      Coeff. for the computation of thrust */
  /*  Ktrq = [Ktrq];      Coeff. for the computation of the torque */
  /* ================================Constants================================= */
  /* '<S44>:1:16' */
  /*  coefficient of drag                                       */
  /* ==============================Actuator Mixer============================== */
  /*  [x(roll), y(pitch), z(yaw)] */
  /*  QUAD [X] */
  /* '<S44>:1:23' */
  /*  MIX = [  0   1  -1   ;      % QUAD [+] */
  /*          -1   0   1   ; */
  /*           0  -1  -1   ; */
  /*           1   0   1   ]; */
  /* ==================================Forces==================================  */
  /*  We are evaluating in Body frame     */
  /* '<S44>:1:35' */
  /*  gravitational force */
  /* '<S44>:1:36' */
  /*  drag force */
  /* --------------------------------Thrust Model------------------------------ */
  /* '<S44>:1:39' */
  DynModel_Y.Thursts[0] = DynModel_Y.Rotor_Speed[0] * DynModel_Y.Rotor_Speed[0];
  DynModel_Y.Thursts[1] = DynModel_Y.Rotor_Speed[1] * DynModel_Y.Rotor_Speed[1];
  DynModel_Y.Thursts[2] = DynModel_Y.Rotor_Speed[2] * DynModel_Y.Rotor_Speed[2];
  DynModel_Y.Thursts[3] = DynModel_Y.Rotor_Speed[3] * DynModel_Y.Rotor_Speed[3];
  rtb_Switch_i = 1.2247084269789534E-5 * DynModel_B.Memory2;
  DynModel_Y.Thursts[0] *= rtb_Switch_i;
  DynModel_Y.Thursts[1] *= rtb_Switch_i;
  DynModel_Y.Thursts[2] *= rtb_Switch_i;
  DynModel_Y.Thursts[3] *= rtb_Switch_i;

  /*  rotor thrust */
  /* -------------------------------------------------------------------------- */
  /* '<S44>:1:42' */
  DynModel_Y.Forces[0] = -DynModel_B.ubvbwb[0] * 10.0 * DynModel_B.Memory2 *
    0.016813708498984763 + DynModel_B.VectorConcatenate[6] * 9.81 * 1.2;
  DynModel_Y.Forces[1] = -DynModel_B.ubvbwb[1] * 10.0 * DynModel_B.Memory2 *
    0.018813708498984762 + DynModel_B.VectorConcatenate[7] * 9.81 * 1.2;
  DynModel_Y.Forces[2] = -DynModel_B.ubvbwb[2] * 10.0 * DynModel_B.Memory2 *
    0.18845573684677208 + DynModel_B.VectorConcatenate[8] * 9.81 * 1.2;

  /* '<S44>:1:43' */
  DynModel_Y.Forces[2] -= ((DynModel_Y.Thursts[0] + DynModel_Y.Thursts[1]) +
    DynModel_Y.Thursts[2]) + DynModel_Y.Thursts[3];

  /* ==================================Moments================================= */
  /*  Thrusts contributions to momentum */
  /* '<S44>:1:48' */
  /*  x moment */
  /* '<S44>:1:49' */
  /*  y moment */
  /*  Torques contributions */
  /* momentum_x = sum(abs(MIX(:, 3)) .* rotor_inertia .* rotors) * omega(1);     % x rotor momentum */
  /* momentum_y = sum(abs(MIX(:, 3)) .* rotor_inertia .* rotors) * omega(2);     % y rotor momentum */
  /* momentum_z = sum(MIX(:, 3) .* rotor_inertia .* rotors) * omega(3);          % z rotor momentum */
  /* --------------------------------Torque Model------------------------------ */
  /* '<S44>:1:57' */
  /*  rotor torque */
  /* -------------------------------------------------------------------------- */
  /* '<S44>:1:60' */
  b_y = ((DynModel_Y.Thursts[0] * 0.2 * 1.4142135623730951 / 2.0 +
          -DynModel_Y.Thursts[1] * 0.2 * 1.4142135623730951 / 2.0) +
         -DynModel_Y.Thursts[2] * 0.2 * 1.4142135623730951 / 2.0) +
    DynModel_Y.Thursts[3] * 0.2 * 1.4142135623730951 / 2.0;
  c_y = ((DynModel_Y.Thursts[0] * 0.2 * 1.4142135623730951 / 2.0 +
          DynModel_Y.Thursts[1] * 0.2 * 1.4142135623730951 / 2.0) +
         -DynModel_Y.Thursts[2] * 0.2 * 1.4142135623730951 / 2.0) +
    -DynModel_Y.Thursts[3] * 0.2 * 1.4142135623730951 / 2.0;
  d_y = ((-7.129366502583864E-8 * DynModel_B.Memory2 * (DynModel_Y.Rotor_Speed[0]
           * DynModel_Y.Rotor_Speed[0]) + 7.129366502583864E-8 *
          DynModel_B.Memory2 * (DynModel_Y.Rotor_Speed[1] *
           DynModel_Y.Rotor_Speed[1])) + -7.129366502583864E-8 *
         DynModel_B.Memory2 * (DynModel_Y.Rotor_Speed[2] *
          DynModel_Y.Rotor_Speed[2])) + 7.129366502583864E-8 *
    DynModel_B.Memory2 * (DynModel_Y.Rotor_Speed[3] * DynModel_Y.Rotor_Speed[3]);

  /* Product: '<S4>/Product' incorporates:
   *  Constant: '<S10>/Constant'
   *  Sum: '<S2>/Add'
   */
  /*  - [momentum_x; momentum_y; momentum_z]; */
  /* ========================================================================== */
  DynModel_B.Product_b[0] = (rtb_Add5[0] + DynModel_Y.Forces[0]) / 1.2;
  DynModel_B.Product_b[1] = (rtb_Add5[1] + DynModel_Y.Forces[1]) / 1.2;
  DynModel_B.Product_b[2] = (rtb_Add5[2] + DynModel_Y.Forces[2]) / 1.2;
  if (rtmIsMajorTimeStep(DynModel_M)) {
    /* ZeroOrderHold: '<S127>/Zero-Order Hold1' */
    rtb_Sum_h_idx_0 = DynModel_B.Product_b[0];
    rtb_Sum_h_idx_1 = DynModel_B.Product_b[1];
    rtb_Sum_h_idx_2 = DynModel_B.Product_b[2];
  }

  /* Product: '<S3>/Matrix Multiply1' */
  for (rtb_Compare_0 = 0; rtb_Compare_0 < 3; rtb_Compare_0++) {
    DynModel_B.MatrixMultiply1[rtb_Compare_0] = 0.0;
    DynModel_B.MatrixMultiply1[rtb_Compare_0] +=
      DynModel_B.VectorConcatenate[rtb_Compare_0 + 6] * 9.81;
  }

  /* End of Product: '<S3>/Matrix Multiply1' */
  if (rtmIsMajorTimeStep(DynModel_M)) {
    /* ZeroOrderHold: '<S127>/Zero-Order Hold2' */
    rtb_ZeroOrderHold2_idx_0 = DynModel_B.MatrixMultiply1[0];
    rtb_ZeroOrderHold2_idx_1 = DynModel_B.MatrixMultiply1[1];
    rtb_ZeroOrderHold2_idx_2 = DynModel_B.MatrixMultiply1[2];
  }

  /* Integrator: '<S4>/p,q,r ' */
  DynModel_B.pqr[0] = DynModel_X.pqr_CSTATE[0];
  DynModel_B.pqr[1] = DynModel_X.pqr_CSTATE[1];
  DynModel_B.pqr[2] = DynModel_X.pqr_CSTATE[2];
  if (rtmIsMajorTimeStep(DynModel_M)) {
    /* ZeroOrderHold: '<S127>/Zero-Order Hold' */
    rtb_Product_k4[0] = DynModel_B.pqr[0];
    rtb_Product_k4[1] = DynModel_B.pqr[1];
    rtb_Product_k4[2] = DynModel_B.pqr[2];
  }

  /* Sqrt: '<S33>/sqrt' incorporates:
   *  Product: '<S34>/Product'
   *  Product: '<S34>/Product1'
   *  Product: '<S34>/Product2'
   *  Product: '<S34>/Product3'
   *  Sum: '<S34>/Sum'
   */
  rtb_Switch_i = sqrt(((rtb_q0q1q2q3[0] * rtb_q0q1q2q3[0] + rtb_q0q1q2q3[1] *
                        rtb_q0q1q2q3[1]) + rtb_q0q1q2q3[2] * rtb_q0q1q2q3[2]) +
                      rtb_q0q1q2q3[3] * rtb_q0q1q2q3[3]);

  /* Product: '<S32>/Product' */
  rtb_ixk = rtb_q0q1q2q3[0] / rtb_Switch_i;

  /* Product: '<S32>/Product1' */
  rtb_IntegratorSecondOrder_o2_j = rtb_q0q1q2q3[1] / rtb_Switch_i;

  /* Product: '<S32>/Product2' */
  rtb_IntegratorSecondOrder_o1_i = rtb_q0q1q2q3[2] / rtb_Switch_i;

  /* Product: '<S32>/Product3' */
  rtb_Switch_i = rtb_q0q1q2q3[3] / rtb_Switch_i;

  /* Fcn: '<S16>/fcn1' */
  u0 = (rtb_IntegratorSecondOrder_o2_j * rtb_IntegratorSecondOrder_o1_i +
        rtb_ixk * rtb_Switch_i) * 2.0;

  /* Fcn: '<S16>/fcn2' */
  u1 = ((rtb_ixk * rtb_ixk + rtb_IntegratorSecondOrder_o2_j *
         rtb_IntegratorSecondOrder_o2_j) - rtb_IntegratorSecondOrder_o1_i *
        rtb_IntegratorSecondOrder_o1_i) - rtb_Switch_i * rtb_Switch_i;

  /* Fcn: '<S16>/fcn3' */
  u0_0 = (rtb_IntegratorSecondOrder_o2_j * rtb_Switch_i - rtb_ixk *
          rtb_IntegratorSecondOrder_o1_i) * -2.0;

  /* Fcn: '<S16>/fcn4' */
  rtb_jxk = (rtb_IntegratorSecondOrder_o1_i * rtb_Switch_i + rtb_ixk *
             rtb_IntegratorSecondOrder_o2_j) * 2.0;

  /* Fcn: '<S16>/fcn5' */
  rtb_ixk = ((rtb_ixk * rtb_ixk - rtb_IntegratorSecondOrder_o2_j *
              rtb_IntegratorSecondOrder_o2_j) - rtb_IntegratorSecondOrder_o1_i *
             rtb_IntegratorSecondOrder_o1_i) + rtb_Switch_i * rtb_Switch_i;
  if (rtmIsMajorTimeStep(DynModel_M)) {
    /* Gain: '<S127>/Gain' incorporates:
     *  Constant: '<S127>/wl_ins'
     *  Constant: '<S3>/center of gravity'
     *  Sum: '<S127>/Sum7'
     */
    rtb_Saturation_l_idx_0 = 0.0;
    rtb_Saturation_l_idx_1 = 0.0;
    rtb_Saturation_l_idx_2 = 0.0;

    /* Sum: '<S138>/Sum' incorporates:
     *  Product: '<S140>/i x j'
     *  Product: '<S140>/j x k'
     *  Product: '<S140>/k x i'
     *  Product: '<S141>/i x k'
     *  Product: '<S141>/j x i'
     *  Product: '<S141>/k x j'
     */
    rtb_Saturation_kg[0] = 0.0;
    rtb_Saturation_kg[1] = 0.0;
    rtb_Saturation_kg[2] = 0.0;
  }

  /* Product: '<S35>/Product' */
  for (rtb_Compare_0 = 0; rtb_Compare_0 < 3; rtb_Compare_0++) {
    rtb_Sum_p[rtb_Compare_0] = DynModel_ConstB.Selector[rtb_Compare_0 + 6] *
      DynModel_B.pqr[2] + (DynModel_ConstB.Selector[rtb_Compare_0 + 3] *
      DynModel_B.pqr[1] + DynModel_ConstB.Selector[rtb_Compare_0] *
      DynModel_B.pqr[0]);
  }

  /* End of Product: '<S35>/Product' */

  /* Product: '<S9>/Product2' incorporates:
   *  MATLAB Function: '<S6>/multicopter'
   *  Product: '<S38>/i x j'
   *  Product: '<S38>/j x k'
   *  Product: '<S38>/k x i'
   *  Product: '<S39>/i x k'
   *  Product: '<S39>/j x i'
   *  Product: '<S39>/k x j'
   *  Sum: '<S37>/Sum'
   *  Sum: '<S9>/Sum2'
   */
  rtb_Add6_0[0] = b_y - (DynModel_B.pqr[1] * rtb_Sum_p[2] - DynModel_B.pqr[2] *
    rtb_Sum_p[1]);
  rtb_Add6_0[1] = c_y - (DynModel_B.pqr[2] * rtb_Sum_p[0] - DynModel_B.pqr[0] *
    rtb_Sum_p[2]);
  rtb_Add6_0[2] = d_y - (DynModel_B.pqr[0] * rtb_Sum_p[1] - DynModel_B.pqr[1] *
    rtb_Sum_p[0]);
  rt_mrdivide_U1d1x3_U2d3x3_Yd1x3(rtb_Add6_0, DynModel_ConstB.Selector2,
    DynModel_B.Product2_m);
  if (rtmIsMajorTimeStep(DynModel_M)) {
    /* Sum: '<S127>/Sum' */
    rtb_Sum_h_idx_0 = (rtb_Sum_h_idx_0 - rtb_ZeroOrderHold2_idx_0) +
      rtb_Saturation_kg[0];
    rtb_Sum_h_idx_1 = (rtb_Sum_h_idx_1 - rtb_ZeroOrderHold2_idx_1) +
      rtb_Saturation_kg[1];
    rtb_Switch_i = (rtb_Sum_h_idx_2 - rtb_ZeroOrderHold2_idx_2) +
      rtb_Saturation_kg[2];

    /* Product: '<S127>/Product' incorporates:
     *  Constant: '<S127>/Scale factors & Cross-coupling  errors'
     */
    for (rtb_Compare_0 = 0; rtb_Compare_0 < 3; rtb_Compare_0++) {
      rtb_Saturation_kg[rtb_Compare_0] = DynModel_ConstP.pooled23[rtb_Compare_0
        + 6] * rtb_Switch_i + (DynModel_ConstP.pooled23[rtb_Compare_0 + 3] *
        rtb_Sum_h_idx_1 + DynModel_ConstP.pooled23[rtb_Compare_0] *
        rtb_Sum_h_idx_0);
    }

    /* End of Product: '<S127>/Product' */

    /* Sum: '<S127>/Sum4' */
    rtb_Sum_p[1] = rtb_Saturation_kg[1];
    rtb_Sum_p[2] = rtb_Saturation_kg[2];

    /* Saturate: '<S127>/Saturation' incorporates:
     *  Gain: '<S130>/Output'
     *  RandomNumber: '<S130>/White Noise'
     *  Sum: '<S127>/Sum1'
     *  Sum: '<S127>/Sum4'
     */
    rtb_Switch_i = 0.011180339887498949 * DynModel_DW.NextOutput_k[0] +
      rtb_Saturation_kg[0];
    if (rtb_Switch_i > 19.62) {
      rtb_Saturation_kg[0] = 19.62;
    } else if (rtb_Switch_i < -19.62) {
      rtb_Saturation_kg[0] = -19.62;
    } else {
      rtb_Saturation_kg[0] = rtb_Switch_i;
    }

    rtb_Switch_i = 0.011180339887498949 * DynModel_DW.NextOutput_k[1] +
      rtb_Sum_p[1];
    if (rtb_Switch_i > 19.62) {
      rtb_Saturation_kg[1] = 19.62;
    } else if (rtb_Switch_i < -19.62) {
      rtb_Saturation_kg[1] = -19.62;
    } else {
      rtb_Saturation_kg[1] = rtb_Switch_i;
    }

    rtb_Switch_i = 0.011180339887498949 * DynModel_DW.NextOutput_k[2] +
      rtb_Sum_p[2];
    if (rtb_Switch_i > 19.62) {
      rtb_Saturation_kg[2] = 19.62;
    } else if (rtb_Switch_i < -19.62) {
      rtb_Saturation_kg[2] = -19.62;
    } else {
      rtb_Saturation_kg[2] = rtb_Switch_i;
    }

    /* End of Saturate: '<S127>/Saturation' */

    /* ZeroOrderHold: '<S128>/Zero-Order Hold' */
    rtb_Saturation_l_idx_0 = DynModel_B.pqr[0];
    rtb_Saturation_l_idx_1 = DynModel_B.pqr[1];
    rtb_Saturation_l_idx_2 = DynModel_B.pqr[2];

    /* Product: '<S128>/Product' incorporates:
     *  Constant: '<S128>/Scale factors & Cross-coupling  errors '
     *  ZeroOrderHold: '<S128>/Zero-Order Hold'
     */
    for (rtb_Compare_0 = 0; rtb_Compare_0 < 3; rtb_Compare_0++) {
      rtb_Product_k4[rtb_Compare_0] = DynModel_ConstP.pooled23[rtb_Compare_0 + 6]
        * DynModel_B.pqr[2] + (DynModel_ConstP.pooled23[rtb_Compare_0 + 3] *
        DynModel_B.pqr[1] + DynModel_ConstP.pooled23[rtb_Compare_0] *
        DynModel_B.pqr[0]);
    }

    /* End of Product: '<S128>/Product' */
  }

  /* Gain: '<S126>/Unit Conversion' */
  DynModel_B.UnitConversion[0] = 0.10197162129779283 * DynModel_B.Product_b[0];
  DynModel_B.UnitConversion[1] = 0.10197162129779283 * DynModel_B.Product_b[1];
  DynModel_B.UnitConversion[2] = 0.10197162129779283 * DynModel_B.Product_b[2];
  if (rtmIsMajorTimeStep(DynModel_M)) {
    /* Saturate: '<S128>/Saturation' incorporates:
     *  Gain: '<S147>/Output'
     *  RandomNumber: '<S147>/White Noise'
     *  Sum: '<S128>/Sum1'
     *  Sum: '<S128>/Sum4'
     */
    rtb_Saturation_l_idx_0 = 0.00070710678118654762 * DynModel_DW.NextOutput_p[0]
      + rtb_Product_k4[0];
    if (rtb_Saturation_l_idx_0 > 4.36) {
      rtb_Saturation_l_idx_0 = 4.36;
    } else {
      if (rtb_Saturation_l_idx_0 < -4.36) {
        rtb_Saturation_l_idx_0 = -4.36;
      }
    }

    rtb_Saturation_l_idx_1 = 0.00070710678118654762 * DynModel_DW.NextOutput_p[1]
      + rtb_Product_k4[1];
    if (rtb_Saturation_l_idx_1 > 4.36) {
      rtb_Saturation_l_idx_1 = 4.36;
    } else {
      if (rtb_Saturation_l_idx_1 < -4.36) {
        rtb_Saturation_l_idx_1 = -4.36;
      }
    }

    rtb_Saturation_l_idx_2 = 0.00070710678118654762 * DynModel_DW.NextOutput_p[2]
      + rtb_Product_k4[2];
    if (rtb_Saturation_l_idx_2 > 4.36) {
      rtb_Saturation_l_idx_2 = 4.36;
    } else {
      if (rtb_Saturation_l_idx_2 < -4.36) {
        rtb_Saturation_l_idx_2 = -4.36;
      }
    }

    /* End of Saturate: '<S128>/Saturation' */
  }

  /* Sum: '<S66>/Add' */
  rtb_Switch_i = (DynModel_B.VectorConcatenate[0] +
                  DynModel_B.VectorConcatenate[4]) +
    DynModel_B.VectorConcatenate[8];

  /* If: '<S53>/If' incorporates:
   *  Sum: '<S66>/Add'
   */
  if (rtmIsMajorTimeStep(DynModel_M)) {
    DynModel_DW.If_ActiveSubsystem = (int8_T)!(rtb_Switch_i > 0.0);
  }

  switch (DynModel_DW.If_ActiveSubsystem) {
   case 0L:
    /* Outputs for IfAction SubSystem: '<S53>/Positive Trace' incorporates:
     *  ActionPort: '<S65>/Action Port'
     */
    /* Sqrt: '<S65>/sqrt' incorporates:
     *  Constant: '<S65>/Constant'
     *  Sum: '<S65>/Sum'
     *  Sum: '<S66>/Add'
     */
    rtb_Switch_i = sqrt(rtb_Switch_i + 1.0);

    /* Gain: '<S65>/Gain' */
    DynModel_B.Merge[0] = 0.5 * rtb_Switch_i;

    /* Gain: '<S65>/Gain1' */
    rtb_Switch_i *= 2.0;

    /* Product: '<S65>/Product' incorporates:
     *  Sum: '<S86>/Add'
     *  Sum: '<S87>/Add'
     *  Sum: '<S88>/Add'
     */
    DynModel_B.Merge[1] = (DynModel_B.VectorConcatenate[7] -
      DynModel_B.VectorConcatenate[5]) / rtb_Switch_i;
    DynModel_B.Merge[2] = (DynModel_B.VectorConcatenate[2] -
      DynModel_B.VectorConcatenate[6]) / rtb_Switch_i;
    DynModel_B.Merge[3] = (DynModel_B.VectorConcatenate[3] -
      DynModel_B.VectorConcatenate[1]) / rtb_Switch_i;

    /* End of Outputs for SubSystem: '<S53>/Positive Trace' */
    break;

   case 1L:
    /* Outputs for IfAction SubSystem: '<S53>/Negative Trace' incorporates:
     *  ActionPort: '<S64>/Action Port'
     */
    /* If: '<S64>/Find Maximum Diagonal Value' */
    if (rtmIsMajorTimeStep(DynModel_M)) {
      if ((DynModel_B.VectorConcatenate[4] > DynModel_B.VectorConcatenate[0]) &&
          (DynModel_B.VectorConcatenate[4] > DynModel_B.VectorConcatenate[8])) {
        DynModel_DW.FindMaximumDiagonalValue_ActiveSubsystem = 0;
      } else if (DynModel_B.VectorConcatenate[8] > DynModel_B.VectorConcatenate
                 [0]) {
        DynModel_DW.FindMaximumDiagonalValue_ActiveSubsystem = 1;
      } else {
        DynModel_DW.FindMaximumDiagonalValue_ActiveSubsystem = 2;
      }
    }

    switch (DynModel_DW.FindMaximumDiagonalValue_ActiveSubsystem) {
     case 0L:
      /* Outputs for IfAction SubSystem: '<S64>/Maximum Value at DCM(2,2)' incorporates:
       *  ActionPort: '<S68>/Action Port'
       */
      /* Sqrt: '<S68>/sqrt' incorporates:
       *  Constant: '<S80>/Constant'
       *  Sum: '<S80>/Add'
       */
      rtb_Switch_i = sqrt(((DynModel_B.VectorConcatenate[4] -
                            DynModel_B.VectorConcatenate[0]) -
                           DynModel_B.VectorConcatenate[8]) + 1.0);

      /* Gain: '<S68>/Gain' */
      DynModel_B.Merge[2] = 0.5 * rtb_Switch_i;

      /* Switch: '<S79>/Switch' incorporates:
       *  Constant: '<S79>/Constant1'
       */
      if (rtb_Switch_i != 0.0) {
        rtb_IntegratorSecondOrder_o2_j = 0.5;
      } else {
        rtb_IntegratorSecondOrder_o2_j = 0.0;
        rtb_Switch_i = 1.0;
      }

      /* End of Switch: '<S79>/Switch' */

      /* Product: '<S79>/Product' */
      rtb_Switch_i = rtb_IntegratorSecondOrder_o2_j / rtb_Switch_i;

      /* Gain: '<S68>/Gain1' incorporates:
       *  Product: '<S68>/Product'
       *  Sum: '<S78>/Add'
       */
      DynModel_B.Merge[1] = (DynModel_B.VectorConcatenate[1] +
        DynModel_B.VectorConcatenate[3]) * rtb_Switch_i;

      /* Gain: '<S68>/Gain3' incorporates:
       *  Product: '<S68>/Product'
       *  Sum: '<S77>/Add'
       */
      DynModel_B.Merge[3] = (DynModel_B.VectorConcatenate[5] +
        DynModel_B.VectorConcatenate[7]) * rtb_Switch_i;

      /* Gain: '<S68>/Gain4' incorporates:
       *  Product: '<S68>/Product'
       *  Sum: '<S76>/Add'
       */
      DynModel_B.Merge[0] = (DynModel_B.VectorConcatenate[2] -
        DynModel_B.VectorConcatenate[6]) * rtb_Switch_i;

      /* End of Outputs for SubSystem: '<S64>/Maximum Value at DCM(2,2)' */
      break;

     case 1L:
      /* Outputs for IfAction SubSystem: '<S64>/Maximum Value at DCM(3,3)' incorporates:
       *  ActionPort: '<S69>/Action Port'
       */
      /* Sqrt: '<S69>/sqrt' incorporates:
       *  Constant: '<S85>/Constant'
       *  Sum: '<S85>/Add'
       */
      rtb_Switch_i = sqrt(((DynModel_B.VectorConcatenate[8] -
                            DynModel_B.VectorConcatenate[0]) -
                           DynModel_B.VectorConcatenate[4]) + 1.0);

      /* Gain: '<S69>/Gain' */
      DynModel_B.Merge[3] = 0.5 * rtb_Switch_i;

      /* Switch: '<S84>/Switch' incorporates:
       *  Constant: '<S84>/Constant1'
       */
      if (rtb_Switch_i != 0.0) {
        rtb_IntegratorSecondOrder_o2_j = 0.5;
      } else {
        rtb_IntegratorSecondOrder_o2_j = 0.0;
        rtb_Switch_i = 1.0;
      }

      /* End of Switch: '<S84>/Switch' */

      /* Product: '<S84>/Product' */
      rtb_Switch_i = rtb_IntegratorSecondOrder_o2_j / rtb_Switch_i;

      /* Gain: '<S69>/Gain1' incorporates:
       *  Product: '<S69>/Product'
       *  Sum: '<S81>/Add'
       */
      DynModel_B.Merge[1] = (DynModel_B.VectorConcatenate[2] +
        DynModel_B.VectorConcatenate[6]) * rtb_Switch_i;

      /* Gain: '<S69>/Gain2' incorporates:
       *  Product: '<S69>/Product'
       *  Sum: '<S82>/Add'
       */
      DynModel_B.Merge[2] = (DynModel_B.VectorConcatenate[5] +
        DynModel_B.VectorConcatenate[7]) * rtb_Switch_i;

      /* Gain: '<S69>/Gain3' incorporates:
       *  Product: '<S69>/Product'
       *  Sum: '<S83>/Add'
       */
      DynModel_B.Merge[0] = (DynModel_B.VectorConcatenate[3] -
        DynModel_B.VectorConcatenate[1]) * rtb_Switch_i;

      /* End of Outputs for SubSystem: '<S64>/Maximum Value at DCM(3,3)' */
      break;

     case 2L:
      /* Outputs for IfAction SubSystem: '<S64>/Maximum Value at DCM(1,1)' incorporates:
       *  ActionPort: '<S67>/Action Port'
       */
      /* Sqrt: '<S67>/sqrt' incorporates:
       *  Constant: '<S75>/Constant'
       *  Sum: '<S75>/Add'
       */
      rtb_Switch_i = sqrt(((DynModel_B.VectorConcatenate[0] -
                            DynModel_B.VectorConcatenate[4]) -
                           DynModel_B.VectorConcatenate[8]) + 1.0);

      /* Gain: '<S67>/Gain' */
      DynModel_B.Merge[1] = 0.5 * rtb_Switch_i;

      /* Switch: '<S74>/Switch' incorporates:
       *  Constant: '<S74>/Constant1'
       */
      if (rtb_Switch_i != 0.0) {
        rtb_IntegratorSecondOrder_o2_j = 0.5;
      } else {
        rtb_IntegratorSecondOrder_o2_j = 0.0;
        rtb_Switch_i = 1.0;
      }

      /* End of Switch: '<S74>/Switch' */

      /* Product: '<S74>/Product' */
      rtb_Switch_i = rtb_IntegratorSecondOrder_o2_j / rtb_Switch_i;

      /* Gain: '<S67>/Gain1' incorporates:
       *  Product: '<S67>/Product'
       *  Sum: '<S73>/Add'
       */
      DynModel_B.Merge[2] = (DynModel_B.VectorConcatenate[1] +
        DynModel_B.VectorConcatenate[3]) * rtb_Switch_i;

      /* Gain: '<S67>/Gain2' incorporates:
       *  Product: '<S67>/Product'
       *  Sum: '<S71>/Add'
       */
      DynModel_B.Merge[3] = (DynModel_B.VectorConcatenate[2] +
        DynModel_B.VectorConcatenate[6]) * rtb_Switch_i;

      /* Gain: '<S67>/Gain3' incorporates:
       *  Product: '<S67>/Product'
       *  Sum: '<S72>/Add'
       */
      DynModel_B.Merge[0] = (DynModel_B.VectorConcatenate[7] -
        DynModel_B.VectorConcatenate[5]) * rtb_Switch_i;

      /* End of Outputs for SubSystem: '<S64>/Maximum Value at DCM(1,1)' */
      break;
    }

    /* End of If: '<S64>/Find Maximum Diagonal Value' */
    /* End of Outputs for SubSystem: '<S53>/Negative Trace' */
    break;
  }

  /* End of If: '<S53>/If' */

  /* DotProduct: '<S18>/Dot Product' */
  rtb_Switch_i = ((rtb_q0q1q2q3[0] * rtb_q0q1q2q3[0] + rtb_q0q1q2q3[1] *
                   rtb_q0q1q2q3[1]) + rtb_q0q1q2q3[2] * rtb_q0q1q2q3[2]) +
    rtb_q0q1q2q3[3] * rtb_q0q1q2q3[3];

  /* Fcn: '<S18>/q0dot' incorporates:
   *  Constant: '<S18>/Constant'
   *  DotProduct: '<S18>/Dot Product'
   *  Sum: '<S18>/Sum'
   */
  DynModel_B.q0dot = ((rtb_q0q1q2q3[1] * DynModel_B.pqr[0] + rtb_q0q1q2q3[2] *
                       DynModel_B.pqr[1]) + rtb_q0q1q2q3[3] * DynModel_B.pqr[2])
    * -0.5 + (1.0 - rtb_Switch_i) * rtb_q0q1q2q3[0];

  /* Fcn: '<S18>/q1dot' incorporates:
   *  Constant: '<S18>/Constant'
   *  DotProduct: '<S18>/Dot Product'
   *  Sum: '<S18>/Sum'
   */
  DynModel_B.q1dot = ((rtb_q0q1q2q3[0] * DynModel_B.pqr[0] + rtb_q0q1q2q3[2] *
                       DynModel_B.pqr[2]) - rtb_q0q1q2q3[3] * DynModel_B.pqr[1])
    * 0.5 + (1.0 - rtb_Switch_i) * rtb_q0q1q2q3[1];

  /* Fcn: '<S18>/q2dot' incorporates:
   *  Constant: '<S18>/Constant'
   *  DotProduct: '<S18>/Dot Product'
   *  Sum: '<S18>/Sum'
   */
  DynModel_B.q2dot = ((rtb_q0q1q2q3[0] * DynModel_B.pqr[1] + rtb_q0q1q2q3[3] *
                       DynModel_B.pqr[0]) - rtb_q0q1q2q3[1] * DynModel_B.pqr[2])
    * 0.5 + (1.0 - rtb_Switch_i) * rtb_q0q1q2q3[2];

  /* Fcn: '<S18>/q3dot' incorporates:
   *  Constant: '<S18>/Constant'
   *  DotProduct: '<S18>/Dot Product'
   *  Sum: '<S18>/Sum'
   */
  DynModel_B.q3dot = ((rtb_q0q1q2q3[0] * DynModel_B.pqr[2] + rtb_q0q1q2q3[1] *
                       DynModel_B.pqr[1]) - rtb_q0q1q2q3[2] * DynModel_B.pqr[0])
    * 0.5 + (1.0 - rtb_Switch_i) * rtb_q0q1q2q3[3];

  /* Sum: '<S4>/Sum' incorporates:
   *  Product: '<S40>/i x j'
   *  Product: '<S40>/j x k'
   *  Product: '<S40>/k x i'
   *  Product: '<S41>/i x k'
   *  Product: '<S41>/j x i'
   *  Product: '<S41>/k x j'
   *  Sum: '<S11>/Sum'
   */
  DynModel_B.Sum[0] = (DynModel_B.ubvbwb[1] * DynModel_B.pqr[2] -
                       DynModel_B.ubvbwb[2] * DynModel_B.pqr[1]) +
    DynModel_B.Product_b[0];
  DynModel_B.Sum[1] = (DynModel_B.ubvbwb[2] * DynModel_B.pqr[0] -
                       DynModel_B.ubvbwb[0] * DynModel_B.pqr[2]) +
    DynModel_B.Product_b[1];
  DynModel_B.Sum[2] = (DynModel_B.ubvbwb[0] * DynModel_B.pqr[1] -
                       DynModel_B.ubvbwb[1] * DynModel_B.pqr[0]) +
    DynModel_B.Product_b[2];

  /* Saturate: '<S2>/Saturation' incorporates:
   *  Inport: '<Root>/PWM1'
   */
  if (DynModel_U.PWM1 > 1.0) {
    DynModel_B.Saturation = 1.0;
  } else if (DynModel_U.PWM1 < 0.0) {
    DynModel_B.Saturation = 0.0;
  } else {
    DynModel_B.Saturation = DynModel_U.PWM1;
  }

  /* End of Saturate: '<S2>/Saturation' */

  /* Saturate: '<S2>/Saturation1' incorporates:
   *  Inport: '<Root>/PWM2'
   */
  if (DynModel_U.PWM2 > 1.0) {
    DynModel_B.Saturation1 = 1.0;
  } else if (DynModel_U.PWM2 < 0.0) {
    DynModel_B.Saturation1 = 0.0;
  } else {
    DynModel_B.Saturation1 = DynModel_U.PWM2;
  }

  /* End of Saturate: '<S2>/Saturation1' */

  /* Saturate: '<S2>/Saturation2' incorporates:
   *  Inport: '<Root>/PWM3'
   */
  if (DynModel_U.PWM3 > 1.0) {
    DynModel_B.Saturation2 = 1.0;
  } else if (DynModel_U.PWM3 < 0.0) {
    DynModel_B.Saturation2 = 0.0;
  } else {
    DynModel_B.Saturation2 = DynModel_U.PWM3;
  }

  /* End of Saturate: '<S2>/Saturation2' */

  /* Saturate: '<S2>/Saturation3' incorporates:
   *  Inport: '<Root>/PWM4'
   */
  if (DynModel_U.PWM4 > 1.0) {
    DynModel_B.Saturation3 = 1.0;
  } else if (DynModel_U.PWM4 < 0.0) {
    DynModel_B.Saturation3 = 0.0;
  } else {
    DynModel_B.Saturation3 = DynModel_U.PWM4;
  }

  /* End of Saturate: '<S2>/Saturation3' */
  if (rtmIsMajorTimeStep(DynModel_M)) {
    /* SignalConversion: '<S5>/TmpSignal ConversionAt SFunction Inport1' incorporates:
     *  MATLAB Function: '<S2>/LiPo Battery'
     */
    DynModel_B.voltage[0] = DynModel_B.Saturation;
    DynModel_B.voltage[1] = DynModel_B.Saturation1;
    DynModel_B.voltage[2] = DynModel_B.Saturation2;
    DynModel_B.voltage[3] = DynModel_B.Saturation3;

    /* MATLAB Function: '<S2>/LiPo Battery' */
    /* MATLAB Function 'DynModel/Dynamics/LiPo Battery': '<S5>:1' */
    /* ================================Constants================================= */
    /* '<S5>:1:6' */
    /*  maximum current per motor             (A) */
    /* '<S5>:1:7' */
    /*  current of auxillary components       (A)         */
    /*  LiPo current capacity                 (Ah)                                  */
    /*  LiPo cell count */
    /* ========================================================================== */
    /* '<S5>:1:19' */
    /* '<S5>:1:20' */
    DynModel_DW.discharge += ((((DynModel_B.voltage[0] * DynModel_B.voltage[0] *
      5.5 + DynModel_B.voltage[1] * DynModel_B.voltage[1] * 5.5) +
      DynModel_B.voltage[2] * DynModel_B.voltage[2] * 5.5) + DynModel_B.voltage
      [3] * DynModel_B.voltage[3] * 5.5) + 0.1) * 1.111111111111111E-6;

    /* '<S5>:1:22' */
    if ((0.0 < DynModel_DW.discharge) && (DynModel_DW.discharge <= 0.2)) {
      /* '<S5>:1:24' */
      /* '<S5>:1:25' */
      rtb_Switch_i = ((DynModel_DW.discharge * DynModel_DW.discharge * 16.975 +
                       -14.029 * pow(DynModel_DW.discharge, 3.0)) - 5.3339 *
                      DynModel_DW.discharge) + 4.2;
    } else if ((0.2 < DynModel_DW.discharge) && (DynModel_DW.discharge < 0.7)) {
      /* '<S5>:1:26' */
      /* '<S5>:1:27' */
      rtb_Switch_i = -0.2 * DynModel_DW.discharge + 3.74;
    } else {
      /* '<S5>:1:29' */
      rtb_Switch_i = ((DynModel_DW.discharge * DynModel_DW.discharge * 89.6 +
                       -48.0 * pow(DynModel_DW.discharge, 3.0)) - 55.08 *
                      DynModel_DW.discharge) + 14.716;
    }

    if (rtb_Switch_i < 2.5) {
      /* '<S5>:1:32' */
      /* '<S5>:1:33' */
      rtb_Switch_i = 0.0;
    }

    /* '<S5>:1:36' */
    rtb_Switch_i *= 2.0;

    /* '<S5>:1:37' */
    DynModel_B.voltage[0] *= rtb_Switch_i;
    DynModel_B.voltage[1] *= rtb_Switch_i;
    DynModel_B.voltage[2] *= rtb_Switch_i;
    DynModel_B.voltage[3] *= rtb_Switch_i;

    /* DataTypeConversion: '<S1>/Data Type Conversion14' */
    /* ========================================================================== */
    rtb_DataTypeConversion14_idx_0 = (real32_T)rtb_Saturation_kg[0];
    rtb_DataTypeConversion14_idx_1 = (real32_T)rtb_Saturation_kg[1];
    rtb_DataTypeConversion14_idx_2 = (real32_T)rtb_Saturation_kg[2];

    /* DataTypeConversion: '<S1>/Data Type Conversion15' */
    rtb_DataTypeConversion15_idx_0 = (real32_T)rtb_Saturation_l_idx_0;
    rtb_DataTypeConversion15_idx_1 = (real32_T)rtb_Saturation_l_idx_1;
    rtb_DataTypeConversion15_idx_2 = (real32_T)rtb_Saturation_l_idx_2;

    /* DataTypeConversion: '<S1>/Data Type Conversion12' */
    rtb_DataTypeConversion12_idx_0 = (real32_T)rtb_Saturation_g_idx_0;
    rtb_DataTypeConversion12_idx_1 = (real32_T)rtb_Saturation_g_idx_1;
    rtb_DataTypeConversion12_idx_2 = (real32_T)rtb_Saturation_g_idx_2;

    /* Outport: '<Root>/Temp' incorporates:
     *  DataTypeConversion: '<S1>/Data Type Conversion'
     */
    DynModel_Y.Temp = (real32_T)rtb_Saturation;

    /* Outport: '<Root>/Press' incorporates:
     *  DataTypeConversion: '<S1>/Data Type Conversion1'
     */
    DynModel_Y.Press = (real32_T)rtb_Saturation_i;

    /* Outport: '<Root>/diff_Pres' incorporates:
     *  DataTypeConversion: '<S1>/Data Type Conversion3'
     */
    DynModel_Y.diff_Pres = (real32_T)rtb_Saturation_gu;

    /* Outport: '<Root>/Baro_Alt' incorporates:
     *  DataTypeConversion: '<S1>/Data Type Conversion4'
     */
    DynModel_Y.Baro_Alt = (real32_T)rtb_Saturation1;

    /* Outport: '<Root>/Gps_Lat' incorporates:
     *  DataTypeConversion: '<S1>/Data Type Conversion5'
     */
    DynModel_Y.Gps_Lat = (real32_T)rtb_Switch_fk;

    /* Outport: '<Root>/Gps_Lon' incorporates:
     *  DataTypeConversion: '<S1>/Data Type Conversion6'
     */
    DynModel_Y.Gps_Lon = (real32_T)rtb_Sum_j;

    /* Outport: '<Root>/Gps_Alt' incorporates:
     *  DataTypeConversion: '<S1>/Data Type Conversion7'
     */
    DynModel_Y.Gps_Alt = (real32_T)rtb_Sum1_p;

    /* Outport: '<Root>/Gps_V' incorporates:
     *  DataTypeConversion: '<S1>/Data Type Conversion8'
     */
    DynModel_Y.Gps_V[0] = (real32_T)rtb_Add2_idx_0;
    DynModel_Y.Gps_V[1] = (real32_T)rtb_Add2_idx_1;
    DynModel_Y.Gps_V[2] = (real32_T)rtb_Add2_idx_2;

    /* Outport: '<Root>/Gps_V_Mod' incorporates:
     *  DataTypeConversion: '<S1>/Data Type Conversion18'
     *  DotProduct: '<S1>/Dot Product'
     *  Sqrt: '<S1>/Sqrt'
     */
    DynModel_Y.Gps_V_Mod = (real32_T)sqrt((rtb_Add2_idx_0 * rtb_Add2_idx_0 +
      rtb_Add2_idx_1 * rtb_Add2_idx_1) + rtb_Add2_idx_2 * rtb_Add2_idx_2);

    /* Outport: '<Root>/COG' incorporates:
     *  DataTypeConversion: '<S1>/Data Type Conversion9'
     */
    DynModel_Y.COG = (real32_T)rtb_TrigonometricFunction;
  }

  /* Outport: '<Root>/Lat_Lon_Alt' incorporates:
   *  DataTypeConversion: '<S1>/Data Type Conversion10'
   */
  DynModel_Y.Lat_Lon_Alt[0] = (real32_T)rtb_Switch_d;

  /* Switch: '<S94>/Switch' incorporates:
   *  Abs: '<S94>/Abs'
   */
  if (fabs(rtb_Sum_e) > 180.0) {
    /* Outport: '<Root>/Lat_Lon_Alt' incorporates:
     *  Bias: '<S94>/Bias'
     *  Bias: '<S94>/Bias1'
     *  Constant: '<S94>/Constant2'
     *  DataTypeConversion: '<S1>/Data Type Conversion10'
     *  Math: '<S94>/Math Function1'
     */
    DynModel_Y.Lat_Lon_Alt[1] = (real32_T)(rt_modd(rtb_Sum_e + 180.0, 360.0) +
      -180.0);
  } else {
    /* Outport: '<Root>/Lat_Lon_Alt' incorporates:
     *  DataTypeConversion: '<S1>/Data Type Conversion10'
     */
    DynModel_Y.Lat_Lon_Alt[1] = (real32_T)rtb_Sum_e;
  }

  /* End of Switch: '<S94>/Switch' */

  /* Outport: '<Root>/Lat_Lon_Alt' incorporates:
   *  DataTypeConversion: '<S1>/Data Type Conversion11'
   */
  DynModel_Y.Lat_Lon_Alt[2] = (real32_T)DynModel_B.Sum1;
  if (rtmIsMajorTimeStep(DynModel_M)) {
    /* Outport: '<Root>/Magn' */
    DynModel_Y.Magn[0] = rtb_DataTypeConversion12_idx_0;
    DynModel_Y.Magn[1] = rtb_DataTypeConversion12_idx_1;
    DynModel_Y.Magn[2] = rtb_DataTypeConversion12_idx_2;
  }

  /* Outport: '<Root>/RPY' incorporates:
   *  DataTypeConversion: '<S1>/Data Type Conversion13'
   *  Sum: '<S3>/Add5'
   *  Trigonometry: '<S16>/Trigonometric Function3'
   */
  DynModel_Y.RPY[0] = (real32_T)(atan2(rtb_jxk, rtb_ixk) + DynModel_B.Output);

  /* Trigonometry: '<S16>/trigFcn' */
  if (u0_0 > 1.0) {
    u0_0 = 1.0;
  } else {
    if (u0_0 < -1.0) {
      u0_0 = -1.0;
    }
  }

  /* Outport: '<Root>/RPY' incorporates:
   *  DataTypeConversion: '<S1>/Data Type Conversion13'
   *  Sum: '<S3>/Add5'
   *  Trigonometry: '<S16>/Trigonometric Function1'
   *  Trigonometry: '<S16>/trigFcn'
   */
  DynModel_Y.RPY[1] = (real32_T)(asin(u0_0) + DynModel_B.Output);
  DynModel_Y.RPY[2] = (real32_T)(atan2(u0, u1) + DynModel_B.Output);
  if (rtmIsMajorTimeStep(DynModel_M)) {
    /* Outport: '<Root>/Accelerometer' */
    DynModel_Y.Accelerometer[0] = rtb_DataTypeConversion14_idx_0;
    DynModel_Y.Accelerometer[1] = rtb_DataTypeConversion14_idx_1;
    DynModel_Y.Accelerometer[2] = rtb_DataTypeConversion14_idx_2;

    /* Outport: '<Root>/Gyro' */
    DynModel_Y.Gyro[0] = rtb_DataTypeConversion15_idx_0;
    DynModel_Y.Gyro[1] = rtb_DataTypeConversion15_idx_1;
    DynModel_Y.Gyro[2] = rtb_DataTypeConversion15_idx_2;
  }

  /* Outport: '<Root>/Quaternion' incorporates:
   *  DataTypeConversion: '<S1>/Data Type Conversion16'
   */
  DynModel_Y.Quaternion[0] = (real32_T)DynModel_B.Merge[0];
  DynModel_Y.Quaternion[1] = (real32_T)DynModel_B.Merge[1];
  DynModel_Y.Quaternion[2] = (real32_T)DynModel_B.Merge[2];
  DynModel_Y.Quaternion[3] = (real32_T)DynModel_B.Merge[3];

  /* Outport: '<Root>/Torques' incorporates:
   *  MATLAB Function: '<S6>/multicopter'
   */
  DynModel_Y.Torques[0] = b_y;
  DynModel_Y.Torques[1] = c_y;
  DynModel_Y.Torques[2] = d_y;

  /* Sum: '<S47>/Sum2' incorporates:
   *  Gain: '<S47>/2*zeta * wn'
   *  Gain: '<S47>/wn^2'
   *  Product: '<S7>/Product'
   *  SecondOrderIntegrator: '<S47>/Integrator, Second-Order'
   *  Sum: '<S47>/Sum3'
   */
  DynModel_B.Sum2 = (DynModel_B.Saturation * DynModel_B.voltage[0] -
                     DynModel_X.IntegratorSecondOrder_CSTATE[0]) * 4900.0 -
    140.0 * DynModel_X.IntegratorSecondOrder_CSTATE[1];

  /* Sum: '<S48>/Sum2' incorporates:
   *  Gain: '<S48>/2*zeta * wn'
   *  Gain: '<S48>/wn^2'
   *  Product: '<S7>/Product'
   *  SecondOrderIntegrator: '<S48>/Integrator, Second-Order'
   *  Sum: '<S48>/Sum3'
   */
  DynModel_B.Sum2_j = (DynModel_B.Saturation1 * DynModel_B.voltage[1] -
                       DynModel_X.IntegratorSecondOrder_CSTATE_h[0]) * 4900.0 -
    140.0 * DynModel_X.IntegratorSecondOrder_CSTATE_h[1];

  /* Sum: '<S49>/Sum2' incorporates:
   *  Gain: '<S49>/2*zeta * wn'
   *  Gain: '<S49>/wn^2'
   *  Product: '<S7>/Product'
   *  SecondOrderIntegrator: '<S49>/Integrator, Second-Order'
   *  Sum: '<S49>/Sum3'
   */
  DynModel_B.Sum2_c = (DynModel_B.Saturation2 * DynModel_B.voltage[2] -
                       DynModel_X.IntegratorSecondOrder_CSTATE_n[0]) * 4900.0 -
    140.0 * DynModel_X.IntegratorSecondOrder_CSTATE_n[1];

  /* Sum: '<S50>/Sum2' incorporates:
   *  Gain: '<S50>/2*zeta * wn'
   *  Gain: '<S50>/wn^2'
   *  Product: '<S7>/Product'
   *  SecondOrderIntegrator: '<S50>/Integrator, Second-Order'
   *  Sum: '<S50>/Sum3'
   */
  DynModel_B.Sum2_p = (DynModel_B.Saturation3 * DynModel_B.voltage[3] -
                       DynModel_X.IntegratorSecondOrder_CSTATE_d[0]) * 4900.0 -
    140.0 * DynModel_X.IntegratorSecondOrder_CSTATE_d[1];
  if (rtmIsMajorTimeStep(DynModel_M)) {
    if (rtmIsMajorTimeStep(DynModel_M)) {
      /* Update for RandomNumber: '<S63>/Random Number' */
      DynModel_DW.NextOutput = rt_nrand_Upu32_Yd_f_pw(&DynModel_DW.RandSeed) *
        0.01;

      /* Update for RandomNumber: '<S61>/Random Number' */
      DynModel_DW.NextOutput_a = rt_nrand_Upu32_Yd_f_pw(&DynModel_DW.RandSeed_f)
        * 0.0031622776601683794;

      /* Update for RandomNumber: '<S59>/Random Number' */
      DynModel_DW.NextOutput_l = rt_nrand_Upu32_Yd_f_pw(&DynModel_DW.RandSeed_fw)
        * 0.0031622776601683794;

      /* Update for RandomNumber: '<S51>/Random Number' */
      DynModel_DW.NextOutput_n = rt_nrand_Upu32_Yd_f_pw(&DynModel_DW.RandSeed_fm)
        * 0.01;

      /* Update for RandomNumber: '<S55>/Random Number' */
      DynModel_DW.NextOutput_o[0] = rt_nrand_Upu32_Yd_f_pw
        (&DynModel_DW.RandSeed_e[0]) * 15.0;
      DynModel_DW.NextOutput_o[1] = rt_nrand_Upu32_Yd_f_pw
        (&DynModel_DW.RandSeed_e[1]) * 15.0;
      DynModel_DW.NextOutput_o[2] = rt_nrand_Upu32_Yd_f_pw
        (&DynModel_DW.RandSeed_e[2]) * 20.0;

      /* Update for RandomNumber: '<S55>/Random Number1' */
      DynModel_DW.NextOutput_h[0] = rt_nrand_Upu32_Yd_f_pw
        (&DynModel_DW.RandSeed_i[0]) * 0.031622776601683791;
      DynModel_DW.NextOutput_h[1] = rt_nrand_Upu32_Yd_f_pw
        (&DynModel_DW.RandSeed_i[1]) * 0.031622776601683791;
      DynModel_DW.NextOutput_h[2] = rt_nrand_Upu32_Yd_f_pw
        (&DynModel_DW.RandSeed_i[2]) * 0.031622776601683791;
    }

    /* Update for Integrator: '<S8>/q0 q1 q2 q3' */
    DynModel_DW.q0q1q2q3_IWORK.IcNeedsLoading = 0;
    if (rtmIsMajorTimeStep(DynModel_M)) {
      /* Update for RandomNumber: '<S58>/Random Number' */
      DynModel_DW.NextOutput_am = rt_nrand_Upu32_Yd_f_pw(&DynModel_DW.RandSeed_p)
        * 0.0031622776601683794;

      /* Update for RandomNumber: '<S52>/White Noise' */
      DynModel_DW.NextOutput_lh = rt_nrand_Upu32_Yd_f_pw(&DynModel_DW.RandSeed_l);

      /* Update for Memory: '<S2>/Memory2' */
      DynModel_DW.Memory2_PreviousInput = DynModel_B.Product3;

      /* Update for RandomNumber: '<S130>/White Noise' */
      DynModel_DW.NextOutput_k[0] = rt_nrand_Upu32_Yd_f_pw
        (&DynModel_DW.RandSeed_ls[0]);
      DynModel_DW.NextOutput_k[1] = rt_nrand_Upu32_Yd_f_pw
        (&DynModel_DW.RandSeed_ls[1]);
      DynModel_DW.NextOutput_k[2] = rt_nrand_Upu32_Yd_f_pw
        (&DynModel_DW.RandSeed_ls[2]);

      /* Update for RandomNumber: '<S147>/White Noise' */
      DynModel_DW.NextOutput_p[0] = rt_nrand_Upu32_Yd_f_pw
        (&DynModel_DW.RandSeed_j[0]);
      DynModel_DW.NextOutput_p[1] = rt_nrand_Upu32_Yd_f_pw
        (&DynModel_DW.RandSeed_j[1]);
      DynModel_DW.NextOutput_p[2] = rt_nrand_Upu32_Yd_f_pw
        (&DynModel_DW.RandSeed_j[2]);
    }
  }                                    /* end MajorTimeStep */

  if (rtmIsMajorTimeStep(DynModel_M)) {
    rt_ertODEUpdateContinuousStates(&DynModel_M->solverInfo);

    /* Update absolute time for base rate */
    /* The "clockTick0" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick0"
     * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick0 and the high bits
     * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++DynModel_M->Timing.clockTick0)) {
      ++DynModel_M->Timing.clockTickH0;
    }

    DynModel_M->Timing.t[0] = rtsiGetSolverStopTime(&DynModel_M->solverInfo);

    {
      /* Update absolute timer for sample time: [0.004s, 0.0s] */
      /* The "clockTick1" counts the number of times the code of this task has
       * been executed. The resolution of this integer timer is 0.004, which is the step size
       * of the task. Size of "clockTick1" ensures timer will not overflow during the
       * application lifespan selected.
       * Timer of this task consists of two 32 bit unsigned integers.
       * The two integers represent the low bits Timing.clockTick1 and the high bits
       * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
       */
      DynModel_M->Timing.clockTick1++;
      if (!DynModel_M->Timing.clockTick1) {
        DynModel_M->Timing.clockTickH1++;
      }
    }
  }                                    /* end MajorTimeStep */
}

/* Derivatives for root system: '<Root>' */
void DynModel_derivatives(void)
{
  XDot_DynModel_T *_rtXdot;
  _rtXdot = ((XDot_DynModel_T *) DynModel_M->ModelData.derivs);

  /* Derivatives for Integrator: '<S4>/xe,ye,ze' */
  _rtXdot->xeyeze_CSTATE[0] = DynModel_B.Product[0];
  _rtXdot->xeyeze_CSTATE[1] = DynModel_B.Product[1];
  _rtXdot->xeyeze_CSTATE[2] = DynModel_B.Product[2];

  /* Derivatives for Integrator: '<S4>/ub,vb,wb' */
  _rtXdot->ubvbwb_CSTATE[0] = DynModel_B.Sum[0];
  _rtXdot->ubvbwb_CSTATE[1] = DynModel_B.Sum[1];
  _rtXdot->ubvbwb_CSTATE[2] = DynModel_B.Sum[2];

  /* Derivatives for Integrator: '<S8>/q0 q1 q2 q3' */
  {
    ((XDot_DynModel_T *) DynModel_M->ModelData.derivs)->q0q1q2q3_CSTATE[0] =
      DynModel_B.q0dot;
    ((XDot_DynModel_T *) DynModel_M->ModelData.derivs)->q0q1q2q3_CSTATE[1] =
      DynModel_B.q1dot;
    ((XDot_DynModel_T *) DynModel_M->ModelData.derivs)->q0q1q2q3_CSTATE[2] =
      DynModel_B.q2dot;
    ((XDot_DynModel_T *) DynModel_M->ModelData.derivs)->q0q1q2q3_CSTATE[3] =
      DynModel_B.q3dot;
  }

  /* Derivatives for SecondOrderIntegrator: '<S47>/Integrator, Second-Order' */
  if (DynModel_DW.IntegratorSecondOrder_MODE == 0) {
    _rtXdot->IntegratorSecondOrder_CSTATE[0] =
      DynModel_X.IntegratorSecondOrder_CSTATE[1];
    _rtXdot->IntegratorSecondOrder_CSTATE[1] = DynModel_B.Sum2;
  }

  /* End of Derivatives for SecondOrderIntegrator: '<S47>/Integrator, Second-Order' */

  /* Derivatives for SecondOrderIntegrator: '<S48>/Integrator, Second-Order' */
  if (DynModel_DW.IntegratorSecondOrder_MODE_p == 0) {
    _rtXdot->IntegratorSecondOrder_CSTATE_h[0] =
      DynModel_X.IntegratorSecondOrder_CSTATE_h[1];
    _rtXdot->IntegratorSecondOrder_CSTATE_h[1] = DynModel_B.Sum2_j;
  }

  /* End of Derivatives for SecondOrderIntegrator: '<S48>/Integrator, Second-Order' */

  /* Derivatives for SecondOrderIntegrator: '<S49>/Integrator, Second-Order' */
  if (DynModel_DW.IntegratorSecondOrder_MODE_a == 0) {
    _rtXdot->IntegratorSecondOrder_CSTATE_n[0] =
      DynModel_X.IntegratorSecondOrder_CSTATE_n[1];
    _rtXdot->IntegratorSecondOrder_CSTATE_n[1] = DynModel_B.Sum2_c;
  }

  /* End of Derivatives for SecondOrderIntegrator: '<S49>/Integrator, Second-Order' */

  /* Derivatives for SecondOrderIntegrator: '<S50>/Integrator, Second-Order' */
  if (DynModel_DW.IntegratorSecondOrder_MODE_pu == 0) {
    _rtXdot->IntegratorSecondOrder_CSTATE_d[0] =
      DynModel_X.IntegratorSecondOrder_CSTATE_d[1];
    _rtXdot->IntegratorSecondOrder_CSTATE_d[1] = DynModel_B.Sum2_p;
  }

  /* End of Derivatives for SecondOrderIntegrator: '<S50>/Integrator, Second-Order' */

  /* Derivatives for Integrator: '<S4>/p,q,r ' */
  _rtXdot->pqr_CSTATE[0] = DynModel_B.Product2_m[0];
  _rtXdot->pqr_CSTATE[1] = DynModel_B.Product2_m[1];
  _rtXdot->pqr_CSTATE[2] = DynModel_B.Product2_m[2];
}

/* Model initialize function */
void DynModel_initialize(void)
{
  /* Registration code */

  /* initialize real-time model */
  (void) memset((void *)DynModel_M, 0,
                sizeof(RT_MODEL_DynModel_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&DynModel_M->solverInfo,
                          &DynModel_M->Timing.simTimeStep);
    rtsiSetTPtr(&DynModel_M->solverInfo, &rtmGetTPtr(DynModel_M));
    rtsiSetStepSizePtr(&DynModel_M->solverInfo, &DynModel_M->Timing.stepSize0);
    rtsiSetdXPtr(&DynModel_M->solverInfo, &DynModel_M->ModelData.derivs);
    rtsiSetContStatesPtr(&DynModel_M->solverInfo, (real_T **)
                         &DynModel_M->ModelData.contStates);
    rtsiSetNumContStatesPtr(&DynModel_M->solverInfo,
      &DynModel_M->Sizes.numContStates);
    rtsiSetErrorStatusPtr(&DynModel_M->solverInfo, (&rtmGetErrorStatus
      (DynModel_M)));
    rtsiSetRTModelPtr(&DynModel_M->solverInfo, DynModel_M);
  }

  rtsiSetSimTimeStep(&DynModel_M->solverInfo, MAJOR_TIME_STEP);
  DynModel_M->ModelData.intgData.y = DynModel_M->ModelData.odeY;
  DynModel_M->ModelData.intgData.f[0] = DynModel_M->ModelData.odeF[0];
  DynModel_M->ModelData.intgData.f[1] = DynModel_M->ModelData.odeF[1];
  DynModel_M->ModelData.intgData.f[2] = DynModel_M->ModelData.odeF[2];
  DynModel_M->ModelData.intgData.f[3] = DynModel_M->ModelData.odeF[3];
  DynModel_M->ModelData.contStates = ((X_DynModel_T *) &DynModel_X);
  rtsiSetSolverData(&DynModel_M->solverInfo, (void *)
                    &DynModel_M->ModelData.intgData);
  rtsiSetSolverName(&DynModel_M->solverInfo,"ode4");
  rtmSetTPtr(DynModel_M, &DynModel_M->Timing.tArray[0]);
  DynModel_M->Timing.stepSize0 = 0.004;
  rtmSetFirstInitCond(DynModel_M, 1);

  /* block I/O */
  (void) memset(((void *) &DynModel_B), 0,
                sizeof(B_DynModel_T));

  /* states (continuous) */
  {
    (void) memset((void *)&DynModel_X, 0,
                  sizeof(X_DynModel_T));
  }

  /* states (dwork) */
  (void) memset((void *)&DynModel_DW, 0,
                sizeof(DW_DynModel_T));

  /* external inputs */
  (void) memset((void *)&DynModel_U, 0,
                sizeof(ExtU_DynModel_T));

  /* external outputs */
  (void) memset((void *)&DynModel_Y, 0,
                sizeof(ExtY_DynModel_T));

  /* Start for If: '<S53>/If' */
  DynModel_DW.If_ActiveSubsystem = -1;

  /* Start for IfAction SubSystem: '<S53>/Negative Trace' */
  /* Start for If: '<S64>/Find Maximum Diagonal Value' */
  DynModel_DW.FindMaximumDiagonalValue_ActiveSubsystem = -1;

  /* End of Start for SubSystem: '<S53>/Negative Trace' */

  /* Start for Merge: '<S53>/Merge' */
  DynModel_B.Merge[0] = 1.0;
  DynModel_B.Merge[1] = 0.0;
  DynModel_B.Merge[2] = 0.0;
  DynModel_B.Merge[3] = 0.0;

  /* ConstCode for Outport: '<Root>/Sonar' */
  DynModel_Y.Sonar = 0.0F;

  {
    uint32_T y;
    real_T y1;

    /* InitializeConditions for Integrator: '<S4>/xe,ye,ze' */
    DynModel_X.xeyeze_CSTATE[0] = 0.0;
    DynModel_X.xeyeze_CSTATE[1] = 0.0;
    DynModel_X.xeyeze_CSTATE[2] = 0.0;

    /* InitializeConditions for RandomNumber: '<S63>/Random Number' */
    DynModel_DW.RandSeed = 1144108930UL;
    DynModel_DW.NextOutput = rt_nrand_Upu32_Yd_f_pw(&DynModel_DW.RandSeed) *
      0.01;

    /* InitializeConditions for RandomNumber: '<S61>/Random Number' */
    DynModel_DW.RandSeed_f = 1144108930UL;
    DynModel_DW.NextOutput_a = rt_nrand_Upu32_Yd_f_pw(&DynModel_DW.RandSeed_f) *
      0.0031622776601683794;

    /* InitializeConditions for Integrator: '<S4>/ub,vb,wb' */
    DynModel_X.ubvbwb_CSTATE[0] = 0.0;
    DynModel_X.ubvbwb_CSTATE[1] = 0.0;
    DynModel_X.ubvbwb_CSTATE[2] = 0.0;

    /* InitializeConditions for RandomNumber: '<S59>/Random Number' */
    DynModel_DW.RandSeed_fw = 1144108930UL;
    DynModel_DW.NextOutput_l = rt_nrand_Upu32_Yd_f_pw(&DynModel_DW.RandSeed_fw) *
      0.0031622776601683794;

    /* InitializeConditions for RandomNumber: '<S51>/Random Number' */
    DynModel_DW.RandSeed_fm = 1144108930UL;
    DynModel_DW.NextOutput_n = rt_nrand_Upu32_Yd_f_pw(&DynModel_DW.RandSeed_fm) *
      0.01;

    /* InitializeConditions for RandomNumber: '<S55>/Random Number' */
    y = 1144108930UL;
    y1 = rt_nrand_Upu32_Yd_f_pw(&y) * 15.0;
    DynModel_DW.NextOutput_o[0] = y1;
    DynModel_DW.RandSeed_e[0] = y;
    y = 1144108930UL;
    y1 = rt_nrand_Upu32_Yd_f_pw(&y) * 15.0;
    DynModel_DW.NextOutput_o[1] = y1;
    DynModel_DW.RandSeed_e[1] = y;
    y = 1144108930UL;
    y1 = rt_nrand_Upu32_Yd_f_pw(&y) * 20.0;
    DynModel_DW.NextOutput_o[2] = y1;
    DynModel_DW.RandSeed_e[2] = y;

    /* InitializeConditions for RandomNumber: '<S55>/Random Number1' */
    y = 1144108930UL;
    y1 = rt_nrand_Upu32_Yd_f_pw(&y) * 0.031622776601683791;
    DynModel_DW.NextOutput_h[0] = y1;
    DynModel_DW.RandSeed_i[0] = y;
    y = 1144108930UL;
    y1 = rt_nrand_Upu32_Yd_f_pw(&y) * 0.031622776601683791;
    DynModel_DW.NextOutput_h[1] = y1;
    DynModel_DW.RandSeed_i[1] = y;
    y = 1144108930UL;
    y1 = rt_nrand_Upu32_Yd_f_pw(&y) * 0.031622776601683791;
    DynModel_DW.NextOutput_h[2] = y1;
    DynModel_DW.RandSeed_i[2] = y;

    /* InitializeConditions for Integrator: '<S8>/q0 q1 q2 q3' */
    if (rtmIsFirstInitCond(DynModel_M)) {
      DynModel_X.q0q1q2q3_CSTATE[0] = 0.0;
      DynModel_X.q0q1q2q3_CSTATE[1] = 0.0;
      DynModel_X.q0q1q2q3_CSTATE[2] = 0.0;
      DynModel_X.q0q1q2q3_CSTATE[3] = 0.0;
    }

    DynModel_DW.q0q1q2q3_IWORK.IcNeedsLoading = 1;

    /* InitializeConditions for RandomNumber: '<S58>/Random Number' */
    DynModel_DW.RandSeed_p = 1144108930UL;
    DynModel_DW.NextOutput_am = rt_nrand_Upu32_Yd_f_pw(&DynModel_DW.RandSeed_p) *
      0.0031622776601683794;

    /* InitializeConditions for RandomNumber: '<S52>/White Noise' */
    DynModel_DW.RandSeed_l = 931168259UL;
    DynModel_DW.NextOutput_lh = rt_nrand_Upu32_Yd_f_pw(&DynModel_DW.RandSeed_l);

    /* InitializeConditions for Memory: '<S2>/Memory2' */
    DynModel_DW.Memory2_PreviousInput = 0.0;

    /* InitializeConditions for SecondOrderIntegrator: '<S47>/Integrator, Second-Order' */
    DynModel_X.IntegratorSecondOrder_CSTATE[0] = 0.0;
    DynModel_X.IntegratorSecondOrder_CSTATE[1] = 0.0;
    DynModel_DW.IntegratorSecondOrder_MODE = 0;

    /* InitializeConditions for SecondOrderIntegrator: '<S48>/Integrator, Second-Order' */
    DynModel_X.IntegratorSecondOrder_CSTATE_h[0] = 0.0;
    DynModel_X.IntegratorSecondOrder_CSTATE_h[1] = 0.0;
    DynModel_DW.IntegratorSecondOrder_MODE_p = 0;

    /* InitializeConditions for SecondOrderIntegrator: '<S49>/Integrator, Second-Order' */
    DynModel_X.IntegratorSecondOrder_CSTATE_n[0] = 0.0;
    DynModel_X.IntegratorSecondOrder_CSTATE_n[1] = 0.0;
    DynModel_DW.IntegratorSecondOrder_MODE_a = 0;

    /* InitializeConditions for SecondOrderIntegrator: '<S50>/Integrator, Second-Order' */
    DynModel_X.IntegratorSecondOrder_CSTATE_d[0] = 0.0;
    DynModel_X.IntegratorSecondOrder_CSTATE_d[1] = 0.0;
    DynModel_DW.IntegratorSecondOrder_MODE_pu = 0;

    /* InitializeConditions for Integrator: '<S4>/p,q,r ' */
    DynModel_X.pqr_CSTATE[0] = 0.0;
    DynModel_X.pqr_CSTATE[1] = 0.0;
    DynModel_X.pqr_CSTATE[2] = 0.0;

    /* InitializeConditions for RandomNumber: '<S130>/White Noise' */
    y = 1373044741UL;
    y1 = rt_nrand_Upu32_Yd_f_pw(&y);
    DynModel_DW.NextOutput_k[0] = y1;
    DynModel_DW.RandSeed_ls[0] = y;
    y = 411009029UL;
    y1 = rt_nrand_Upu32_Yd_f_pw(&y);
    DynModel_DW.NextOutput_k[1] = y1;
    DynModel_DW.RandSeed_ls[1] = y;
    y = 1845526542UL;
    y1 = rt_nrand_Upu32_Yd_f_pw(&y);
    DynModel_DW.NextOutput_k[2] = y1;
    DynModel_DW.RandSeed_ls[2] = y;

    /* InitializeConditions for RandomNumber: '<S147>/White Noise' */
    y = 1689616386UL;
    y1 = rt_nrand_Upu32_Yd_f_pw(&y);
    DynModel_DW.NextOutput_p[0] = y1;
    DynModel_DW.RandSeed_j[0] = y;
    y = 1998225409UL;
    y1 = rt_nrand_Upu32_Yd_f_pw(&y);
    DynModel_DW.NextOutput_p[1] = y1;
    DynModel_DW.RandSeed_j[1] = y;
    y = 1181220867UL;
    y1 = rt_nrand_Upu32_Yd_f_pw(&y);
    DynModel_DW.NextOutput_p[2] = y1;
    DynModel_DW.RandSeed_j[2] = y;

    /* InitializeConditions for MATLAB Function: '<S2>/LiPo Battery' */
    DynModel_DW.discharge = 0.0;

    /* set "at time zero" to false */
    if (rtmIsFirstInitCond(DynModel_M)) {
      rtmSetFirstInitCond(DynModel_M, 0);
    }
  }
}

/* Model terminate function */
void DynModel_terminate(void)
{
  /* (no terminate code required) */
}
