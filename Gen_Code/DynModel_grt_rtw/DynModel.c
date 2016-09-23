/*
 * DynModel.c
 *
 * Code generation for model "DynModel".
 *
 * Model version              : 1.372
 * Simulink Coder version : 8.8 (R2015a) 09-Feb-2015
 * C source code generated on : Tue Sep 06 19:45:50 2016
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
  int_T nXc = 13;
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
  real_T rtb_TSamp;
  real_T rtb_Saturation_ph;
  real_T rtb_Product3_nw;
  real_T rtb_jxi;
  real_T rtb_Switch_i;
  real_T rtb_Saturation_f;
  real_T rtb_Saturation_mk;
  real_T rtb_Saturation1_m;
  real_T rtb_Switch_fk;
  real_T rtb_Sum_j;
  real_T rtb_Switch_d;
  real_T rtb_TrigonometricFunction;
  real_T rtb_Sum_e;
  real_T rtb_Saturation;
  real_T rtb_Saturation1;
  real_T rtb_Saturation2;
  real_T rtb_Saturation3;
  real_T rtb_forces[3];
  real_T rtb_kxj;
  real_T rtb_ixj;
  real_T rtb_Product_k4[3];
  real_T rtb_Saturation_j[3];
  real_T rtb_Sum1_p;
  real_T tmp[3];
  int16_T rtb_Compare_0;
  real_T rtb_VectorConcatenate_idx_0;
  real_T rtb_VectorConcatenate_idx_1;
  real_T rtb_Sum_h_idx_0;
  real_T rtb_Sum_h_idx_1;
  real_T rtb_Sum_h_idx_2;
  real_T rtb_Add2_idx_0;
  real_T rtb_Add2_idx_1;
  real_T rtb_Add2_idx_2;
  real_T rtb_Saturation_l3_idx_2;
  real_T rtb_Saturation_l3_idx_1;
  real_T rtb_Saturation_l3_idx_0;
  real_T rtb_Saturation_g_idx_0;
  real_T rtb_Saturation_g_idx_1;
  real_T rtb_Saturation_g_idx_2;
  real_T rtb_ZeroOrderHold2_idx_0;
  real_T rtb_ZeroOrderHold2_idx_1;
  real_T rtb_ZeroOrderHold2_idx_2;
  real32_T rtb_DataTypeConversion14_idx_0;
  real32_T rtb_DataTypeConversion14_idx_1;
  real32_T rtb_DataTypeConversion14_idx_2;
  real32_T rtb_DataTypeConversion15_idx_0;
  real32_T rtb_DataTypeConversion15_idx_1;
  real32_T rtb_DataTypeConversion15_idx_2;
  real32_T rtb_DataTypeConversion12_idx_0;
  real32_T rtb_DataTypeConversion12_idx_1;
  real32_T rtb_DataTypeConversion12_idx_2;
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

  /* Sum: '<S56>/Sum1' incorporates:
   *  UnaryMinus: '<S56>/Ze2height'
   */
  DynModel_B.Sum1 = -DynModel_B.xeyeze[2];

  /* Saturate: '<S58>/Limit  altitude  to troposhere' */
  if (DynModel_B.Sum1 > 11000.0) {
    rtb_Switch_i = 11000.0;
  } else if (DynModel_B.Sum1 < 0.0) {
    rtb_Switch_i = 0.0;
  } else {
    rtb_Switch_i = DynModel_B.Sum1;
  }

  /* Sum: '<S58>/Sum1' incorporates:
   *  Constant: '<S58>/Sea Level  Temperature'
   *  Gain: '<S58>/Lapse Rate'
   *  Saturate: '<S58>/Limit  altitude  to troposhere'
   */
  DynModel_B.Sum1_e = 288.15 - 0.0065 * rtb_Switch_i;
  if (rtmIsMajorTimeStep(DynModel_M)) {
    /* Sum: '<S65>/Sum2' incorporates:
     *  Constant: '<S65>/K2C'
     *  RandomNumber: '<S65>/Random Number'
     *  Sum: '<S65>/Add1'
     */
    rtb_Saturation_ph = (273.15 + DynModel_B.Sum1_e) + DynModel_DW.NextOutput;

    /* Saturate: '<S65>/Saturation' */
    if (rtb_Saturation_ph > 85.0) {
      rtb_Saturation_ph = 85.0;
    } else {
      if (rtb_Saturation_ph < 40.0) {
        rtb_Saturation_ph = 40.0;
      }
    }

    /* End of Saturate: '<S65>/Saturation' */
  }

  /* Gain: '<S58>/1//T0' */
  rtb_Product3_nw = 0.00347041471455839 * DynModel_B.Sum1_e;

  /* Math: '<S58>/(T//T0)^(g//LR) ' */
  if (rtb_Product3_nw < 0.0) {
    rtb_jxi = -pow(-rtb_Product3_nw, 5.2558756014667134);
  } else {
    rtb_jxi = pow(rtb_Product3_nw, 5.2558756014667134);
  }

  /* End of Math: '<S58>/(T//T0)^(g//LR) ' */

  /* Saturate: '<S58>/Limit  altitude  to Stratosphere' incorporates:
   *  Constant: '<S58>/Altitude of Troposphere'
   *  Sum: '<S58>/Sum'
   */
  if (11000.0 - DynModel_B.Sum1 > 0.0) {
    rtb_Switch_i = 0.0;
  } else if (11000.0 - DynModel_B.Sum1 < -9000.0) {
    rtb_Switch_i = -9000.0;
  } else {
    rtb_Switch_i = 11000.0 - DynModel_B.Sum1;
  }

  /* Math: '<S58>/Stratosphere Model' incorporates:
   *  Gain: '<S58>/g//R'
   *  Product: '<S58>/Product1'
   *  Saturate: '<S58>/Limit  altitude  to Stratosphere'
   *
   * About '<S58>/Stratosphere Model':
   *  Operator: exp
   */
  rtb_Switch_i = exp(1.0 / DynModel_B.Sum1_e * (0.034163191409533639 *
    rtb_Switch_i));

  /* Product: '<S58>/Product2' incorporates:
   *  Gain: '<S58>/P0'
   */
  DynModel_B.Product2 = 101325.0 * rtb_jxi * rtb_Switch_i;
  if (rtmIsMajorTimeStep(DynModel_M)) {
    /* Sum: '<S63>/Sum2' incorporates:
     *  Gain: '<S63>/Bar2mBar'
     *  Gain: '<S63>/Pa2Bar'
     *  RandomNumber: '<S63>/Random Number'
     */
    rtb_Saturation_f = 1.0E-5 * DynModel_B.Product2 * 1000.0 +
      DynModel_DW.NextOutput_a;

    /* Saturate: '<S63>/Saturation' */
    if (rtb_Saturation_f > 1200.0) {
      rtb_Saturation_f = 1200.0;
    } else {
      if (rtb_Saturation_f < 10.0) {
        rtb_Saturation_f = 10.0;
      }
    }

    /* End of Saturate: '<S63>/Saturation' */
  }

  /* Product: '<S58>/Product3' incorporates:
   *  Gain: '<S58>/rho0'
   *  Product: '<S58>/Product'
   */
  DynModel_B.Product3 = rtb_jxi / rtb_Product3_nw * 1.225 * rtb_Switch_i;

  /* Integrator: '<S4>/ub,vb,wb' */
  DynModel_B.ubvbwb[0] = DynModel_X.ubvbwb_CSTATE[0];
  DynModel_B.ubvbwb[1] = DynModel_X.ubvbwb_CSTATE[1];
  DynModel_B.ubvbwb[2] = DynModel_X.ubvbwb_CSTATE[2];
  if (rtmIsMajorTimeStep(DynModel_M)) {
    /* Sum: '<S61>/Sum2' incorporates:
     *  DotProduct: '<S61>/Dot Product'
     *  Gain: '<S61>/Bar2mBar'
     *  Gain: '<S61>/Gain'
     *  Gain: '<S61>/Pa2Bar'
     *  Product: '<S61>/Product'
     *  RandomNumber: '<S61>/Random Number'
     */
    rtb_Saturation_mk = ((DynModel_B.ubvbwb[0] * DynModel_B.ubvbwb[0] +
                          DynModel_B.ubvbwb[1] * DynModel_B.ubvbwb[1]) +
                         DynModel_B.ubvbwb[2] * DynModel_B.ubvbwb[2]) *
      DynModel_B.Product3 * 0.5 * 1.0E-5 * 1000.0 + DynModel_DW.NextOutput_l;

    /* Saturate: '<S61>/Saturation' */
    if (rtb_Saturation_mk > 1000.0) {
      rtb_Saturation_mk = 1000.0;
    } else {
      if (rtb_Saturation_mk < 0.0) {
        rtb_Saturation_mk = 0.0;
      }
    }

    /* End of Saturate: '<S61>/Saturation' */

    /* Sum: '<S53>/Add1' incorporates:
     *  RandomNumber: '<S53>/Random Number'
     */
    rtb_Saturation1_m = DynModel_B.Sum1 + DynModel_DW.NextOutput_n;

    /* Saturate: '<S53>/Saturation1' */
    if (!(rtb_Saturation1_m >= 0.0)) {
      rtb_Saturation1_m = 0.0;
    }

    /* End of Saturate: '<S53>/Saturation1' */

    /* Sum: '<S109>/Sum' incorporates:
     *  Gain: '<S113>/Unit Conversion'
     *  Product: '<S112>/rad lat'
     *  RandomNumber: '<S57>/Random Number'
     *  Sum: '<S57>/Add1'
     */
    rtb_jxi = (DynModel_DW.NextOutput_o[0] + DynModel_B.xeyeze[0]) *
      1.5708579706943943E-7 * 57.295779513082323 + 43.718691;

    /* Switch: '<S117>/Switch' incorporates:
     *  Abs: '<S117>/Abs'
     *  Bias: '<S117>/Bias'
     *  Bias: '<S117>/Bias1'
     *  Constant: '<S117>/Constant2'
     *  Math: '<S117>/Math Function1'
     */
    if (fabs(rtb_jxi) > 180.0) {
      rtb_jxi = rt_modd(rtb_jxi + 180.0, 360.0) + -180.0;
    }

    /* End of Switch: '<S117>/Switch' */

    /* Abs: '<S114>/Abs1' */
    rtb_Switch_i = fabs(rtb_jxi);

    /* Switch: '<S114>/Switch' incorporates:
     *  Bias: '<S114>/Bias'
     *  Bias: '<S114>/Bias1'
     *  Constant: '<S110>/Constant'
     *  Constant: '<S110>/Constant1'
     *  Constant: '<S116>/Constant'
     *  Gain: '<S114>/Gain'
     *  Product: '<S114>/Divide1'
     *  RelationalOperator: '<S116>/Compare'
     *  Signum: '<S114>/Sign1'
     *  Switch: '<S110>/Switch1'
     */
    if ((rtb_Switch_i > 90.0) > 0) {
      /* Signum: '<S114>/Sign1' */
      if (rtb_jxi < 0.0) {
        rtb_jxi = -1.0;
      } else {
        if (rtb_jxi > 0.0) {
          rtb_jxi = 1.0;
        }
      }

      rtb_Switch_fk = (-(rtb_Switch_i + -90.0) + 90.0) * rtb_jxi;
      rtb_Compare_0 = 180;
    } else {
      rtb_Switch_fk = rtb_jxi;
      rtb_Compare_0 = 0;
    }

    /* End of Switch: '<S114>/Switch' */

    /* Sum: '<S110>/Sum' incorporates:
     *  Gain: '<S113>/Unit Conversion'
     *  Product: '<S112>/rad long '
     *  RandomNumber: '<S57>/Random Number'
     *  Sum: '<S109>/Sum'
     *  Sum: '<S57>/Add1'
     */
    rtb_Sum_j = ((DynModel_DW.NextOutput_o[1] + DynModel_B.xeyeze[1]) *
                 2.1658460268129011E-7 * 57.295779513082323 +
                 DynModel_ConstB.Switch_d) + (real_T)rtb_Compare_0;

    /* Switch: '<S115>/Switch' incorporates:
     *  Abs: '<S115>/Abs'
     *  Bias: '<S115>/Bias'
     *  Bias: '<S115>/Bias1'
     *  Constant: '<S115>/Constant2'
     *  Math: '<S115>/Math Function1'
     */
    if (fabs(rtb_Sum_j) > 180.0) {
      rtb_Sum_j = rt_modd(rtb_Sum_j + 180.0, 360.0) + -180.0;
    }

    /* End of Switch: '<S115>/Switch' */

    /* Sum: '<S109>/Sum1' incorporates:
     *  RandomNumber: '<S57>/Random Number'
     *  Sum: '<S57>/Add1'
     *  UnaryMinus: '<S109>/Ze2height'
     */
    rtb_Sum1_p = -(DynModel_DW.NextOutput_o[2] + DynModel_B.xeyeze[2]);

    /* RandomNumber: '<S57>/Random Number1' */
    rtb_Add2_idx_0 = DynModel_DW.NextOutput_h[0];
    rtb_Add2_idx_1 = DynModel_DW.NextOutput_h[1];
    rtb_Add2_idx_2 = DynModel_DW.NextOutput_h[2];
  }

  /* Integrator: '<S9>/q0 q1 q2 q3' */
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

  /* Sqrt: '<S31>/sqrt' incorporates:
   *  Product: '<S32>/Product'
   *  Product: '<S32>/Product1'
   *  Product: '<S32>/Product2'
   *  Product: '<S32>/Product3'
   *  Sum: '<S32>/Sum'
   */
  rtb_Switch_i = sqrt(((rtb_q0q1q2q3[0] * rtb_q0q1q2q3[0] + rtb_q0q1q2q3[1] *
                        rtb_q0q1q2q3[1]) + rtb_q0q1q2q3[2] * rtb_q0q1q2q3[2]) +
                      rtb_q0q1q2q3[3] * rtb_q0q1q2q3[3]);

  /* Product: '<S30>/Product' */
  rtb_jxi = rtb_q0q1q2q3[0] / rtb_Switch_i;

  /* Product: '<S30>/Product1' */
  rtb_Product3_nw = rtb_q0q1q2q3[1] / rtb_Switch_i;

  /* Product: '<S30>/Product2' */
  rtb_Switch_d = rtb_q0q1q2q3[2] / rtb_Switch_i;

  /* Product: '<S30>/Product3' */
  rtb_Switch_i = rtb_q0q1q2q3[3] / rtb_Switch_i;

  /* Sum: '<S20>/Sum' incorporates:
   *  Product: '<S20>/Product'
   *  Product: '<S20>/Product1'
   *  Product: '<S20>/Product2'
   *  Product: '<S20>/Product3'
   */
  DynModel_B.VectorConcatenate[0] = ((rtb_jxi * rtb_jxi + rtb_Product3_nw *
    rtb_Product3_nw) - rtb_Switch_d * rtb_Switch_d) - rtb_Switch_i *
    rtb_Switch_i;

  /* Gain: '<S23>/Gain' incorporates:
   *  Product: '<S23>/Product2'
   *  Product: '<S23>/Product3'
   *  Sum: '<S23>/Sum'
   */
  DynModel_B.VectorConcatenate[1] = (rtb_Product3_nw * rtb_Switch_d -
    rtb_Switch_i * rtb_jxi) * 2.0;

  /* Gain: '<S26>/Gain' incorporates:
   *  Product: '<S26>/Product1'
   *  Product: '<S26>/Product2'
   *  Sum: '<S26>/Sum'
   */
  DynModel_B.VectorConcatenate[2] = (rtb_jxi * rtb_Switch_d + rtb_Product3_nw *
    rtb_Switch_i) * 2.0;

  /* Gain: '<S21>/Gain' incorporates:
   *  Product: '<S21>/Product2'
   *  Product: '<S21>/Product3'
   *  Sum: '<S21>/Sum'
   */
  DynModel_B.VectorConcatenate[3] = (rtb_Switch_i * rtb_jxi + rtb_Product3_nw *
    rtb_Switch_d) * 2.0;

  /* Sum: '<S24>/Sum' incorporates:
   *  Product: '<S24>/Product'
   *  Product: '<S24>/Product1'
   *  Product: '<S24>/Product2'
   *  Product: '<S24>/Product3'
   */
  DynModel_B.VectorConcatenate[4] = ((rtb_jxi * rtb_jxi - rtb_Product3_nw *
    rtb_Product3_nw) + rtb_Switch_d * rtb_Switch_d) - rtb_Switch_i *
    rtb_Switch_i;

  /* Gain: '<S27>/Gain' incorporates:
   *  Product: '<S27>/Product1'
   *  Product: '<S27>/Product2'
   *  Sum: '<S27>/Sum'
   */
  DynModel_B.VectorConcatenate[5] = (rtb_Switch_d * rtb_Switch_i - rtb_jxi *
    rtb_Product3_nw) * 2.0;

  /* Gain: '<S22>/Gain' incorporates:
   *  Product: '<S22>/Product1'
   *  Product: '<S22>/Product2'
   *  Sum: '<S22>/Sum'
   */
  DynModel_B.VectorConcatenate[6] = (rtb_Product3_nw * rtb_Switch_i - rtb_jxi *
    rtb_Switch_d) * 2.0;

  /* Gain: '<S25>/Gain' incorporates:
   *  Product: '<S25>/Product1'
   *  Product: '<S25>/Product2'
   *  Sum: '<S25>/Sum'
   */
  DynModel_B.VectorConcatenate[7] = (rtb_jxi * rtb_Product3_nw + rtb_Switch_d *
    rtb_Switch_i) * 2.0;

  /* Sum: '<S28>/Sum' incorporates:
   *  Product: '<S28>/Product'
   *  Product: '<S28>/Product1'
   *  Product: '<S28>/Product2'
   *  Product: '<S28>/Product3'
   */
  DynModel_B.VectorConcatenate[8] = ((rtb_jxi * rtb_jxi - rtb_Product3_nw *
    rtb_Product3_nw) - rtb_Switch_d * rtb_Switch_d) + rtb_Switch_i *
    rtb_Switch_i;

  /* Product: '<S15>/Product' incorporates:
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

  /* End of Product: '<S15>/Product' */
  if (rtmIsMajorTimeStep(DynModel_M)) {
    /* Sum: '<S57>/Add2' */
    rtb_Add2_idx_0 += DynModel_B.Product[0];
    rtb_Add2_idx_1 += DynModel_B.Product[1];
    rtb_Add2_idx_2 += DynModel_B.Product[2];

    /* Trigonometry: '<S57>/Trigonometric Function' */
    rtb_TrigonometricFunction = atan2(rtb_Add2_idx_1, rtb_Add2_idx_0);
  }

  /* Sum: '<S56>/Sum' incorporates:
   *  Gain: '<S94>/Unit Conversion'
   *  Product: '<S93>/rad lat'
   *  Product: '<S93>/x*cos'
   */
  rtb_jxi = DynModel_B.xeyeze[0] * 1.5708579706943943E-7 * 57.295779513082323 +
    43.718691;

  /* Switch: '<S98>/Switch' incorporates:
   *  Abs: '<S98>/Abs'
   *  Bias: '<S98>/Bias'
   *  Bias: '<S98>/Bias1'
   *  Constant: '<S98>/Constant2'
   *  Math: '<S98>/Math Function1'
   */
  if (fabs(rtb_jxi) > 180.0) {
    rtb_jxi = rt_modd(rtb_jxi + 180.0, 360.0) + -180.0;
  }

  /* End of Switch: '<S98>/Switch' */

  /* Abs: '<S95>/Abs1' */
  rtb_Switch_i = fabs(rtb_jxi);

  /* Switch: '<S95>/Switch' incorporates:
   *  Bias: '<S95>/Bias'
   *  Bias: '<S95>/Bias1'
   *  Constant: '<S91>/Constant'
   *  Constant: '<S91>/Constant1'
   *  Constant: '<S97>/Constant'
   *  Gain: '<S95>/Gain'
   *  Product: '<S95>/Divide1'
   *  RelationalOperator: '<S97>/Compare'
   *  Signum: '<S95>/Sign1'
   *  Switch: '<S91>/Switch1'
   */
  if ((rtb_Switch_i > 90.0) > 0) {
    /* Signum: '<S95>/Sign1' */
    if (rtb_jxi < 0.0) {
      rtb_jxi = -1.0;
    } else {
      if (rtb_jxi > 0.0) {
        rtb_jxi = 1.0;
      }
    }

    rtb_Switch_d = (-(rtb_Switch_i + -90.0) + 90.0) * rtb_jxi;
    rtb_Compare_0 = 180;
  } else {
    rtb_Switch_d = rtb_jxi;
    rtb_Compare_0 = 0;
  }

  /* End of Switch: '<S95>/Switch' */

  /* Sum: '<S91>/Sum' incorporates:
   *  Gain: '<S94>/Unit Conversion'
   *  Product: '<S93>/rad long '
   *  Product: '<S93>/y*cos'
   *  Sum: '<S56>/Sum'
   */
  rtb_Sum_e = (2.1658460268129011E-7 * DynModel_B.xeyeze[1] * 57.295779513082323
               + DynModel_ConstB.Switch_b) + (real_T)rtb_Compare_0;
  if (rtmIsMajorTimeStep(DynModel_M)) {
    /* Product: '<S60>/Product2' */
    rtb_Sum_h_idx_0 = 0.86602540378443871;
    rtb_Sum_h_idx_1 = 0.49999999999999994;
    rtb_Sum_h_idx_2 = 0.0;

    /* Sum: '<S60>/Sum2' incorporates:
     *  Product: '<S60>/Matrix Multiply2'
     *  RandomNumber: '<S60>/Random Number'
     *  Saturate: '<S60>/Saturation'
     */
    for (rtb_Compare_0 = 0; rtb_Compare_0 < 3; rtb_Compare_0++) {
      rtb_forces[rtb_Compare_0] = (DynModel_B.VectorConcatenate[rtb_Compare_0 +
        3] * 0.49999999999999994 + DynModel_B.VectorConcatenate[rtb_Compare_0] *
        0.86602540378443871) + DynModel_DW.NextOutput_am;
    }

    /* End of Sum: '<S60>/Sum2' */

    /* Saturate: '<S60>/Saturation' incorporates:
     *  Product: '<S60>/Matrix Multiply2'
     */
    if (rtb_forces[0] > 2.0) {
      rtb_Saturation_g_idx_0 = 2.0;
    } else if (rtb_forces[0] < -2.0) {
      rtb_Saturation_g_idx_0 = -2.0;
    } else {
      rtb_Saturation_g_idx_0 = rtb_forces[0];
    }

    if (rtb_forces[1] > 2.0) {
      rtb_Saturation_g_idx_1 = 2.0;
    } else if (rtb_forces[1] < -2.0) {
      rtb_Saturation_g_idx_1 = -2.0;
    } else {
      rtb_Saturation_g_idx_1 = rtb_forces[1];
    }

    if (rtb_forces[2] > 2.0) {
      rtb_Saturation_g_idx_2 = 2.0;
    } else if (rtb_forces[2] < -2.0) {
      rtb_Saturation_g_idx_2 = -2.0;
    } else {
      rtb_Saturation_g_idx_2 = rtb_forces[2];
    }

    /* Gain: '<S54>/Output' incorporates:
     *  RandomNumber: '<S54>/White Noise'
     */
    DynModel_B.Output = 0.00019364916731037085 * DynModel_DW.NextOutput_lh;

    /* SampleTimeMath: '<S43>/TSamp'
     *
     * About '<S43>/TSamp':
     *  y = u * K where K = 1 / ( w * Ts )
     */
    rtb_TSamp = DynModel_B.Sum1 * 250.0;

    /* Gain: '<S5>/Gain1' incorporates:
     *  Sum: '<S43>/Diff'
     *  UnitDelay: '<S43>/UD'
     */
    DynModel_B.Gain1[0] = 0.0;
    DynModel_B.Gain1[1] = 0.0;
    DynModel_B.Gain1[2] = (rtb_TSamp - DynModel_DW.UD_DSTATE) * 50.0;
  }

  /* Switch: '<S5>/Switch' incorporates:
   *  Constant: '<S5>/Constant2'
   *  Gain: '<S5>/Gain4'
   *  Product: '<S5>/Matrix Multiply'
   *  Sum: '<S5>/Add8'
   */
  if (-DynModel_B.Sum1 >= 0.0) {
    /* Sum: '<S5>/Add8' incorporates:
     *  Constant: '<S5>/Constant2'
     *  Gain: '<S5>/Gain3'
     *  Product: '<S5>/Matrix Multiply'
     */
    rtb_Switch_i = (10000.0 * DynModel_B.Sum1 + -11.772) + DynModel_B.Gain1[2];
    for (rtb_Compare_0 = 0; rtb_Compare_0 < 3; rtb_Compare_0++) {
      DynModel_Y.Forces[rtb_Compare_0] = 0.0;
      DynModel_Y.Forces[rtb_Compare_0] +=
        DynModel_B.VectorConcatenate[rtb_Compare_0] * DynModel_B.Gain1[0];
      DynModel_Y.Forces[rtb_Compare_0] +=
        DynModel_B.VectorConcatenate[rtb_Compare_0 + 3] * DynModel_B.Gain1[1];
      DynModel_Y.Forces[rtb_Compare_0] +=
        DynModel_B.VectorConcatenate[rtb_Compare_0 + 6] * rtb_Switch_i;
    }
  } else {
    DynModel_Y.Forces[0] = 0.0;
    DynModel_Y.Forces[1] = 0.0;
    DynModel_Y.Forces[2] = 0.0;
  }

  /* End of Switch: '<S5>/Switch' */
  if (rtmIsMajorTimeStep(DynModel_M)) {
    /* Memory: '<S2>/Memory2' */
    DynModel_B.Memory2 = DynModel_DW.Memory2_PreviousInput;
  }

  /* Saturate: '<S2>/Saturation' incorporates:
   *  Inport: '<Root>/PWM1'
   */
  if (DynModel_U.PWM1 > 1.0) {
    rtb_Saturation = 1.0;
  } else if (DynModel_U.PWM1 < 0.0) {
    rtb_Saturation = 0.0;
  } else {
    rtb_Saturation = DynModel_U.PWM1;
  }

  /* End of Saturate: '<S2>/Saturation' */

  /* Saturate: '<S2>/Saturation1' incorporates:
   *  Inport: '<Root>/PWM2'
   */
  if (DynModel_U.PWM2 > 1.0) {
    rtb_Saturation1 = 1.0;
  } else if (DynModel_U.PWM2 < 0.0) {
    rtb_Saturation1 = 0.0;
  } else {
    rtb_Saturation1 = DynModel_U.PWM2;
  }

  /* End of Saturate: '<S2>/Saturation1' */

  /* Saturate: '<S2>/Saturation2' incorporates:
   *  Inport: '<Root>/PWM3'
   */
  if (DynModel_U.PWM3 > 1.0) {
    rtb_Saturation2 = 1.0;
  } else if (DynModel_U.PWM3 < 0.0) {
    rtb_Saturation2 = 0.0;
  } else {
    rtb_Saturation2 = DynModel_U.PWM3;
  }

  /* End of Saturate: '<S2>/Saturation2' */

  /* Saturate: '<S2>/Saturation3' incorporates:
   *  Inport: '<Root>/PWM4'
   */
  if (DynModel_U.PWM4 > 1.0) {
    rtb_Saturation3 = 1.0;
  } else if (DynModel_U.PWM4 < 0.0) {
    rtb_Saturation3 = 0.0;
  } else {
    rtb_Saturation3 = DynModel_U.PWM4;
  }

  /* End of Saturate: '<S2>/Saturation3' */

  /* MATLAB Function: '<S7>/multicopter' incorporates:
   *  SignalConversion: '<S46>/TmpSignal ConversionAt SFunction Inport4'
   */
  /* MATLAB Function 'DynModel/Dynamics/Subsystem/multicopter': '<S46>:1' */
  /* ===============================Parameters================================= */
  /*  */
  /*  a = [ax;ay;az];     Vector with the cross-sectional areas */
  /*  M = [M];            Frame Mass */
  /*  l = [l];            Lenght of the Quadcopter arm */
  /*  Kthr = [Kthr];      Coeff. for the computation of thrust */
  /*  Ktrq = [Ktrq];      Coeff. for the computation of the torque */
  /* ================================Constants================================= */
  /* '<S46>:1:16' */
  /*  coefficient of drag                                       */
  /* ==============================Actuator Mixer============================== */
  /*  [x(roll), y(pitch), z(yaw)] */
  /*  QUAD [X] */
  /* '<S46>:1:23' */
  /*  MIX = [  0   1  -1   ;      % QUAD [+] */
  /*          -1   0   1   ; */
  /*           0  -1  -1   ; */
  /*           1   0   1   ]; */
  /* ==================================Forces==================================  */
  /*  We are evaluating in Body frame     */
  /* '<S46>:1:35' */
  /*  gravitational force */
  /*  drag force */
  /* '<S46>:1:36' */
  /* --------------------------------Thrust Model------------------------------ */
  /* thrust = Kt * density * rotors.^2;                  % rotor thrust */
  /* '<S46>:1:40' */
  DynModel_Y.Thrusts[0] = 5.0 * rtb_Saturation;
  DynModel_Y.Thrusts[1] = 5.0 * rtb_Saturation1;
  DynModel_Y.Thrusts[2] = 5.0 * rtb_Saturation2;
  DynModel_Y.Thrusts[3] = 5.0 * rtb_Saturation3;

  /* Sum: '<S2>/Add7' incorporates:
   *  Constant: '<S7>/h_ref12'
   *  MATLAB Function: '<S7>/multicopter'
   */
  /* -------------------------------------------------------------------------- */
  /* '<S46>:1:43' */
  /* '<S46>:1:44' */
  /* ==================================Moments================================= */
  /*  Thrusts contributions to momentum */
  /* '<S46>:1:49' */
  /*  x moment */
  /* '<S46>:1:50' */
  /*  y moment */
  /*  Torques contributions */
  /* momentum_x = sum(abs(MIX(:, 3)) .* rotor_inertia .* rotors) * omega(1);     % x rotor momentum */
  /* momentum_y = sum(abs(MIX(:, 3)) .* rotor_inertia .* rotors) * omega(2);     % y rotor momentum */
  /* momentum_z = sum(MIX(:, 3) .* rotor_inertia .* rotors) * omega(3);          % z rotor momentum */
  /* --------------------------------Torque Model------------------------------ */
  /* mz = MIX(:,3) .* Kq * density .* rotors.^2;              % rotor torque */
  /* '<S46>:1:59' */
  /* -------------------------------------------------------------------------- */
  /* '<S46>:1:62' */
  /*  - [momentum_x; momentum_y; momentum_z]; */
  /* ========================================================================== */
  DynModel_Y.Forces[0] += DynModel_B.VectorConcatenate[6] * 9.81 * 1.2 +
    -DynModel_B.ubvbwb[0] * DynModel_B.Memory2 * 0.016813708498984763;
  DynModel_Y.Forces[1] += DynModel_B.VectorConcatenate[7] * 9.81 * 1.2 +
    -DynModel_B.ubvbwb[1] * DynModel_B.Memory2 * 0.018813708498984762;
  DynModel_Y.Forces[2] += (DynModel_B.VectorConcatenate[8] * 9.81 * 1.2 +
    -DynModel_B.ubvbwb[2] * DynModel_B.Memory2 * 0.18845573684677208) -
    (((DynModel_Y.Thrusts[0] + DynModel_Y.Thrusts[1]) + DynModel_Y.Thrusts[2]) +
     DynModel_Y.Thrusts[3]);

  /* Product: '<S4>/Product' incorporates:
   *  Constant: '<S11>/Constant'
   */
  DynModel_B.Product_b[0] = DynModel_Y.Forces[0] / 1.2;
  DynModel_B.Product_b[1] = DynModel_Y.Forces[1] / 1.2;
  DynModel_B.Product_b[2] = DynModel_Y.Forces[2] / 1.2;
  if (rtmIsMajorTimeStep(DynModel_M)) {
    /* ZeroOrderHold: '<S129>/Zero-Order Hold1' */
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
    /* ZeroOrderHold: '<S129>/Zero-Order Hold2' */
    rtb_ZeroOrderHold2_idx_0 = DynModel_B.MatrixMultiply1[0];
    rtb_ZeroOrderHold2_idx_1 = DynModel_B.MatrixMultiply1[1];
    rtb_ZeroOrderHold2_idx_2 = DynModel_B.MatrixMultiply1[2];
  }

  /* Integrator: '<S4>/p,q,r ' */
  DynModel_B.pqr[0] = DynModel_X.pqr_CSTATE[0];
  DynModel_B.pqr[1] = DynModel_X.pqr_CSTATE[1];
  DynModel_B.pqr[2] = DynModel_X.pqr_CSTATE[2];
  if (rtmIsMajorTimeStep(DynModel_M)) {
    /* ZeroOrderHold: '<S129>/Zero-Order Hold' */
    rtb_Product_k4[0] = DynModel_B.pqr[0];
    rtb_Product_k4[1] = DynModel_B.pqr[1];
    rtb_Product_k4[2] = DynModel_B.pqr[2];
  }

  /* Sqrt: '<S34>/sqrt' incorporates:
   *  Product: '<S35>/Product'
   *  Product: '<S35>/Product1'
   *  Product: '<S35>/Product2'
   *  Product: '<S35>/Product3'
   *  Sum: '<S35>/Sum'
   */
  rtb_Product3_nw = sqrt(((rtb_q0q1q2q3[0] * rtb_q0q1q2q3[0] + rtb_q0q1q2q3[1] *
    rtb_q0q1q2q3[1]) + rtb_q0q1q2q3[2] * rtb_q0q1q2q3[2]) + rtb_q0q1q2q3[3] *
    rtb_q0q1q2q3[3]);

  /* Product: '<S33>/Product' */
  rtb_jxi = rtb_q0q1q2q3[0] / rtb_Product3_nw;

  /* Product: '<S33>/Product1' */
  rtb_kxj = rtb_q0q1q2q3[1] / rtb_Product3_nw;

  /* Product: '<S33>/Product2' */
  rtb_ixj = rtb_q0q1q2q3[2] / rtb_Product3_nw;

  /* Product: '<S33>/Product3' */
  rtb_Product3_nw = rtb_q0q1q2q3[3] / rtb_Product3_nw;

  /* Trigonometry: '<S17>/Trigonometric Function1' incorporates:
   *  Fcn: '<S17>/fcn1'
   *  Fcn: '<S17>/fcn2'
   */
  rtb_VectorConcatenate_idx_0 = atan2((rtb_kxj * rtb_ixj + rtb_jxi *
    rtb_Product3_nw) * 2.0, ((rtb_jxi * rtb_jxi + rtb_kxj * rtb_kxj) - rtb_ixj *
    rtb_ixj) - rtb_Product3_nw * rtb_Product3_nw);

  /* Fcn: '<S17>/fcn3' */
  rtb_Switch_i = (rtb_kxj * rtb_Product3_nw - rtb_jxi * rtb_ixj) * -2.0;

  /* Trigonometry: '<S17>/trigFcn' */
  if (rtb_Switch_i > 1.0) {
    rtb_Switch_i = 1.0;
  } else {
    if (rtb_Switch_i < -1.0) {
      rtb_Switch_i = -1.0;
    }
  }

  rtb_VectorConcatenate_idx_1 = asin(rtb_Switch_i);

  /* End of Trigonometry: '<S17>/trigFcn' */

  /* Fcn: '<S17>/fcn4' */
  rtb_Switch_i = (rtb_ixj * rtb_Product3_nw + rtb_jxi * rtb_kxj) * 2.0;

  /* Fcn: '<S17>/fcn5' */
  rtb_jxi = ((rtb_jxi * rtb_jxi - rtb_kxj * rtb_kxj) - rtb_ixj * rtb_ixj) +
    rtb_Product3_nw * rtb_Product3_nw;

  /* Trigonometry: '<S17>/Trigonometric Function3' */
  rtb_kxj = atan2(rtb_Switch_i, rtb_jxi);
  if (rtmIsMajorTimeStep(DynModel_M)) {
    /* Gain: '<S129>/Gain' incorporates:
     *  Constant: '<S129>/wl_ins'
     *  Constant: '<S3>/center of gravity'
     *  Sum: '<S129>/Sum7'
     */
    rtb_Saturation_l3_idx_0 = 0.0;
    rtb_Saturation_l3_idx_1 = 0.0;
    rtb_Saturation_l3_idx_2 = 0.0;

    /* Sum: '<S140>/Sum' incorporates:
     *  Product: '<S142>/i x j'
     *  Product: '<S142>/j x k'
     *  Product: '<S142>/k x i'
     *  Product: '<S143>/i x k'
     *  Product: '<S143>/j x i'
     *  Product: '<S143>/k x j'
     */
    rtb_Saturation_j[0] = 0.0;
    rtb_Saturation_j[1] = 0.0;
    rtb_Saturation_j[2] = 0.0;
  }

  /* Switch: '<S5>/Switch1' incorporates:
   *  Gain: '<S5>/Gain'
   *  Gain: '<S5>/Gain4'
   */
  if (-DynModel_B.Sum1 >= 0.0) {
    DynModel_Y.Torques[0] = -rtb_kxj;
    DynModel_Y.Torques[1] = -rtb_VectorConcatenate_idx_1;
    DynModel_Y.Torques[2] = -rtb_VectorConcatenate_idx_0;
  } else {
    DynModel_Y.Torques[0] = 0.0;
    DynModel_Y.Torques[1] = 0.0;
    DynModel_Y.Torques[2] = 0.0;
  }

  /* End of Switch: '<S5>/Switch1' */

  /* Sum: '<S2>/Add6' incorporates:
   *  Constant: '<S7>/h_ref10'
   *  Gain: '<S2>/Gain'
   *  MATLAB Function: '<S7>/multicopter'
   *  SignalConversion: '<S46>/TmpSignal ConversionAt SFunction Inport4'
   */
  DynModel_Y.Torques[0] = ((((DynModel_Y.Thrusts[0] * 0.2 * 1.4142135623730951 /
    2.0 + -DynModel_Y.Thrusts[1] * 0.2 * 1.4142135623730951 / 2.0) +
    -DynModel_Y.Thrusts[2] * 0.2 * 1.4142135623730951 / 2.0) +
    DynModel_Y.Thrusts[3] * 0.2 * 1.4142135623730951 / 2.0) +
    DynModel_Y.Torques[0]) - 0.01 * DynModel_B.pqr[0];
  DynModel_Y.Torques[1] = ((((DynModel_Y.Thrusts[0] * 0.2 * 1.4142135623730951 /
    2.0 + DynModel_Y.Thrusts[1] * 0.2 * 1.4142135623730951 / 2.0) +
    -DynModel_Y.Thrusts[2] * 0.2 * 1.4142135623730951 / 2.0) +
    -DynModel_Y.Thrusts[3] * 0.2 * 1.4142135623730951 / 2.0) +
    DynModel_Y.Torques[1]) - 0.01 * DynModel_B.pqr[1];
  DynModel_Y.Torques[2] = ((((-0.05 * rtb_Saturation + 0.05 * rtb_Saturation1) +
    -0.05 * rtb_Saturation2) + 0.05 * rtb_Saturation3) + DynModel_Y.Torques[2])
    - 0.01 * DynModel_B.pqr[2];

  /* Product: '<S36>/Product' */
  for (rtb_Compare_0 = 0; rtb_Compare_0 < 3; rtb_Compare_0++) {
    rtb_forces[rtb_Compare_0] = DynModel_ConstB.Selector[rtb_Compare_0 + 6] *
      DynModel_B.pqr[2] + (DynModel_ConstB.Selector[rtb_Compare_0 + 3] *
      DynModel_B.pqr[1] + DynModel_ConstB.Selector[rtb_Compare_0] *
      DynModel_B.pqr[0]);
  }

  /* End of Product: '<S36>/Product' */

  /* Product: '<S10>/Product2' incorporates:
   *  Product: '<S39>/i x j'
   *  Product: '<S39>/j x k'
   *  Product: '<S39>/k x i'
   *  Product: '<S40>/i x k'
   *  Product: '<S40>/j x i'
   *  Product: '<S40>/k x j'
   *  Sum: '<S10>/Sum2'
   *  Sum: '<S38>/Sum'
   */
  tmp[0] = DynModel_Y.Torques[0] - (DynModel_B.pqr[1] * rtb_forces[2] -
    DynModel_B.pqr[2] * rtb_forces[1]);
  tmp[1] = DynModel_Y.Torques[1] - (DynModel_B.pqr[2] * rtb_forces[0] -
    DynModel_B.pqr[0] * rtb_forces[2]);
  tmp[2] = DynModel_Y.Torques[2] - (DynModel_B.pqr[0] * rtb_forces[1] -
    DynModel_B.pqr[1] * rtb_forces[0]);
  rt_mrdivide_U1d1x3_U2d3x3_Yd1x3(tmp, DynModel_ConstB.Selector2,
    DynModel_B.Product2_m);
  if (rtmIsMajorTimeStep(DynModel_M)) {
    /* Sum: '<S129>/Sum' */
    rtb_Sum_h_idx_0 = (rtb_Sum_h_idx_0 - rtb_ZeroOrderHold2_idx_0) +
      rtb_Saturation_j[0];
    rtb_Sum_h_idx_1 = (rtb_Sum_h_idx_1 - rtb_ZeroOrderHold2_idx_1) +
      rtb_Saturation_j[1];
    rtb_Switch_i = (rtb_Sum_h_idx_2 - rtb_ZeroOrderHold2_idx_2) +
      rtb_Saturation_j[2];

    /* Product: '<S129>/Product' incorporates:
     *  Constant: '<S129>/Scale factors & Cross-coupling  errors'
     */
    for (rtb_Compare_0 = 0; rtb_Compare_0 < 3; rtb_Compare_0++) {
      rtb_Saturation_j[rtb_Compare_0] = DynModel_ConstP.pooled23[rtb_Compare_0 +
        6] * rtb_Switch_i + (DynModel_ConstP.pooled23[rtb_Compare_0 + 3] *
        rtb_Sum_h_idx_1 + DynModel_ConstP.pooled23[rtb_Compare_0] *
        rtb_Sum_h_idx_0);
    }

    /* End of Product: '<S129>/Product' */

    /* Sum: '<S129>/Sum4' */
    rtb_forces[1] = rtb_Saturation_j[1];
    rtb_forces[2] = rtb_Saturation_j[2];

    /* Saturate: '<S129>/Saturation' incorporates:
     *  Gain: '<S132>/Output'
     *  RandomNumber: '<S132>/White Noise'
     *  Sum: '<S129>/Sum1'
     *  Sum: '<S129>/Sum4'
     */
    rtb_Switch_i = 0.011180339887498949 * DynModel_DW.NextOutput_k[0] +
      rtb_Saturation_j[0];
    if (rtb_Switch_i > 19.62) {
      rtb_Saturation_j[0] = 19.62;
    } else if (rtb_Switch_i < -19.62) {
      rtb_Saturation_j[0] = -19.62;
    } else {
      rtb_Saturation_j[0] = rtb_Switch_i;
    }

    rtb_Switch_i = 0.011180339887498949 * DynModel_DW.NextOutput_k[1] +
      rtb_forces[1];
    if (rtb_Switch_i > 19.62) {
      rtb_Saturation_j[1] = 19.62;
    } else if (rtb_Switch_i < -19.62) {
      rtb_Saturation_j[1] = -19.62;
    } else {
      rtb_Saturation_j[1] = rtb_Switch_i;
    }

    rtb_Switch_i = 0.011180339887498949 * DynModel_DW.NextOutput_k[2] +
      rtb_forces[2];
    if (rtb_Switch_i > 19.62) {
      rtb_Saturation_j[2] = 19.62;
    } else if (rtb_Switch_i < -19.62) {
      rtb_Saturation_j[2] = -19.62;
    } else {
      rtb_Saturation_j[2] = rtb_Switch_i;
    }

    /* End of Saturate: '<S129>/Saturation' */

    /* ZeroOrderHold: '<S130>/Zero-Order Hold' */
    rtb_Saturation_l3_idx_0 = DynModel_B.pqr[0];
    rtb_Saturation_l3_idx_1 = DynModel_B.pqr[1];
    rtb_Saturation_l3_idx_2 = DynModel_B.pqr[2];

    /* Product: '<S130>/Product' incorporates:
     *  Constant: '<S130>/Scale factors & Cross-coupling  errors '
     *  ZeroOrderHold: '<S130>/Zero-Order Hold'
     */
    for (rtb_Compare_0 = 0; rtb_Compare_0 < 3; rtb_Compare_0++) {
      rtb_Product_k4[rtb_Compare_0] = DynModel_ConstP.pooled23[rtb_Compare_0 + 6]
        * DynModel_B.pqr[2] + (DynModel_ConstP.pooled23[rtb_Compare_0 + 3] *
        DynModel_B.pqr[1] + DynModel_ConstP.pooled23[rtb_Compare_0] *
        DynModel_B.pqr[0]);
    }

    /* End of Product: '<S130>/Product' */
  }

  /* Gain: '<S128>/Unit Conversion' */
  DynModel_B.UnitConversion[0] = 0.10197162129779283 * DynModel_B.Product_b[0];
  DynModel_B.UnitConversion[1] = 0.10197162129779283 * DynModel_B.Product_b[1];
  DynModel_B.UnitConversion[2] = 0.10197162129779283 * DynModel_B.Product_b[2];
  if (rtmIsMajorTimeStep(DynModel_M)) {
    /* Saturate: '<S130>/Saturation' incorporates:
     *  Gain: '<S149>/Output'
     *  RandomNumber: '<S149>/White Noise'
     *  Sum: '<S130>/Sum1'
     *  Sum: '<S130>/Sum4'
     */
    rtb_Saturation_l3_idx_0 = 0.00070710678118654762 * DynModel_DW.NextOutput_p
      [0] + rtb_Product_k4[0];
    if (rtb_Saturation_l3_idx_0 > 4.36) {
      rtb_Saturation_l3_idx_0 = 4.36;
    } else {
      if (rtb_Saturation_l3_idx_0 < -4.36) {
        rtb_Saturation_l3_idx_0 = -4.36;
      }
    }

    rtb_Saturation_l3_idx_1 = 0.00070710678118654762 * DynModel_DW.NextOutput_p
      [1] + rtb_Product_k4[1];
    if (rtb_Saturation_l3_idx_1 > 4.36) {
      rtb_Saturation_l3_idx_1 = 4.36;
    } else {
      if (rtb_Saturation_l3_idx_1 < -4.36) {
        rtb_Saturation_l3_idx_1 = -4.36;
      }
    }

    rtb_Saturation_l3_idx_2 = 0.00070710678118654762 * DynModel_DW.NextOutput_p
      [2] + rtb_Product_k4[2];
    if (rtb_Saturation_l3_idx_2 > 4.36) {
      rtb_Saturation_l3_idx_2 = 4.36;
    } else {
      if (rtb_Saturation_l3_idx_2 < -4.36) {
        rtb_Saturation_l3_idx_2 = -4.36;
      }
    }

    /* End of Saturate: '<S130>/Saturation' */
  }

  /* Sum: '<S68>/Add' */
  rtb_Switch_i = (DynModel_B.VectorConcatenate[0] +
                  DynModel_B.VectorConcatenate[4]) +
    DynModel_B.VectorConcatenate[8];

  /* If: '<S55>/If' incorporates:
   *  Sum: '<S68>/Add'
   */
  if (rtmIsMajorTimeStep(DynModel_M)) {
    DynModel_DW.If_ActiveSubsystem = (int8_T)!(rtb_Switch_i > 0.0);
  }

  switch (DynModel_DW.If_ActiveSubsystem) {
   case 0L:
    /* Outputs for IfAction SubSystem: '<S55>/Positive Trace' incorporates:
     *  ActionPort: '<S67>/Action Port'
     */
    /* Sqrt: '<S67>/sqrt' incorporates:
     *  Constant: '<S67>/Constant'
     *  Sum: '<S67>/Sum'
     *  Sum: '<S68>/Add'
     */
    rtb_Switch_i = sqrt(rtb_Switch_i + 1.0);

    /* Gain: '<S67>/Gain' */
    DynModel_B.Merge[0] = 0.5 * rtb_Switch_i;

    /* Gain: '<S67>/Gain1' */
    rtb_Switch_i *= 2.0;

    /* Product: '<S67>/Product' incorporates:
     *  Sum: '<S88>/Add'
     *  Sum: '<S89>/Add'
     *  Sum: '<S90>/Add'
     */
    DynModel_B.Merge[1] = (DynModel_B.VectorConcatenate[7] -
      DynModel_B.VectorConcatenate[5]) / rtb_Switch_i;
    DynModel_B.Merge[2] = (DynModel_B.VectorConcatenate[2] -
      DynModel_B.VectorConcatenate[6]) / rtb_Switch_i;
    DynModel_B.Merge[3] = (DynModel_B.VectorConcatenate[3] -
      DynModel_B.VectorConcatenate[1]) / rtb_Switch_i;

    /* End of Outputs for SubSystem: '<S55>/Positive Trace' */
    break;

   case 1L:
    /* Outputs for IfAction SubSystem: '<S55>/Negative Trace' incorporates:
     *  ActionPort: '<S66>/Action Port'
     */
    /* If: '<S66>/Find Maximum Diagonal Value' */
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
      /* Outputs for IfAction SubSystem: '<S66>/Maximum Value at DCM(2,2)' incorporates:
       *  ActionPort: '<S70>/Action Port'
       */
      /* Sqrt: '<S70>/sqrt' incorporates:
       *  Constant: '<S82>/Constant'
       *  Sum: '<S82>/Add'
       */
      rtb_Switch_i = sqrt(((DynModel_B.VectorConcatenate[4] -
                            DynModel_B.VectorConcatenate[0]) -
                           DynModel_B.VectorConcatenate[8]) + 1.0);

      /* Gain: '<S70>/Gain' */
      DynModel_B.Merge[2] = 0.5 * rtb_Switch_i;

      /* Switch: '<S81>/Switch' incorporates:
       *  Constant: '<S81>/Constant1'
       */
      if (rtb_Switch_i != 0.0) {
        rtb_jxi = 0.5;
      } else {
        rtb_jxi = 0.0;
        rtb_Switch_i = 1.0;
      }

      /* End of Switch: '<S81>/Switch' */

      /* Product: '<S81>/Product' */
      rtb_Switch_i = rtb_jxi / rtb_Switch_i;

      /* Gain: '<S70>/Gain1' incorporates:
       *  Product: '<S70>/Product'
       *  Sum: '<S80>/Add'
       */
      DynModel_B.Merge[1] = (DynModel_B.VectorConcatenate[1] +
        DynModel_B.VectorConcatenate[3]) * rtb_Switch_i;

      /* Gain: '<S70>/Gain3' incorporates:
       *  Product: '<S70>/Product'
       *  Sum: '<S79>/Add'
       */
      DynModel_B.Merge[3] = (DynModel_B.VectorConcatenate[5] +
        DynModel_B.VectorConcatenate[7]) * rtb_Switch_i;

      /* Gain: '<S70>/Gain4' incorporates:
       *  Product: '<S70>/Product'
       *  Sum: '<S78>/Add'
       */
      DynModel_B.Merge[0] = (DynModel_B.VectorConcatenate[2] -
        DynModel_B.VectorConcatenate[6]) * rtb_Switch_i;

      /* End of Outputs for SubSystem: '<S66>/Maximum Value at DCM(2,2)' */
      break;

     case 1L:
      /* Outputs for IfAction SubSystem: '<S66>/Maximum Value at DCM(3,3)' incorporates:
       *  ActionPort: '<S71>/Action Port'
       */
      /* Sqrt: '<S71>/sqrt' incorporates:
       *  Constant: '<S87>/Constant'
       *  Sum: '<S87>/Add'
       */
      rtb_Switch_i = sqrt(((DynModel_B.VectorConcatenate[8] -
                            DynModel_B.VectorConcatenate[0]) -
                           DynModel_B.VectorConcatenate[4]) + 1.0);

      /* Gain: '<S71>/Gain' */
      DynModel_B.Merge[3] = 0.5 * rtb_Switch_i;

      /* Switch: '<S86>/Switch' incorporates:
       *  Constant: '<S86>/Constant1'
       */
      if (rtb_Switch_i != 0.0) {
        rtb_jxi = 0.5;
      } else {
        rtb_jxi = 0.0;
        rtb_Switch_i = 1.0;
      }

      /* End of Switch: '<S86>/Switch' */

      /* Product: '<S86>/Product' */
      rtb_Switch_i = rtb_jxi / rtb_Switch_i;

      /* Gain: '<S71>/Gain1' incorporates:
       *  Product: '<S71>/Product'
       *  Sum: '<S83>/Add'
       */
      DynModel_B.Merge[1] = (DynModel_B.VectorConcatenate[2] +
        DynModel_B.VectorConcatenate[6]) * rtb_Switch_i;

      /* Gain: '<S71>/Gain2' incorporates:
       *  Product: '<S71>/Product'
       *  Sum: '<S84>/Add'
       */
      DynModel_B.Merge[2] = (DynModel_B.VectorConcatenate[5] +
        DynModel_B.VectorConcatenate[7]) * rtb_Switch_i;

      /* Gain: '<S71>/Gain3' incorporates:
       *  Product: '<S71>/Product'
       *  Sum: '<S85>/Add'
       */
      DynModel_B.Merge[0] = (DynModel_B.VectorConcatenate[3] -
        DynModel_B.VectorConcatenate[1]) * rtb_Switch_i;

      /* End of Outputs for SubSystem: '<S66>/Maximum Value at DCM(3,3)' */
      break;

     case 2L:
      /* Outputs for IfAction SubSystem: '<S66>/Maximum Value at DCM(1,1)' incorporates:
       *  ActionPort: '<S69>/Action Port'
       */
      /* Sqrt: '<S69>/sqrt' incorporates:
       *  Constant: '<S77>/Constant'
       *  Sum: '<S77>/Add'
       */
      rtb_Switch_i = sqrt(((DynModel_B.VectorConcatenate[0] -
                            DynModel_B.VectorConcatenate[4]) -
                           DynModel_B.VectorConcatenate[8]) + 1.0);

      /* Gain: '<S69>/Gain' */
      DynModel_B.Merge[1] = 0.5 * rtb_Switch_i;

      /* Switch: '<S76>/Switch' incorporates:
       *  Constant: '<S76>/Constant1'
       */
      if (rtb_Switch_i != 0.0) {
        rtb_jxi = 0.5;
      } else {
        rtb_jxi = 0.0;
        rtb_Switch_i = 1.0;
      }

      /* End of Switch: '<S76>/Switch' */

      /* Product: '<S76>/Product' */
      rtb_Switch_i = rtb_jxi / rtb_Switch_i;

      /* Gain: '<S69>/Gain1' incorporates:
       *  Product: '<S69>/Product'
       *  Sum: '<S75>/Add'
       */
      DynModel_B.Merge[2] = (DynModel_B.VectorConcatenate[1] +
        DynModel_B.VectorConcatenate[3]) * rtb_Switch_i;

      /* Gain: '<S69>/Gain2' incorporates:
       *  Product: '<S69>/Product'
       *  Sum: '<S73>/Add'
       */
      DynModel_B.Merge[3] = (DynModel_B.VectorConcatenate[2] +
        DynModel_B.VectorConcatenate[6]) * rtb_Switch_i;

      /* Gain: '<S69>/Gain3' incorporates:
       *  Product: '<S69>/Product'
       *  Sum: '<S74>/Add'
       */
      DynModel_B.Merge[0] = (DynModel_B.VectorConcatenate[7] -
        DynModel_B.VectorConcatenate[5]) * rtb_Switch_i;

      /* End of Outputs for SubSystem: '<S66>/Maximum Value at DCM(1,1)' */
      break;
    }

    /* End of If: '<S66>/Find Maximum Diagonal Value' */
    /* End of Outputs for SubSystem: '<S55>/Negative Trace' */
    break;
  }

  /* End of If: '<S55>/If' */

  /* DotProduct: '<S19>/Dot Product' */
  rtb_Switch_i = ((rtb_q0q1q2q3[0] * rtb_q0q1q2q3[0] + rtb_q0q1q2q3[1] *
                   rtb_q0q1q2q3[1]) + rtb_q0q1q2q3[2] * rtb_q0q1q2q3[2]) +
    rtb_q0q1q2q3[3] * rtb_q0q1q2q3[3];

  /* Fcn: '<S19>/q0dot' incorporates:
   *  Constant: '<S19>/Constant'
   *  DotProduct: '<S19>/Dot Product'
   *  Sum: '<S19>/Sum'
   */
  DynModel_B.q0dot = ((rtb_q0q1q2q3[1] * DynModel_B.pqr[0] + rtb_q0q1q2q3[2] *
                       DynModel_B.pqr[1]) + rtb_q0q1q2q3[3] * DynModel_B.pqr[2])
    * -0.5 + (1.0 - rtb_Switch_i) * rtb_q0q1q2q3[0];

  /* Fcn: '<S19>/q1dot' incorporates:
   *  Constant: '<S19>/Constant'
   *  DotProduct: '<S19>/Dot Product'
   *  Sum: '<S19>/Sum'
   */
  DynModel_B.q1dot = ((rtb_q0q1q2q3[0] * DynModel_B.pqr[0] + rtb_q0q1q2q3[2] *
                       DynModel_B.pqr[2]) - rtb_q0q1q2q3[3] * DynModel_B.pqr[1])
    * 0.5 + (1.0 - rtb_Switch_i) * rtb_q0q1q2q3[1];

  /* Fcn: '<S19>/q2dot' incorporates:
   *  Constant: '<S19>/Constant'
   *  DotProduct: '<S19>/Dot Product'
   *  Sum: '<S19>/Sum'
   */
  DynModel_B.q2dot = ((rtb_q0q1q2q3[0] * DynModel_B.pqr[1] + rtb_q0q1q2q3[3] *
                       DynModel_B.pqr[0]) - rtb_q0q1q2q3[1] * DynModel_B.pqr[2])
    * 0.5 + (1.0 - rtb_Switch_i) * rtb_q0q1q2q3[2];

  /* Fcn: '<S19>/q3dot' incorporates:
   *  Constant: '<S19>/Constant'
   *  DotProduct: '<S19>/Dot Product'
   *  Sum: '<S19>/Sum'
   */
  DynModel_B.q3dot = ((rtb_q0q1q2q3[0] * DynModel_B.pqr[2] + rtb_q0q1q2q3[1] *
                       DynModel_B.pqr[1]) - rtb_q0q1q2q3[2] * DynModel_B.pqr[0])
    * 0.5 + (1.0 - rtb_Switch_i) * rtb_q0q1q2q3[3];

  /* Sum: '<S4>/Sum' incorporates:
   *  Product: '<S41>/i x j'
   *  Product: '<S41>/j x k'
   *  Product: '<S41>/k x i'
   *  Product: '<S42>/i x k'
   *  Product: '<S42>/j x i'
   *  Product: '<S42>/k x j'
   *  Sum: '<S12>/Sum'
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
  if (rtmIsMajorTimeStep(DynModel_M)) {
    /* DataTypeConversion: '<S1>/Data Type Conversion14' */
    rtb_DataTypeConversion14_idx_0 = (real32_T)rtb_Saturation_j[0];
    rtb_DataTypeConversion14_idx_1 = (real32_T)rtb_Saturation_j[1];
    rtb_DataTypeConversion14_idx_2 = (real32_T)rtb_Saturation_j[2];

    /* DataTypeConversion: '<S1>/Data Type Conversion15' */
    rtb_DataTypeConversion15_idx_0 = (real32_T)rtb_Saturation_l3_idx_0;
    rtb_DataTypeConversion15_idx_1 = (real32_T)rtb_Saturation_l3_idx_1;
    rtb_DataTypeConversion15_idx_2 = (real32_T)rtb_Saturation_l3_idx_2;

    /* DataTypeConversion: '<S1>/Data Type Conversion12' */
    rtb_DataTypeConversion12_idx_0 = (real32_T)rtb_Saturation_g_idx_0;
    rtb_DataTypeConversion12_idx_1 = (real32_T)rtb_Saturation_g_idx_1;
    rtb_DataTypeConversion12_idx_2 = (real32_T)rtb_Saturation_g_idx_2;

    /* Outport: '<Root>/Temp' incorporates:
     *  DataTypeConversion: '<S1>/Data Type Conversion'
     */
    DynModel_Y.Temp = (real32_T)rtb_Saturation_ph;

    /* Outport: '<Root>/Press' incorporates:
     *  DataTypeConversion: '<S1>/Data Type Conversion1'
     */
    DynModel_Y.Press = (real32_T)rtb_Saturation_f;

    /* Outport: '<Root>/diff_Pres' incorporates:
     *  DataTypeConversion: '<S1>/Data Type Conversion3'
     */
    DynModel_Y.diff_Pres = (real32_T)rtb_Saturation_mk;

    /* Outport: '<Root>/Baro_Alt' incorporates:
     *  DataTypeConversion: '<S1>/Data Type Conversion4'
     */
    DynModel_Y.Baro_Alt = (real32_T)rtb_Saturation1_m;

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

  /* Switch: '<S96>/Switch' incorporates:
   *  Abs: '<S96>/Abs'
   */
  if (fabs(rtb_Sum_e) > 180.0) {
    /* Outport: '<Root>/Lat_Lon_Alt' incorporates:
     *  Bias: '<S96>/Bias'
     *  Bias: '<S96>/Bias1'
     *  Constant: '<S96>/Constant2'
     *  DataTypeConversion: '<S1>/Data Type Conversion10'
     *  Math: '<S96>/Math Function1'
     */
    DynModel_Y.Lat_Lon_Alt[1] = (real32_T)(rt_modd(rtb_Sum_e + 180.0, 360.0) +
      -180.0);
  } else {
    /* Outport: '<Root>/Lat_Lon_Alt' incorporates:
     *  DataTypeConversion: '<S1>/Data Type Conversion10'
     */
    DynModel_Y.Lat_Lon_Alt[1] = (real32_T)rtb_Sum_e;
  }

  /* End of Switch: '<S96>/Switch' */

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
   */
  DynModel_Y.RPY[0] = (real32_T)(rtb_kxj + DynModel_B.Output);
  DynModel_Y.RPY[1] = (real32_T)(rtb_VectorConcatenate_idx_1 + DynModel_B.Output);
  DynModel_Y.RPY[2] = (real32_T)(rtb_VectorConcatenate_idx_0 + DynModel_B.Output);
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

  /* Outport: '<Root>/Rotor_Speed' */
  DynModel_Y.Rotor_Speed[0] = rtb_Saturation;
  DynModel_Y.Rotor_Speed[1] = rtb_Saturation1;
  DynModel_Y.Rotor_Speed[2] = rtb_Saturation2;
  DynModel_Y.Rotor_Speed[3] = rtb_Saturation3;
  if (rtmIsMajorTimeStep(DynModel_M)) {
    if (rtmIsMajorTimeStep(DynModel_M)) {
      /* Update for RandomNumber: '<S65>/Random Number' */
      DynModel_DW.NextOutput = rt_nrand_Upu32_Yd_f_pw(&DynModel_DW.RandSeed) *
        0.01;

      /* Update for RandomNumber: '<S63>/Random Number' */
      DynModel_DW.NextOutput_a = rt_nrand_Upu32_Yd_f_pw(&DynModel_DW.RandSeed_f)
        * 0.0031622776601683794;

      /* Update for RandomNumber: '<S61>/Random Number' */
      DynModel_DW.NextOutput_l = rt_nrand_Upu32_Yd_f_pw(&DynModel_DW.RandSeed_fw)
        * 0.0031622776601683794;

      /* Update for RandomNumber: '<S53>/Random Number' */
      DynModel_DW.NextOutput_n = rt_nrand_Upu32_Yd_f_pw(&DynModel_DW.RandSeed_fm)
        * 0.01;

      /* Update for RandomNumber: '<S57>/Random Number' */
      DynModel_DW.NextOutput_o[0] = rt_nrand_Upu32_Yd_f_pw
        (&DynModel_DW.RandSeed_e[0]) * 0.1;
      DynModel_DW.NextOutput_o[1] = rt_nrand_Upu32_Yd_f_pw
        (&DynModel_DW.RandSeed_e[1]) * 0.1;
      DynModel_DW.NextOutput_o[2] = rt_nrand_Upu32_Yd_f_pw
        (&DynModel_DW.RandSeed_e[2]) * 0.1;

      /* Update for RandomNumber: '<S57>/Random Number1' */
      DynModel_DW.NextOutput_h[0] = rt_nrand_Upu32_Yd_f_pw
        (&DynModel_DW.RandSeed_i[0]) * 0.00031622776601683794;
      DynModel_DW.NextOutput_h[1] = rt_nrand_Upu32_Yd_f_pw
        (&DynModel_DW.RandSeed_i[1]) * 0.00031622776601683794;
      DynModel_DW.NextOutput_h[2] = rt_nrand_Upu32_Yd_f_pw
        (&DynModel_DW.RandSeed_i[2]) * 0.00031622776601683794;
    }

    /* Update for Integrator: '<S9>/q0 q1 q2 q3' */
    DynModel_DW.q0q1q2q3_IWORK.IcNeedsLoading = 0;
    if (rtmIsMajorTimeStep(DynModel_M)) {
      /* Update for RandomNumber: '<S60>/Random Number' */
      DynModel_DW.NextOutput_am = rt_nrand_Upu32_Yd_f_pw(&DynModel_DW.RandSeed_p)
        * 0.0031622776601683794;

      /* Update for RandomNumber: '<S54>/White Noise' */
      DynModel_DW.NextOutput_lh = rt_nrand_Upu32_Yd_f_pw(&DynModel_DW.RandSeed_l);

      /* Update for UnitDelay: '<S43>/UD' */
      DynModel_DW.UD_DSTATE = rtb_TSamp;

      /* Update for Memory: '<S2>/Memory2' */
      DynModel_DW.Memory2_PreviousInput = DynModel_B.Product3;

      /* Update for RandomNumber: '<S132>/White Noise' */
      DynModel_DW.NextOutput_k[0] = rt_nrand_Upu32_Yd_f_pw
        (&DynModel_DW.RandSeed_ls[0]);
      DynModel_DW.NextOutput_k[1] = rt_nrand_Upu32_Yd_f_pw
        (&DynModel_DW.RandSeed_ls[1]);
      DynModel_DW.NextOutput_k[2] = rt_nrand_Upu32_Yd_f_pw
        (&DynModel_DW.RandSeed_ls[2]);

      /* Update for RandomNumber: '<S149>/White Noise' */
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

  /* Derivatives for Integrator: '<S9>/q0 q1 q2 q3' */
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

  /* Start for If: '<S55>/If' */
  DynModel_DW.If_ActiveSubsystem = -1;

  /* Start for IfAction SubSystem: '<S55>/Negative Trace' */
  /* Start for If: '<S66>/Find Maximum Diagonal Value' */
  DynModel_DW.FindMaximumDiagonalValue_ActiveSubsystem = -1;

  /* End of Start for SubSystem: '<S55>/Negative Trace' */

  /* Start for Merge: '<S55>/Merge' */
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

    /* InitializeConditions for RandomNumber: '<S65>/Random Number' */
    DynModel_DW.RandSeed = 1144108930UL;
    DynModel_DW.NextOutput = rt_nrand_Upu32_Yd_f_pw(&DynModel_DW.RandSeed) *
      0.01;

    /* InitializeConditions for RandomNumber: '<S63>/Random Number' */
    DynModel_DW.RandSeed_f = 1144108930UL;
    DynModel_DW.NextOutput_a = rt_nrand_Upu32_Yd_f_pw(&DynModel_DW.RandSeed_f) *
      0.0031622776601683794;

    /* InitializeConditions for Integrator: '<S4>/ub,vb,wb' */
    DynModel_X.ubvbwb_CSTATE[0] = 0.0;
    DynModel_X.ubvbwb_CSTATE[1] = 0.0;
    DynModel_X.ubvbwb_CSTATE[2] = 0.0;

    /* InitializeConditions for RandomNumber: '<S61>/Random Number' */
    DynModel_DW.RandSeed_fw = 1144108930UL;
    DynModel_DW.NextOutput_l = rt_nrand_Upu32_Yd_f_pw(&DynModel_DW.RandSeed_fw) *
      0.0031622776601683794;

    /* InitializeConditions for RandomNumber: '<S53>/Random Number' */
    DynModel_DW.RandSeed_fm = 1144108930UL;
    DynModel_DW.NextOutput_n = rt_nrand_Upu32_Yd_f_pw(&DynModel_DW.RandSeed_fm) *
      0.01;

    /* InitializeConditions for RandomNumber: '<S57>/Random Number' */
    y = 1144108930UL;
    y1 = rt_nrand_Upu32_Yd_f_pw(&y) * 0.1;
    DynModel_DW.NextOutput_o[0] = y1;
    DynModel_DW.RandSeed_e[0] = y;
    y = 1144108930UL;
    y1 = rt_nrand_Upu32_Yd_f_pw(&y) * 0.1;
    DynModel_DW.NextOutput_o[1] = y1;
    DynModel_DW.RandSeed_e[1] = y;
    y = 1144108930UL;
    y1 = rt_nrand_Upu32_Yd_f_pw(&y) * 0.1;
    DynModel_DW.NextOutput_o[2] = y1;
    DynModel_DW.RandSeed_e[2] = y;

    /* InitializeConditions for RandomNumber: '<S57>/Random Number1' */
    y = 1144108930UL;
    y1 = rt_nrand_Upu32_Yd_f_pw(&y) * 0.00031622776601683794;
    DynModel_DW.NextOutput_h[0] = y1;
    DynModel_DW.RandSeed_i[0] = y;
    y = 1144108930UL;
    y1 = rt_nrand_Upu32_Yd_f_pw(&y) * 0.00031622776601683794;
    DynModel_DW.NextOutput_h[1] = y1;
    DynModel_DW.RandSeed_i[1] = y;
    y = 1144108930UL;
    y1 = rt_nrand_Upu32_Yd_f_pw(&y) * 0.00031622776601683794;
    DynModel_DW.NextOutput_h[2] = y1;
    DynModel_DW.RandSeed_i[2] = y;

    /* InitializeConditions for Integrator: '<S9>/q0 q1 q2 q3' */
    if (rtmIsFirstInitCond(DynModel_M)) {
      DynModel_X.q0q1q2q3_CSTATE[0] = 0.0;
      DynModel_X.q0q1q2q3_CSTATE[1] = 0.0;
      DynModel_X.q0q1q2q3_CSTATE[2] = 0.0;
      DynModel_X.q0q1q2q3_CSTATE[3] = 0.0;
    }

    DynModel_DW.q0q1q2q3_IWORK.IcNeedsLoading = 1;

    /* InitializeConditions for RandomNumber: '<S60>/Random Number' */
    DynModel_DW.RandSeed_p = 1144108930UL;
    DynModel_DW.NextOutput_am = rt_nrand_Upu32_Yd_f_pw(&DynModel_DW.RandSeed_p) *
      0.0031622776601683794;

    /* InitializeConditions for RandomNumber: '<S54>/White Noise' */
    DynModel_DW.RandSeed_l = 931168259UL;
    DynModel_DW.NextOutput_lh = rt_nrand_Upu32_Yd_f_pw(&DynModel_DW.RandSeed_l);

    /* InitializeConditions for UnitDelay: '<S43>/UD' */
    DynModel_DW.UD_DSTATE = 0.0;

    /* InitializeConditions for Memory: '<S2>/Memory2' */
    DynModel_DW.Memory2_PreviousInput = 0.0;

    /* InitializeConditions for Integrator: '<S4>/p,q,r ' */
    DynModel_X.pqr_CSTATE[0] = 0.0;
    DynModel_X.pqr_CSTATE[1] = 0.0;
    DynModel_X.pqr_CSTATE[2] = 0.0;

    /* InitializeConditions for RandomNumber: '<S132>/White Noise' */
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

    /* InitializeConditions for RandomNumber: '<S149>/White Noise' */
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
