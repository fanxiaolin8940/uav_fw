/*
 * DynModel.c
 *
 * Code generation for model "DynModel".
 *
 * Model version              : 1.411
 * Simulink Coder version : 8.8 (R2015a) 09-Feb-2015
 * C source code generated on : Wed Oct 12 08:49:20 2016
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

/* Exported block parameters */
real_T Kvreact = -50.0;                  /* Variable: Kvreact
                                        * Referenced by: '<S5>/Gain2'
                                        */
real_T Kwreact = -10.0;                  /* Variable: Kwreact
                                        * Referenced by: '<S5>/Gain5'
                                        */

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
  real_T rtb_Saturation;
  real_T rtb_Product3_nw;
  real_T rtb_jxi;
  real_T rtb_Switch_i;
  real_T rtb_Saturation_i;
  real_T rtb_Saturation_gu;
  real_T rtb_Saturation1;
  real_T rtb_Switch_fk;
  real_T rtb_Sum_j;
  real_T rtb_Switch_d;
  real_T rtb_TrigonometricFunction;
  real_T rtb_Sum_e;
  real_T rtb_forces[3];
  real_T rtb_kxj;
  real_T rtb_ixj;
  real_T rtb_Product_k4[3];
  real_T rtb_Saturation_kg[3];
  real_T rtb_Sum1_p;
  real_T rtb_Add6_0[3];
  real_T tmp[3];
  real_T rtb_thrust;
  real_T rtb_Add7;
  int16_T rtb_Compare_0;
  real_T rtb_VectorConcatenate_idx_1;
  real_T rtb_Sum_h_idx_0;
  real_T rtb_Sum_h_idx_1;
  real_T rtb_Sum_h_idx_2;
  real32_T rtb_DataTypeConversion21_idx_0;
  real32_T rtb_DataTypeConversion21_idx_1;
  real32_T rtb_DataTypeConversion21_idx_2;
  real32_T rtb_DataTypeConversion21_idx_3;
  real_T rtb_Add2_h_idx_0;
  real_T rtb_Add2_h_idx_1;
  real_T rtb_Add2_h_idx_2;
  real_T rtb_thrust_idx_0;
  real_T rtb_Saturation_l_idx_2;
  real_T rtb_Saturation_l_idx_1;
  real_T rtb_Saturation_l_idx_0;
  real_T rtb_Add6_idx_1;
  real_T rtb_Add7_idx_0;
  real_T rtb_Add7_idx_1;
  real_T rtb_Saturation_g_idx_0;
  real_T rtb_Saturation_g_idx_1;
  real_T rtb_Saturation_g_idx_2;
  real_T rtb_thrust_idx_1;
  real_T rtb_thrust_idx_2;
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
  real_T u0;
  real_T u1;
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

  /* Sum: '<S61>/Sum1' incorporates:
   *  UnaryMinus: '<S61>/Ze2height'
   */
  DynModel_B.Sum1 = -DynModel_B.xeyeze[2];

  /* Saturate: '<S63>/Limit  altitude  to troposhere' */
  if (DynModel_B.Sum1 > 11000.0) {
    rtb_Switch_i = 11000.0;
  } else if (DynModel_B.Sum1 < 0.0) {
    rtb_Switch_i = 0.0;
  } else {
    rtb_Switch_i = DynModel_B.Sum1;
  }

  /* Sum: '<S63>/Sum1' incorporates:
   *  Constant: '<S63>/Sea Level  Temperature'
   *  Gain: '<S63>/Lapse Rate'
   *  Saturate: '<S63>/Limit  altitude  to troposhere'
   */
  DynModel_B.Sum1_e = 288.15 - 0.0065 * rtb_Switch_i;
  if (rtmIsMajorTimeStep(DynModel_M)) {
    /* Sum: '<S70>/Sum2' incorporates:
     *  Constant: '<S70>/K2C'
     *  RandomNumber: '<S70>/Random Number'
     *  Sum: '<S70>/Add1'
     */
    rtb_Saturation = (273.15 + DynModel_B.Sum1_e) + DynModel_DW.NextOutput;

    /* Saturate: '<S70>/Saturation' */
    if (rtb_Saturation > 85.0) {
      rtb_Saturation = 85.0;
    } else {
      if (rtb_Saturation < 40.0) {
        rtb_Saturation = 40.0;
      }
    }

    /* End of Saturate: '<S70>/Saturation' */
  }

  /* Gain: '<S63>/1//T0' */
  rtb_Product3_nw = 0.00347041471455839 * DynModel_B.Sum1_e;

  /* Math: '<S63>/(T//T0)^(g//LR) ' */
  if (rtb_Product3_nw < 0.0) {
    rtb_jxi = -pow(-rtb_Product3_nw, 5.2558756014667134);
  } else {
    rtb_jxi = pow(rtb_Product3_nw, 5.2558756014667134);
  }

  /* End of Math: '<S63>/(T//T0)^(g//LR) ' */

  /* Saturate: '<S63>/Limit  altitude  to Stratosphere' incorporates:
   *  Constant: '<S63>/Altitude of Troposphere'
   *  Sum: '<S63>/Sum'
   */
  if (11000.0 - DynModel_B.Sum1 > 0.0) {
    rtb_Switch_i = 0.0;
  } else if (11000.0 - DynModel_B.Sum1 < -9000.0) {
    rtb_Switch_i = -9000.0;
  } else {
    rtb_Switch_i = 11000.0 - DynModel_B.Sum1;
  }

  /* Math: '<S63>/Stratosphere Model' incorporates:
   *  Gain: '<S63>/g//R'
   *  Product: '<S63>/Product1'
   *  Saturate: '<S63>/Limit  altitude  to Stratosphere'
   *
   * About '<S63>/Stratosphere Model':
   *  Operator: exp
   */
  rtb_Switch_i = exp(1.0 / DynModel_B.Sum1_e * (0.034163191409533639 *
    rtb_Switch_i));

  /* Product: '<S63>/Product2' incorporates:
   *  Gain: '<S63>/P0'
   */
  DynModel_B.Product2 = 101325.0 * rtb_jxi * rtb_Switch_i;
  if (rtmIsMajorTimeStep(DynModel_M)) {
    /* Sum: '<S68>/Sum2' incorporates:
     *  Gain: '<S68>/Bar2mBar'
     *  Gain: '<S68>/Pa2Bar'
     *  RandomNumber: '<S68>/Random Number'
     */
    rtb_Saturation_i = 1.0E-5 * DynModel_B.Product2 * 1000.0 +
      DynModel_DW.NextOutput_a;

    /* Saturate: '<S68>/Saturation' */
    if (rtb_Saturation_i > 1200.0) {
      rtb_Saturation_i = 1200.0;
    } else {
      if (rtb_Saturation_i < 10.0) {
        rtb_Saturation_i = 10.0;
      }
    }

    /* End of Saturate: '<S68>/Saturation' */
  }

  /* Product: '<S63>/Product3' incorporates:
   *  Gain: '<S63>/rho0'
   *  Product: '<S63>/Product'
   */
  DynModel_B.Product3 = rtb_jxi / rtb_Product3_nw * 1.225 * rtb_Switch_i;

  /* Integrator: '<S4>/ub,vb,wb' */
  DynModel_B.ubvbwb[0] = DynModel_X.ubvbwb_CSTATE[0];
  DynModel_B.ubvbwb[1] = DynModel_X.ubvbwb_CSTATE[1];
  DynModel_B.ubvbwb[2] = DynModel_X.ubvbwb_CSTATE[2];
  if (rtmIsMajorTimeStep(DynModel_M)) {
    /* Sum: '<S66>/Sum2' incorporates:
     *  DotProduct: '<S66>/Dot Product'
     *  Gain: '<S66>/Bar2mBar'
     *  Gain: '<S66>/Gain'
     *  Gain: '<S66>/Pa2Bar'
     *  Product: '<S66>/Product'
     *  RandomNumber: '<S66>/Random Number'
     */
    rtb_Saturation_gu = ((DynModel_B.ubvbwb[0] * DynModel_B.ubvbwb[0] +
                          DynModel_B.ubvbwb[1] * DynModel_B.ubvbwb[1]) +
                         DynModel_B.ubvbwb[2] * DynModel_B.ubvbwb[2]) *
      DynModel_B.Product3 * 0.5 * 1.0E-5 * 1000.0 + DynModel_DW.NextOutput_l;

    /* Saturate: '<S66>/Saturation' */
    if (rtb_Saturation_gu > 1000.0) {
      rtb_Saturation_gu = 1000.0;
    } else {
      if (rtb_Saturation_gu < 0.0) {
        rtb_Saturation_gu = 0.0;
      }
    }

    /* End of Saturate: '<S66>/Saturation' */

    /* Sum: '<S58>/Add1' incorporates:
     *  RandomNumber: '<S58>/Random Number'
     */
    rtb_Saturation1 = DynModel_B.Sum1 + DynModel_DW.NextOutput_n;

    /* Saturate: '<S58>/Saturation1' */
    if (!(rtb_Saturation1 >= 0.0)) {
      rtb_Saturation1 = 0.0;
    }

    /* End of Saturate: '<S58>/Saturation1' */

    /* Sum: '<S114>/Sum' incorporates:
     *  Gain: '<S118>/Unit Conversion'
     *  Product: '<S117>/rad lat'
     *  RandomNumber: '<S62>/Random Number'
     *  Sum: '<S62>/Add1'
     */
    rtb_jxi = (DynModel_DW.NextOutput_o[0] + DynModel_B.xeyeze[0]) *
      1.5708579706943943E-7 * 57.295779513082323 + 43.718691;

    /* Switch: '<S122>/Switch' incorporates:
     *  Abs: '<S122>/Abs'
     *  Bias: '<S122>/Bias'
     *  Bias: '<S122>/Bias1'
     *  Constant: '<S122>/Constant2'
     *  Math: '<S122>/Math Function1'
     */
    if (fabs(rtb_jxi) > 180.0) {
      rtb_jxi = rt_modd(rtb_jxi + 180.0, 360.0) + -180.0;
    }

    /* End of Switch: '<S122>/Switch' */

    /* Abs: '<S119>/Abs1' */
    rtb_Switch_i = fabs(rtb_jxi);

    /* Switch: '<S119>/Switch' incorporates:
     *  Bias: '<S119>/Bias'
     *  Bias: '<S119>/Bias1'
     *  Constant: '<S115>/Constant'
     *  Constant: '<S115>/Constant1'
     *  Constant: '<S121>/Constant'
     *  Gain: '<S119>/Gain'
     *  Product: '<S119>/Divide1'
     *  RelationalOperator: '<S121>/Compare'
     *  Signum: '<S119>/Sign1'
     *  Switch: '<S115>/Switch1'
     */
    if ((rtb_Switch_i > 90.0) > 0) {
      /* Signum: '<S119>/Sign1' */
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

    /* End of Switch: '<S119>/Switch' */

    /* Sum: '<S115>/Sum' incorporates:
     *  Gain: '<S118>/Unit Conversion'
     *  Product: '<S117>/rad long '
     *  RandomNumber: '<S62>/Random Number'
     *  Sum: '<S114>/Sum'
     *  Sum: '<S62>/Add1'
     */
    rtb_Sum_j = ((DynModel_DW.NextOutput_o[1] + DynModel_B.xeyeze[1]) *
                 2.1658460268129011E-7 * 57.295779513082323 +
                 DynModel_ConstB.Switch_d) + (real_T)rtb_Compare_0;

    /* Switch: '<S120>/Switch' incorporates:
     *  Abs: '<S120>/Abs'
     *  Bias: '<S120>/Bias'
     *  Bias: '<S120>/Bias1'
     *  Constant: '<S120>/Constant2'
     *  Math: '<S120>/Math Function1'
     */
    if (fabs(rtb_Sum_j) > 180.0) {
      rtb_Sum_j = rt_modd(rtb_Sum_j + 180.0, 360.0) + -180.0;
    }

    /* End of Switch: '<S120>/Switch' */

    /* Sum: '<S114>/Sum1' incorporates:
     *  RandomNumber: '<S62>/Random Number'
     *  Sum: '<S62>/Add1'
     *  UnaryMinus: '<S114>/Ze2height'
     */
    rtb_Sum1_p = -(DynModel_DW.NextOutput_o[2] + DynModel_B.xeyeze[2]);

    /* RandomNumber: '<S62>/Random Number1' */
    rtb_Add2_h_idx_0 = DynModel_DW.NextOutput_h[0];
    rtb_Add2_h_idx_1 = DynModel_DW.NextOutput_h[1];
    rtb_Add2_h_idx_2 = DynModel_DW.NextOutput_h[2];
  }

  /* Integrator: '<S10>/q0 q1 q2 q3' */
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

  /* Sqrt: '<S32>/sqrt' incorporates:
   *  Product: '<S33>/Product'
   *  Product: '<S33>/Product1'
   *  Product: '<S33>/Product2'
   *  Product: '<S33>/Product3'
   *  Sum: '<S33>/Sum'
   */
  rtb_Switch_i = sqrt(((rtb_q0q1q2q3[0] * rtb_q0q1q2q3[0] + rtb_q0q1q2q3[1] *
                        rtb_q0q1q2q3[1]) + rtb_q0q1q2q3[2] * rtb_q0q1q2q3[2]) +
                      rtb_q0q1q2q3[3] * rtb_q0q1q2q3[3]);

  /* Product: '<S31>/Product' */
  rtb_jxi = rtb_q0q1q2q3[0] / rtb_Switch_i;

  /* Product: '<S31>/Product1' */
  rtb_Product3_nw = rtb_q0q1q2q3[1] / rtb_Switch_i;

  /* Product: '<S31>/Product2' */
  rtb_Switch_d = rtb_q0q1q2q3[2] / rtb_Switch_i;

  /* Product: '<S31>/Product3' */
  rtb_Switch_i = rtb_q0q1q2q3[3] / rtb_Switch_i;

  /* Sum: '<S21>/Sum' incorporates:
   *  Product: '<S21>/Product'
   *  Product: '<S21>/Product1'
   *  Product: '<S21>/Product2'
   *  Product: '<S21>/Product3'
   */
  DynModel_B.VectorConcatenate[0] = ((rtb_jxi * rtb_jxi + rtb_Product3_nw *
    rtb_Product3_nw) - rtb_Switch_d * rtb_Switch_d) - rtb_Switch_i *
    rtb_Switch_i;

  /* Gain: '<S24>/Gain' incorporates:
   *  Product: '<S24>/Product2'
   *  Product: '<S24>/Product3'
   *  Sum: '<S24>/Sum'
   */
  DynModel_B.VectorConcatenate[1] = (rtb_Product3_nw * rtb_Switch_d -
    rtb_Switch_i * rtb_jxi) * 2.0;

  /* Gain: '<S27>/Gain' incorporates:
   *  Product: '<S27>/Product1'
   *  Product: '<S27>/Product2'
   *  Sum: '<S27>/Sum'
   */
  DynModel_B.VectorConcatenate[2] = (rtb_jxi * rtb_Switch_d + rtb_Product3_nw *
    rtb_Switch_i) * 2.0;

  /* Gain: '<S22>/Gain' incorporates:
   *  Product: '<S22>/Product2'
   *  Product: '<S22>/Product3'
   *  Sum: '<S22>/Sum'
   */
  DynModel_B.VectorConcatenate[3] = (rtb_Switch_i * rtb_jxi + rtb_Product3_nw *
    rtb_Switch_d) * 2.0;

  /* Sum: '<S25>/Sum' incorporates:
   *  Product: '<S25>/Product'
   *  Product: '<S25>/Product1'
   *  Product: '<S25>/Product2'
   *  Product: '<S25>/Product3'
   */
  DynModel_B.VectorConcatenate[4] = ((rtb_jxi * rtb_jxi - rtb_Product3_nw *
    rtb_Product3_nw) + rtb_Switch_d * rtb_Switch_d) - rtb_Switch_i *
    rtb_Switch_i;

  /* Gain: '<S28>/Gain' incorporates:
   *  Product: '<S28>/Product1'
   *  Product: '<S28>/Product2'
   *  Sum: '<S28>/Sum'
   */
  DynModel_B.VectorConcatenate[5] = (rtb_Switch_d * rtb_Switch_i - rtb_jxi *
    rtb_Product3_nw) * 2.0;

  /* Gain: '<S23>/Gain' incorporates:
   *  Product: '<S23>/Product1'
   *  Product: '<S23>/Product2'
   *  Sum: '<S23>/Sum'
   */
  DynModel_B.VectorConcatenate[6] = (rtb_Product3_nw * rtb_Switch_i - rtb_jxi *
    rtb_Switch_d) * 2.0;

  /* Gain: '<S26>/Gain' incorporates:
   *  Product: '<S26>/Product1'
   *  Product: '<S26>/Product2'
   *  Sum: '<S26>/Sum'
   */
  DynModel_B.VectorConcatenate[7] = (rtb_jxi * rtb_Product3_nw + rtb_Switch_d *
    rtb_Switch_i) * 2.0;

  /* Sum: '<S29>/Sum' incorporates:
   *  Product: '<S29>/Product'
   *  Product: '<S29>/Product1'
   *  Product: '<S29>/Product2'
   *  Product: '<S29>/Product3'
   */
  DynModel_B.VectorConcatenate[8] = ((rtb_jxi * rtb_jxi - rtb_Product3_nw *
    rtb_Product3_nw) - rtb_Switch_d * rtb_Switch_d) + rtb_Switch_i *
    rtb_Switch_i;

  /* Product: '<S16>/Product' incorporates:
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

  /* End of Product: '<S16>/Product' */
  if (rtmIsMajorTimeStep(DynModel_M)) {
    /* Sum: '<S62>/Add2' */
    rtb_Add2_h_idx_0 += DynModel_B.Product[0];
    rtb_Add2_h_idx_1 += DynModel_B.Product[1];
    rtb_Add2_h_idx_2 += DynModel_B.Product[2];

    /* Trigonometry: '<S62>/Trigonometric Function' */
    rtb_TrigonometricFunction = atan2(rtb_Add2_h_idx_1, rtb_Add2_h_idx_0);
  }

  /* Sum: '<S61>/Sum' incorporates:
   *  Gain: '<S99>/Unit Conversion'
   *  Product: '<S98>/rad lat'
   *  Product: '<S98>/x*cos'
   */
  rtb_jxi = DynModel_B.xeyeze[0] * 1.5708579706943943E-7 * 57.295779513082323 +
    43.718691;

  /* Switch: '<S103>/Switch' incorporates:
   *  Abs: '<S103>/Abs'
   *  Bias: '<S103>/Bias'
   *  Bias: '<S103>/Bias1'
   *  Constant: '<S103>/Constant2'
   *  Math: '<S103>/Math Function1'
   */
  if (fabs(rtb_jxi) > 180.0) {
    rtb_jxi = rt_modd(rtb_jxi + 180.0, 360.0) + -180.0;
  }

  /* End of Switch: '<S103>/Switch' */

  /* Abs: '<S100>/Abs1' */
  rtb_Switch_i = fabs(rtb_jxi);

  /* Switch: '<S100>/Switch' incorporates:
   *  Bias: '<S100>/Bias'
   *  Bias: '<S100>/Bias1'
   *  Constant: '<S102>/Constant'
   *  Constant: '<S96>/Constant'
   *  Constant: '<S96>/Constant1'
   *  Gain: '<S100>/Gain'
   *  Product: '<S100>/Divide1'
   *  RelationalOperator: '<S102>/Compare'
   *  Signum: '<S100>/Sign1'
   *  Switch: '<S96>/Switch1'
   */
  if ((rtb_Switch_i > 90.0) > 0) {
    /* Signum: '<S100>/Sign1' */
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

  /* End of Switch: '<S100>/Switch' */

  /* Sum: '<S96>/Sum' incorporates:
   *  Gain: '<S99>/Unit Conversion'
   *  Product: '<S98>/rad long '
   *  Product: '<S98>/y*cos'
   *  Sum: '<S61>/Sum'
   */
  rtb_Sum_e = (2.1658460268129011E-7 * DynModel_B.xeyeze[1] * 57.295779513082323
               + DynModel_ConstB.Switch_b) + (real_T)rtb_Compare_0;
  if (rtmIsMajorTimeStep(DynModel_M)) {
    /* Product: '<S65>/Product2' */
    rtb_Sum_h_idx_0 = 1.0;
    rtb_Sum_h_idx_1 = 0.0;
    rtb_Sum_h_idx_2 = 0.0;

    /* Sum: '<S65>/Sum2' incorporates:
     *  Product: '<S65>/Matrix Multiply2'
     *  RandomNumber: '<S65>/Random Number'
     *  Saturate: '<S65>/Saturation'
     */
    for (rtb_Compare_0 = 0; rtb_Compare_0 < 3; rtb_Compare_0++) {
      rtb_forces[rtb_Compare_0] = DynModel_B.VectorConcatenate[rtb_Compare_0] +
        DynModel_DW.NextOutput_am;
    }

    /* End of Sum: '<S65>/Sum2' */

    /* Saturate: '<S65>/Saturation' incorporates:
     *  Product: '<S65>/Matrix Multiply2'
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

    /* Gain: '<S59>/Output' incorporates:
     *  RandomNumber: '<S59>/White Noise'
     */
    DynModel_B.Output = 0.00019364916731037085 * DynModel_DW.NextOutput_lh;

    /* SampleTimeMath: '<S44>/TSamp'
     *
     * About '<S44>/TSamp':
     *  y = u * K where K = 1 / ( w * Ts )
     */
    rtb_TSamp = DynModel_B.Sum1 * 250.0;

    /* Gain: '<S5>/Gain1' incorporates:
     *  Sum: '<S44>/Diff'
     *  UnitDelay: '<S44>/UD'
     */
    DynModel_B.Gain1[0] = 0.0;
    DynModel_B.Gain1[1] = 0.0;
    DynModel_B.Gain1[2] = (rtb_TSamp - DynModel_DW.UD_DSTATE) * 50.0;
  }

  /* Switch: '<S5>/Switch' incorporates:
   *  Gain: '<S5>/Gain4'
   *  Sum: '<S5>/Add1'
   */
  if (-DynModel_B.Sum1 >= 0.0) {
    /* Sum: '<S5>/Add8' incorporates:
     *  Constant: '<S5>/Constant2'
     *  Gain: '<S5>/Gain3'
     *  Product: '<S5>/Matrix Multiply'
     */
    rtb_Switch_i = (10000.0 * DynModel_B.Sum1 + -11.772) + DynModel_B.Gain1[2];

    /* Product: '<S5>/Matrix Multiply1' incorporates:
     *  Gain: '<S5>/Gain2'
     *  Sum: '<S5>/Add1'
     */
    for (rtb_Compare_0 = 0; rtb_Compare_0 < 3; rtb_Compare_0++) {
      rtb_forces[rtb_Compare_0] = DynModel_B.VectorConcatenate[rtb_Compare_0 + 6]
        * (Kvreact * DynModel_B.Product[2]) +
        (DynModel_B.VectorConcatenate[rtb_Compare_0 + 3] * (Kvreact *
          DynModel_B.Product[1]) + Kvreact * DynModel_B.Product[0] *
         DynModel_B.VectorConcatenate[rtb_Compare_0]);
    }

    /* End of Product: '<S5>/Matrix Multiply1' */

    /* Product: '<S5>/Matrix Multiply' incorporates:
     *  Constant: '<S5>/Constant2'
     *  Sum: '<S5>/Add1'
     *  Sum: '<S5>/Add8'
     */
    for (rtb_Compare_0 = 0; rtb_Compare_0 < 3; rtb_Compare_0++) {
      tmp[rtb_Compare_0] = DynModel_B.VectorConcatenate[rtb_Compare_0 + 6] *
        rtb_Switch_i + (DynModel_B.VectorConcatenate[rtb_Compare_0 + 3] *
                        DynModel_B.Gain1[1] +
                        DynModel_B.VectorConcatenate[rtb_Compare_0] *
                        DynModel_B.Gain1[0]);
    }

    rtb_Add7_idx_0 = rtb_forces[0] + tmp[0];
    rtb_Add7_idx_1 = rtb_forces[1] + tmp[1];
    rtb_jxi = rtb_forces[2] + tmp[2];
  } else {
    rtb_Add7_idx_0 = 0.0;
    rtb_Add7_idx_1 = 0.0;
    rtb_jxi = 0.0;
  }

  /* End of Switch: '<S5>/Switch' */
  if (rtmIsMajorTimeStep(DynModel_M)) {
    /* Memory: '<S2>/Memory2' */
    DynModel_B.Memory2 = DynModel_DW.Memory2_PreviousInput;

    /* Saturate: '<S2>/Saturation' incorporates:
     *  Inport: '<Root>/PWM1'
     */
    if (DynModel_U.PWM1 > 1.0) {
      rtb_Switch_i = 1.0;
    } else if (DynModel_U.PWM1 < 0.0) {
      rtb_Switch_i = 0.0;
    } else {
      rtb_Switch_i = DynModel_U.PWM1;
    }

    /* DiscreteTransferFcn: '<S45>/Discrete Transfer Fcn' incorporates:
     *  Gain: '<S45>/Gain1'
     *  Gain: '<S7>/PWM2V'
     *  Gain: '<S7>/RPM2RADS'
     *  Gain: '<S7>/V2RPM'
     *  Saturate: '<S2>/Saturation'
     */
    DynModel_DW.DiscreteTransferFcn_tmp = 12.0 * rtb_Switch_i * 950.0 *
      0.10471975511965977 * 0.029126213592233007 - -0.9417475728155339 *
      DynModel_DW.DiscreteTransferFcn_states;
    DynModel_B.DiscreteTransferFcn = DynModel_DW.DiscreteTransferFcn_tmp +
      DynModel_DW.DiscreteTransferFcn_states;

    /* Saturate: '<S2>/Saturation1' incorporates:
     *  Inport: '<Root>/PWM2'
     */
    if (DynModel_U.PWM2 > 1.0) {
      rtb_Switch_i = 1.0;
    } else if (DynModel_U.PWM2 < 0.0) {
      rtb_Switch_i = 0.0;
    } else {
      rtb_Switch_i = DynModel_U.PWM2;
    }

    /* DiscreteTransferFcn: '<S46>/Discrete Transfer Fcn' incorporates:
     *  Gain: '<S46>/Gain1'
     *  Gain: '<S7>/PWM2V'
     *  Gain: '<S7>/RPM2RADS'
     *  Gain: '<S7>/V2RPM'
     *  Saturate: '<S2>/Saturation1'
     */
    DynModel_DW.DiscreteTransferFcn_tmp_i = 12.0 * rtb_Switch_i * 950.0 *
      0.10471975511965977 * 0.029126213592233007 - -0.9417475728155339 *
      DynModel_DW.DiscreteTransferFcn_states_m;
    DynModel_B.DiscreteTransferFcn_o = DynModel_DW.DiscreteTransferFcn_tmp_i +
      DynModel_DW.DiscreteTransferFcn_states_m;

    /* Saturate: '<S2>/Saturation2' incorporates:
     *  Inport: '<Root>/PWM3'
     */
    if (DynModel_U.PWM3 > 1.0) {
      rtb_Switch_i = 1.0;
    } else if (DynModel_U.PWM3 < 0.0) {
      rtb_Switch_i = 0.0;
    } else {
      rtb_Switch_i = DynModel_U.PWM3;
    }

    /* DiscreteTransferFcn: '<S47>/Discrete Transfer Fcn' incorporates:
     *  Gain: '<S47>/Gain1'
     *  Gain: '<S7>/PWM2V'
     *  Gain: '<S7>/RPM2RADS'
     *  Gain: '<S7>/V2RPM'
     *  Saturate: '<S2>/Saturation2'
     */
    DynModel_DW.DiscreteTransferFcn_tmp_b = 12.0 * rtb_Switch_i * 950.0 *
      0.10471975511965977 * 0.029126213592233007 - -0.9417475728155339 *
      DynModel_DW.DiscreteTransferFcn_states_c;
    DynModel_B.DiscreteTransferFcn_i = DynModel_DW.DiscreteTransferFcn_tmp_b +
      DynModel_DW.DiscreteTransferFcn_states_c;

    /* Saturate: '<S2>/Saturation3' incorporates:
     *  Inport: '<Root>/PWM4'
     */
    if (DynModel_U.PWM4 > 1.0) {
      rtb_Switch_i = 1.0;
    } else if (DynModel_U.PWM4 < 0.0) {
      rtb_Switch_i = 0.0;
    } else {
      rtb_Switch_i = DynModel_U.PWM4;
    }

    /* DiscreteTransferFcn: '<S48>/Discrete Transfer Fcn' incorporates:
     *  Gain: '<S48>/Gain1'
     *  Gain: '<S7>/PWM2V'
     *  Gain: '<S7>/RPM2RADS'
     *  Gain: '<S7>/V2RPM'
     *  Saturate: '<S2>/Saturation3'
     */
    DynModel_DW.DiscreteTransferFcn_tmp_e = 12.0 * rtb_Switch_i * 950.0 *
      0.10471975511965977 * 0.029126213592233007 - -0.9417475728155339 *
      DynModel_DW.DiscreteTransferFcn_states_a;
    DynModel_B.DiscreteTransferFcn_ib = DynModel_DW.DiscreteTransferFcn_tmp_e +
      DynModel_DW.DiscreteTransferFcn_states_a;
  }

  /* MATLAB Function: '<S8>/multicopter' incorporates:
   *  Constant: '<S8>/h_ref8'
   *  SignalConversion: '<S51>/TmpSignal ConversionAt SFunction Inport4'
   */
  /* MATLAB Function 'DynModel/Dynamics/Subsystem/multicopter': '<S51>:1' */
  /* ===============================Parameters================================= */
  /*  */
  /*  a = [ax;ay;az];     Vector with the cross-sectional areas */
  /*  M = [M];            Frame Mass */
  /*  l = [l];            Lenght of the Quadcopter arm */
  /*  Kthr = [Kthr];      Coeff. for the computation of thrust */
  /*  Ktrq = [Ktrq];      Coeff. for the computation of the torque */
  /*  Cd = [Cd];          Coeff. of Drag */
  /* ==============================Actuator Mixer============================== */
  /*  [x(roll), y(pitch), z(yaw)] */
  /*  QUAD [X] */
  /* '<S51>:1:18' */
  /*  MIX = [  0   1  -1   ;      % QUAD [+] */
  /*          -1   0   1   ; */
  /*           0  -1  -1   ; */
  /*           1   0   1   ]; */
  /* ==============================  Forces  ==================================  */
  /*                      We are evaluating in Body frame */
  /*  gravitational force */
  /* '<S51>:1:31' */
  /*  drag force (Square) */
  /* '<S51>:1:33' */
  /* --------------------------------Thrust Model------------------------------ */
  /* '<S51>:1:37' */
  rtb_Switch_i = 6.1235421348947671E-6 * DynModel_B.Memory2;
  rtb_thrust_idx_0 = DynModel_B.DiscreteTransferFcn *
    DynModel_B.DiscreteTransferFcn * rtb_Switch_i;
  rtb_thrust_idx_1 = DynModel_B.DiscreteTransferFcn_o *
    DynModel_B.DiscreteTransferFcn_o * rtb_Switch_i;
  rtb_thrust_idx_2 = DynModel_B.DiscreteTransferFcn_i *
    DynModel_B.DiscreteTransferFcn_i * rtb_Switch_i;
  rtb_thrust = DynModel_B.DiscreteTransferFcn_ib *
    DynModel_B.DiscreteTransferFcn_ib * rtb_Switch_i;

  /* Sum: '<S2>/Add7' incorporates:
   *  Constant: '<S8>/h_ref12'
   *  MATLAB Function: '<S8>/multicopter'
   */
  /*  rotor thrust */
  /* thrust = 5 * rotors; */
  /* -------------------------------------------------------------------------- */
  /* '<S51>:1:41' */
  /* '<S51>:1:42' */
  /* ==================================Moments================================= */
  /*  Thrusts contributions to momentum */
  /* '<S51>:1:47' */
  /*  x moment */
  /* '<S51>:1:48' */
  /*  y moment */
  /*  Torques contributions */
  /*  x rotor momentum */
  /* momentum_x = sum(abs(MIX(:, 3)) .* rotor_inertia .* rotors) * omega(1); */
  /*  y rotor momentum */
  /* momentum_y = sum(abs(MIX(:, 3)) .* rotor_inertia .* rotors) * omega(2); */
  /*  z rotor momentum */
  /* momentum_z = sum(MIX(:, 3) .* rotor_inertia .* rotors) * omega(3);           */
  /* --------------------------------Torque Model------------------------------ */
  /* '<S51>:1:59' */
  /*  rotor torque */
  /* mz = MIX(:,3).*0.05.*rotors; */
  /* -------------------------------------------------------------------------- */
  /* '<S51>:1:63' */
  /*  - [momentum_x; momentum_y; momentum_z]; */
  /* ========================================================================== */
  rtb_Add7_idx_0 += (-(DynModel_B.ubvbwb[0] * fabs(DynModel_B.ubvbwb[0])) *
                     DynModel_B.Memory2 * 0.016813708498984763 +
                     -DynModel_B.ubvbwb[0] * DynModel_B.Memory2 *
                     0.016813708498984763) + DynModel_B.VectorConcatenate[6] *
    9.81 * 1.2;
  rtb_Add7_idx_1 += (-(DynModel_B.ubvbwb[1] * fabs(DynModel_B.ubvbwb[1])) *
                     DynModel_B.Memory2 * 0.018813708498984762 +
                     -DynModel_B.ubvbwb[1] * DynModel_B.Memory2 *
                     0.018813708498984762) + DynModel_B.VectorConcatenate[7] *
    9.81 * 1.2;
  rtb_Add7 = rtb_jxi + (((-(DynModel_B.ubvbwb[2] * fabs(DynModel_B.ubvbwb[2])) *
    DynModel_B.Memory2 * 0.18845573684677208 + -DynModel_B.ubvbwb[2] *
    DynModel_B.Memory2 * 0.18845573684677208) + DynModel_B.VectorConcatenate[8] *
    9.81 * 1.2) - (((rtb_thrust_idx_0 + rtb_thrust_idx_1) + rtb_thrust_idx_2) +
                   rtb_thrust));

  /* Product: '<S4>/Product' incorporates:
   *  Constant: '<S12>/Constant'
   */
  DynModel_B.Product_b[0] = rtb_Add7_idx_0 / 1.2;
  DynModel_B.Product_b[1] = rtb_Add7_idx_1 / 1.2;
  DynModel_B.Product_b[2] = rtb_Add7 / 1.2;
  if (rtmIsMajorTimeStep(DynModel_M)) {
    /* ZeroOrderHold: '<S134>/Zero-Order Hold1' */
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
    /* ZeroOrderHold: '<S134>/Zero-Order Hold2' */
    rtb_ZeroOrderHold2_idx_0 = DynModel_B.MatrixMultiply1[0];
    rtb_ZeroOrderHold2_idx_1 = DynModel_B.MatrixMultiply1[1];
    rtb_ZeroOrderHold2_idx_2 = DynModel_B.MatrixMultiply1[2];
  }

  /* Integrator: '<S4>/p,q,r ' */
  DynModel_B.pqr[0] = DynModel_X.pqr_CSTATE[0];
  DynModel_B.pqr[1] = DynModel_X.pqr_CSTATE[1];
  DynModel_B.pqr[2] = DynModel_X.pqr_CSTATE[2];
  if (rtmIsMajorTimeStep(DynModel_M)) {
    /* ZeroOrderHold: '<S134>/Zero-Order Hold' */
    rtb_Product_k4[0] = DynModel_B.pqr[0];
    rtb_Product_k4[1] = DynModel_B.pqr[1];
    rtb_Product_k4[2] = DynModel_B.pqr[2];
  }

  /* Sqrt: '<S35>/sqrt' incorporates:
   *  Product: '<S36>/Product'
   *  Product: '<S36>/Product1'
   *  Product: '<S36>/Product2'
   *  Product: '<S36>/Product3'
   *  Sum: '<S36>/Sum'
   */
  rtb_Product3_nw = sqrt(((rtb_q0q1q2q3[0] * rtb_q0q1q2q3[0] + rtb_q0q1q2q3[1] *
    rtb_q0q1q2q3[1]) + rtb_q0q1q2q3[2] * rtb_q0q1q2q3[2]) + rtb_q0q1q2q3[3] *
    rtb_q0q1q2q3[3]);

  /* Product: '<S34>/Product' */
  rtb_jxi = rtb_q0q1q2q3[0] / rtb_Product3_nw;

  /* Product: '<S34>/Product1' */
  rtb_kxj = rtb_q0q1q2q3[1] / rtb_Product3_nw;

  /* Product: '<S34>/Product2' */
  rtb_ixj = rtb_q0q1q2q3[2] / rtb_Product3_nw;

  /* Product: '<S34>/Product3' */
  rtb_Product3_nw = rtb_q0q1q2q3[3] / rtb_Product3_nw;

  /* Fcn: '<S18>/fcn1' */
  u0 = (rtb_kxj * rtb_ixj + rtb_jxi * rtb_Product3_nw) * 2.0;

  /* Fcn: '<S18>/fcn2' */
  u1 = ((rtb_jxi * rtb_jxi + rtb_kxj * rtb_kxj) - rtb_ixj * rtb_ixj) -
    rtb_Product3_nw * rtb_Product3_nw;

  /* Fcn: '<S18>/fcn3' */
  rtb_Switch_i = (rtb_kxj * rtb_Product3_nw - rtb_jxi * rtb_ixj) * -2.0;

  /* Trigonometry: '<S18>/trigFcn' */
  if (rtb_Switch_i > 1.0) {
    rtb_Switch_i = 1.0;
  } else {
    if (rtb_Switch_i < -1.0) {
      rtb_Switch_i = -1.0;
    }
  }

  rtb_VectorConcatenate_idx_1 = asin(rtb_Switch_i);

  /* End of Trigonometry: '<S18>/trigFcn' */

  /* Fcn: '<S18>/fcn4' */
  rtb_Switch_i = (rtb_ixj * rtb_Product3_nw + rtb_jxi * rtb_kxj) * 2.0;

  /* Fcn: '<S18>/fcn5' */
  rtb_jxi = ((rtb_jxi * rtb_jxi - rtb_kxj * rtb_kxj) - rtb_ixj * rtb_ixj) +
    rtb_Product3_nw * rtb_Product3_nw;

  /* Trigonometry: '<S18>/Trigonometric Function3' */
  rtb_ixj = atan2(rtb_Switch_i, rtb_jxi);
  if (rtmIsMajorTimeStep(DynModel_M)) {
    /* Gain: '<S134>/Gain' incorporates:
     *  Constant: '<S134>/wl_ins'
     *  Constant: '<S3>/center of gravity'
     *  Sum: '<S134>/Sum7'
     */
    rtb_Saturation_l_idx_0 = 0.0;
    rtb_Saturation_l_idx_1 = 0.0;
    rtb_Saturation_l_idx_2 = 0.0;

    /* Sum: '<S145>/Sum' incorporates:
     *  Product: '<S147>/i x j'
     *  Product: '<S147>/j x k'
     *  Product: '<S147>/k x i'
     *  Product: '<S148>/i x k'
     *  Product: '<S148>/j x i'
     *  Product: '<S148>/k x j'
     */
    rtb_Saturation_kg[0] = 0.0;
    rtb_Saturation_kg[1] = 0.0;
    rtb_Saturation_kg[2] = 0.0;
  }

  /* Switch: '<S5>/Switch1' incorporates:
   *  Gain: '<S5>/Gain'
   *  Gain: '<S5>/Gain4'
   *  Gain: '<S5>/Gain5'
   *  Sum: '<S5>/Add2'
   */
  if (-DynModel_B.Sum1 >= 0.0) {
    rtb_Product3_nw = Kwreact * DynModel_B.pqr[0] + -rtb_ixj;
    rtb_Add6_idx_1 = Kwreact * DynModel_B.pqr[1] + -rtb_VectorConcatenate_idx_1;
    rtb_Switch_i = Kwreact * DynModel_B.pqr[2];
  } else {
    rtb_Product3_nw = 0.0;
    rtb_Add6_idx_1 = 0.0;
    rtb_Switch_i = 0.0;
  }

  /* End of Switch: '<S5>/Switch1' */

  /* Sum: '<S2>/Add6' incorporates:
   *  Constant: '<S8>/h_ref10'
   *  Gain: '<S2>/Gain'
   *  MATLAB Function: '<S8>/multicopter'
   *  SignalConversion: '<S51>/TmpSignal ConversionAt SFunction Inport4'
   */
  rtb_Product3_nw = ((((rtb_thrust_idx_0 * 0.2 * 1.4142135623730951 / 2.0 +
                        -rtb_thrust_idx_1 * 0.2 * 1.4142135623730951 / 2.0) +
                       -rtb_thrust_idx_2 * 0.2 * 1.4142135623730951 / 2.0) +
                      rtb_thrust * 0.2 * 1.4142135623730951 / 2.0) +
                     rtb_Product3_nw) - 0.1 * DynModel_B.pqr[0];
  rtb_Add6_idx_1 = ((((rtb_thrust_idx_0 * 0.2 * 1.4142135623730951 / 2.0 +
                       rtb_thrust_idx_1 * 0.2 * 1.4142135623730951 / 2.0) +
                      -rtb_thrust_idx_2 * 0.2 * 1.4142135623730951 / 2.0) +
                     -rtb_thrust * 0.2 * 1.4142135623730951 / 2.0) +
                    rtb_Add6_idx_1) - 0.1 * DynModel_B.pqr[1];
  rtb_kxj = ((((-3.564683251291932E-8 * DynModel_B.Memory2 *
                (DynModel_B.DiscreteTransferFcn * DynModel_B.DiscreteTransferFcn)
                + 3.564683251291932E-8 * DynModel_B.Memory2 *
                (DynModel_B.DiscreteTransferFcn_o *
                 DynModel_B.DiscreteTransferFcn_o)) + -3.564683251291932E-8 *
               DynModel_B.Memory2 * (DynModel_B.DiscreteTransferFcn_i *
    DynModel_B.DiscreteTransferFcn_i)) + 3.564683251291932E-8 *
              DynModel_B.Memory2 * (DynModel_B.DiscreteTransferFcn_ib *
    DynModel_B.DiscreteTransferFcn_ib)) + rtb_Switch_i) - 0.1 * DynModel_B.pqr[2];

  /* Product: '<S37>/Product' */
  for (rtb_Compare_0 = 0; rtb_Compare_0 < 3; rtb_Compare_0++) {
    rtb_forces[rtb_Compare_0] = DynModel_ConstB.Selector[rtb_Compare_0 + 6] *
      DynModel_B.pqr[2] + (DynModel_ConstB.Selector[rtb_Compare_0 + 3] *
      DynModel_B.pqr[1] + DynModel_ConstB.Selector[rtb_Compare_0] *
      DynModel_B.pqr[0]);
  }

  /* End of Product: '<S37>/Product' */

  /* Product: '<S11>/Product2' incorporates:
   *  Product: '<S40>/i x j'
   *  Product: '<S40>/j x k'
   *  Product: '<S40>/k x i'
   *  Product: '<S41>/i x k'
   *  Product: '<S41>/j x i'
   *  Product: '<S41>/k x j'
   *  Sum: '<S11>/Sum2'
   *  Sum: '<S39>/Sum'
   */
  rtb_Add6_0[0] = rtb_Product3_nw - (DynModel_B.pqr[1] * rtb_forces[2] -
    DynModel_B.pqr[2] * rtb_forces[1]);
  rtb_Add6_0[1] = rtb_Add6_idx_1 - (DynModel_B.pqr[2] * rtb_forces[0] -
    DynModel_B.pqr[0] * rtb_forces[2]);
  rtb_Add6_0[2] = rtb_kxj - (DynModel_B.pqr[0] * rtb_forces[1] - DynModel_B.pqr
    [1] * rtb_forces[0]);
  rt_mrdivide_U1d1x3_U2d3x3_Yd1x3(rtb_Add6_0, DynModel_ConstB.Selector2,
    DynModel_B.Product2_m);
  if (rtmIsMajorTimeStep(DynModel_M)) {
    /* Sum: '<S134>/Sum' */
    rtb_Sum_h_idx_0 = (rtb_Sum_h_idx_0 - rtb_ZeroOrderHold2_idx_0) +
      rtb_Saturation_kg[0];
    rtb_Sum_h_idx_1 = (rtb_Sum_h_idx_1 - rtb_ZeroOrderHold2_idx_1) +
      rtb_Saturation_kg[1];
    rtb_Switch_i = (rtb_Sum_h_idx_2 - rtb_ZeroOrderHold2_idx_2) +
      rtb_Saturation_kg[2];

    /* Product: '<S134>/Product' incorporates:
     *  Constant: '<S134>/Scale factors & Cross-coupling  errors'
     */
    for (rtb_Compare_0 = 0; rtb_Compare_0 < 3; rtb_Compare_0++) {
      rtb_Saturation_kg[rtb_Compare_0] = DynModel_ConstP.pooled23[rtb_Compare_0
        + 6] * rtb_Switch_i + (DynModel_ConstP.pooled23[rtb_Compare_0 + 3] *
        rtb_Sum_h_idx_1 + DynModel_ConstP.pooled23[rtb_Compare_0] *
        rtb_Sum_h_idx_0);
    }

    /* End of Product: '<S134>/Product' */

    /* Sum: '<S134>/Sum4' */
    rtb_forces[1] = rtb_Saturation_kg[1];
    rtb_forces[2] = rtb_Saturation_kg[2];

    /* Saturate: '<S134>/Saturation' incorporates:
     *  Gain: '<S137>/Output'
     *  RandomNumber: '<S137>/White Noise'
     *  Sum: '<S134>/Sum1'
     *  Sum: '<S134>/Sum4'
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
      rtb_forces[1];
    if (rtb_Switch_i > 19.62) {
      rtb_Saturation_kg[1] = 19.62;
    } else if (rtb_Switch_i < -19.62) {
      rtb_Saturation_kg[1] = -19.62;
    } else {
      rtb_Saturation_kg[1] = rtb_Switch_i;
    }

    rtb_Switch_i = 0.011180339887498949 * DynModel_DW.NextOutput_k[2] +
      rtb_forces[2];
    if (rtb_Switch_i > 19.62) {
      rtb_Saturation_kg[2] = 19.62;
    } else if (rtb_Switch_i < -19.62) {
      rtb_Saturation_kg[2] = -19.62;
    } else {
      rtb_Saturation_kg[2] = rtb_Switch_i;
    }

    /* End of Saturate: '<S134>/Saturation' */

    /* ZeroOrderHold: '<S135>/Zero-Order Hold' */
    rtb_Saturation_l_idx_0 = DynModel_B.pqr[0];
    rtb_Saturation_l_idx_1 = DynModel_B.pqr[1];
    rtb_Saturation_l_idx_2 = DynModel_B.pqr[2];

    /* Product: '<S135>/Product' incorporates:
     *  Constant: '<S135>/Scale factors & Cross-coupling  errors '
     *  ZeroOrderHold: '<S135>/Zero-Order Hold'
     */
    for (rtb_Compare_0 = 0; rtb_Compare_0 < 3; rtb_Compare_0++) {
      rtb_Product_k4[rtb_Compare_0] = DynModel_ConstP.pooled23[rtb_Compare_0 + 6]
        * DynModel_B.pqr[2] + (DynModel_ConstP.pooled23[rtb_Compare_0 + 3] *
        DynModel_B.pqr[1] + DynModel_ConstP.pooled23[rtb_Compare_0] *
        DynModel_B.pqr[0]);
    }

    /* End of Product: '<S135>/Product' */
  }

  /* Gain: '<S133>/Unit Conversion' */
  DynModel_B.UnitConversion[0] = 0.10197162129779283 * DynModel_B.Product_b[0];
  DynModel_B.UnitConversion[1] = 0.10197162129779283 * DynModel_B.Product_b[1];
  DynModel_B.UnitConversion[2] = 0.10197162129779283 * DynModel_B.Product_b[2];
  if (rtmIsMajorTimeStep(DynModel_M)) {
    /* Saturate: '<S135>/Saturation' incorporates:
     *  Gain: '<S154>/Output'
     *  RandomNumber: '<S154>/White Noise'
     *  Sum: '<S135>/Sum1'
     *  Sum: '<S135>/Sum4'
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

    /* End of Saturate: '<S135>/Saturation' */
  }

  /* Sum: '<S73>/Add' */
  rtb_Switch_i = (DynModel_B.VectorConcatenate[0] +
                  DynModel_B.VectorConcatenate[4]) +
    DynModel_B.VectorConcatenate[8];

  /* If: '<S60>/If' incorporates:
   *  Sum: '<S73>/Add'
   */
  if (rtmIsMajorTimeStep(DynModel_M)) {
    DynModel_DW.If_ActiveSubsystem = (int8_T)!(rtb_Switch_i > 0.0);
  }

  switch (DynModel_DW.If_ActiveSubsystem) {
   case 0L:
    /* Outputs for IfAction SubSystem: '<S60>/Positive Trace' incorporates:
     *  ActionPort: '<S72>/Action Port'
     */
    /* Sqrt: '<S72>/sqrt' incorporates:
     *  Constant: '<S72>/Constant'
     *  Sum: '<S72>/Sum'
     *  Sum: '<S73>/Add'
     */
    rtb_Switch_i = sqrt(rtb_Switch_i + 1.0);

    /* Gain: '<S72>/Gain' */
    DynModel_B.Merge[0] = 0.5 * rtb_Switch_i;

    /* Gain: '<S72>/Gain1' */
    rtb_Switch_i *= 2.0;

    /* Product: '<S72>/Product' incorporates:
     *  Sum: '<S93>/Add'
     *  Sum: '<S94>/Add'
     *  Sum: '<S95>/Add'
     */
    DynModel_B.Merge[1] = (DynModel_B.VectorConcatenate[7] -
      DynModel_B.VectorConcatenate[5]) / rtb_Switch_i;
    DynModel_B.Merge[2] = (DynModel_B.VectorConcatenate[2] -
      DynModel_B.VectorConcatenate[6]) / rtb_Switch_i;
    DynModel_B.Merge[3] = (DynModel_B.VectorConcatenate[3] -
      DynModel_B.VectorConcatenate[1]) / rtb_Switch_i;

    /* End of Outputs for SubSystem: '<S60>/Positive Trace' */
    break;

   case 1L:
    /* Outputs for IfAction SubSystem: '<S60>/Negative Trace' incorporates:
     *  ActionPort: '<S71>/Action Port'
     */
    /* If: '<S71>/Find Maximum Diagonal Value' */
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
      /* Outputs for IfAction SubSystem: '<S71>/Maximum Value at DCM(2,2)' incorporates:
       *  ActionPort: '<S75>/Action Port'
       */
      /* Sqrt: '<S75>/sqrt' incorporates:
       *  Constant: '<S87>/Constant'
       *  Sum: '<S87>/Add'
       */
      rtb_Switch_i = sqrt(((DynModel_B.VectorConcatenate[4] -
                            DynModel_B.VectorConcatenate[0]) -
                           DynModel_B.VectorConcatenate[8]) + 1.0);

      /* Gain: '<S75>/Gain' */
      DynModel_B.Merge[2] = 0.5 * rtb_Switch_i;

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

      /* Gain: '<S75>/Gain1' incorporates:
       *  Product: '<S75>/Product'
       *  Sum: '<S85>/Add'
       */
      DynModel_B.Merge[1] = (DynModel_B.VectorConcatenate[1] +
        DynModel_B.VectorConcatenate[3]) * rtb_Switch_i;

      /* Gain: '<S75>/Gain3' incorporates:
       *  Product: '<S75>/Product'
       *  Sum: '<S84>/Add'
       */
      DynModel_B.Merge[3] = (DynModel_B.VectorConcatenate[5] +
        DynModel_B.VectorConcatenate[7]) * rtb_Switch_i;

      /* Gain: '<S75>/Gain4' incorporates:
       *  Product: '<S75>/Product'
       *  Sum: '<S83>/Add'
       */
      DynModel_B.Merge[0] = (DynModel_B.VectorConcatenate[2] -
        DynModel_B.VectorConcatenate[6]) * rtb_Switch_i;

      /* End of Outputs for SubSystem: '<S71>/Maximum Value at DCM(2,2)' */
      break;

     case 1L:
      /* Outputs for IfAction SubSystem: '<S71>/Maximum Value at DCM(3,3)' incorporates:
       *  ActionPort: '<S76>/Action Port'
       */
      /* Sqrt: '<S76>/sqrt' incorporates:
       *  Constant: '<S92>/Constant'
       *  Sum: '<S92>/Add'
       */
      rtb_Switch_i = sqrt(((DynModel_B.VectorConcatenate[8] -
                            DynModel_B.VectorConcatenate[0]) -
                           DynModel_B.VectorConcatenate[4]) + 1.0);

      /* Gain: '<S76>/Gain' */
      DynModel_B.Merge[3] = 0.5 * rtb_Switch_i;

      /* Switch: '<S91>/Switch' incorporates:
       *  Constant: '<S91>/Constant1'
       */
      if (rtb_Switch_i != 0.0) {
        rtb_jxi = 0.5;
      } else {
        rtb_jxi = 0.0;
        rtb_Switch_i = 1.0;
      }

      /* End of Switch: '<S91>/Switch' */

      /* Product: '<S91>/Product' */
      rtb_Switch_i = rtb_jxi / rtb_Switch_i;

      /* Gain: '<S76>/Gain1' incorporates:
       *  Product: '<S76>/Product'
       *  Sum: '<S88>/Add'
       */
      DynModel_B.Merge[1] = (DynModel_B.VectorConcatenate[2] +
        DynModel_B.VectorConcatenate[6]) * rtb_Switch_i;

      /* Gain: '<S76>/Gain2' incorporates:
       *  Product: '<S76>/Product'
       *  Sum: '<S89>/Add'
       */
      DynModel_B.Merge[2] = (DynModel_B.VectorConcatenate[5] +
        DynModel_B.VectorConcatenate[7]) * rtb_Switch_i;

      /* Gain: '<S76>/Gain3' incorporates:
       *  Product: '<S76>/Product'
       *  Sum: '<S90>/Add'
       */
      DynModel_B.Merge[0] = (DynModel_B.VectorConcatenate[3] -
        DynModel_B.VectorConcatenate[1]) * rtb_Switch_i;

      /* End of Outputs for SubSystem: '<S71>/Maximum Value at DCM(3,3)' */
      break;

     case 2L:
      /* Outputs for IfAction SubSystem: '<S71>/Maximum Value at DCM(1,1)' incorporates:
       *  ActionPort: '<S74>/Action Port'
       */
      /* Sqrt: '<S74>/sqrt' incorporates:
       *  Constant: '<S82>/Constant'
       *  Sum: '<S82>/Add'
       */
      rtb_Switch_i = sqrt(((DynModel_B.VectorConcatenate[0] -
                            DynModel_B.VectorConcatenate[4]) -
                           DynModel_B.VectorConcatenate[8]) + 1.0);

      /* Gain: '<S74>/Gain' */
      DynModel_B.Merge[1] = 0.5 * rtb_Switch_i;

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

      /* Gain: '<S74>/Gain1' incorporates:
       *  Product: '<S74>/Product'
       *  Sum: '<S80>/Add'
       */
      DynModel_B.Merge[2] = (DynModel_B.VectorConcatenate[1] +
        DynModel_B.VectorConcatenate[3]) * rtb_Switch_i;

      /* Gain: '<S74>/Gain2' incorporates:
       *  Product: '<S74>/Product'
       *  Sum: '<S78>/Add'
       */
      DynModel_B.Merge[3] = (DynModel_B.VectorConcatenate[2] +
        DynModel_B.VectorConcatenate[6]) * rtb_Switch_i;

      /* Gain: '<S74>/Gain3' incorporates:
       *  Product: '<S74>/Product'
       *  Sum: '<S79>/Add'
       */
      DynModel_B.Merge[0] = (DynModel_B.VectorConcatenate[7] -
        DynModel_B.VectorConcatenate[5]) * rtb_Switch_i;

      /* End of Outputs for SubSystem: '<S71>/Maximum Value at DCM(1,1)' */
      break;
    }

    /* End of If: '<S71>/Find Maximum Diagonal Value' */
    /* End of Outputs for SubSystem: '<S60>/Negative Trace' */
    break;
  }

  /* End of If: '<S60>/If' */

  /* DotProduct: '<S20>/Dot Product' */
  rtb_Switch_i = ((rtb_q0q1q2q3[0] * rtb_q0q1q2q3[0] + rtb_q0q1q2q3[1] *
                   rtb_q0q1q2q3[1]) + rtb_q0q1q2q3[2] * rtb_q0q1q2q3[2]) +
    rtb_q0q1q2q3[3] * rtb_q0q1q2q3[3];

  /* Fcn: '<S20>/q0dot' incorporates:
   *  Constant: '<S20>/Constant'
   *  DotProduct: '<S20>/Dot Product'
   *  Sum: '<S20>/Sum'
   */
  DynModel_B.q0dot = ((rtb_q0q1q2q3[1] * DynModel_B.pqr[0] + rtb_q0q1q2q3[2] *
                       DynModel_B.pqr[1]) + rtb_q0q1q2q3[3] * DynModel_B.pqr[2])
    * -0.5 + (1.0 - rtb_Switch_i) * rtb_q0q1q2q3[0];

  /* Fcn: '<S20>/q1dot' incorporates:
   *  Constant: '<S20>/Constant'
   *  DotProduct: '<S20>/Dot Product'
   *  Sum: '<S20>/Sum'
   */
  DynModel_B.q1dot = ((rtb_q0q1q2q3[0] * DynModel_B.pqr[0] + rtb_q0q1q2q3[2] *
                       DynModel_B.pqr[2]) - rtb_q0q1q2q3[3] * DynModel_B.pqr[1])
    * 0.5 + (1.0 - rtb_Switch_i) * rtb_q0q1q2q3[1];

  /* Fcn: '<S20>/q2dot' incorporates:
   *  Constant: '<S20>/Constant'
   *  DotProduct: '<S20>/Dot Product'
   *  Sum: '<S20>/Sum'
   */
  DynModel_B.q2dot = ((rtb_q0q1q2q3[0] * DynModel_B.pqr[1] + rtb_q0q1q2q3[3] *
                       DynModel_B.pqr[0]) - rtb_q0q1q2q3[1] * DynModel_B.pqr[2])
    * 0.5 + (1.0 - rtb_Switch_i) * rtb_q0q1q2q3[2];

  /* Fcn: '<S20>/q3dot' incorporates:
   *  Constant: '<S20>/Constant'
   *  DotProduct: '<S20>/Dot Product'
   *  Sum: '<S20>/Sum'
   */
  DynModel_B.q3dot = ((rtb_q0q1q2q3[0] * DynModel_B.pqr[2] + rtb_q0q1q2q3[1] *
                       DynModel_B.pqr[1]) - rtb_q0q1q2q3[2] * DynModel_B.pqr[0])
    * 0.5 + (1.0 - rtb_Switch_i) * rtb_q0q1q2q3[3];

  /* Sum: '<S4>/Sum' incorporates:
   *  Product: '<S42>/i x j'
   *  Product: '<S42>/j x k'
   *  Product: '<S42>/k x i'
   *  Product: '<S43>/i x k'
   *  Product: '<S43>/j x i'
   *  Product: '<S43>/k x j'
   *  Sum: '<S13>/Sum'
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
    rtb_DataTypeConversion14_idx_0 = (real32_T)rtb_Saturation_kg[0];
    rtb_DataTypeConversion14_idx_1 = (real32_T)rtb_Saturation_kg[1];
    rtb_DataTypeConversion14_idx_2 = (real32_T)rtb_Saturation_kg[2];

    /* DataTypeConversion: '<S1>/Data Type Conversion15' */
    rtb_DataTypeConversion15_idx_0 = (real32_T)rtb_Saturation_l_idx_0;
    rtb_DataTypeConversion15_idx_1 = (real32_T)rtb_Saturation_l_idx_1;
    rtb_DataTypeConversion15_idx_2 = (real32_T)rtb_Saturation_l_idx_2;

    /* DataTypeConversion: '<S1>/Data Type Conversion21' */
    rtb_DataTypeConversion21_idx_0 = (real32_T)DynModel_B.DiscreteTransferFcn;
    rtb_DataTypeConversion21_idx_1 = (real32_T)DynModel_B.DiscreteTransferFcn_o;
    rtb_DataTypeConversion21_idx_2 = (real32_T)DynModel_B.DiscreteTransferFcn_i;
    rtb_DataTypeConversion21_idx_3 = (real32_T)DynModel_B.DiscreteTransferFcn_ib;

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
    DynModel_Y.Gps_V[0] = (real32_T)rtb_Add2_h_idx_0;
    DynModel_Y.Gps_V[1] = (real32_T)rtb_Add2_h_idx_1;
    DynModel_Y.Gps_V[2] = (real32_T)rtb_Add2_h_idx_2;

    /* Outport: '<Root>/Gps_V_Mod' incorporates:
     *  DataTypeConversion: '<S1>/Data Type Conversion18'
     *  DotProduct: '<S1>/Dot Product'
     *  Sqrt: '<S1>/Sqrt'
     */
    DynModel_Y.Gps_V_Mod = (real32_T)sqrt((rtb_Add2_h_idx_0 * rtb_Add2_h_idx_0 +
      rtb_Add2_h_idx_1 * rtb_Add2_h_idx_1) + rtb_Add2_h_idx_2 * rtb_Add2_h_idx_2);

    /* Outport: '<Root>/COG' incorporates:
     *  DataTypeConversion: '<S1>/Data Type Conversion9'
     */
    DynModel_Y.COG = (real32_T)rtb_TrigonometricFunction;
  }

  /* Outport: '<Root>/Lat_Lon_Alt' incorporates:
   *  DataTypeConversion: '<S1>/Data Type Conversion10'
   */
  DynModel_Y.Lat_Lon_Alt[0] = (real32_T)rtb_Switch_d;

  /* Switch: '<S101>/Switch' incorporates:
   *  Abs: '<S101>/Abs'
   */
  if (fabs(rtb_Sum_e) > 180.0) {
    /* Outport: '<Root>/Lat_Lon_Alt' incorporates:
     *  Bias: '<S101>/Bias'
     *  Bias: '<S101>/Bias1'
     *  Constant: '<S101>/Constant2'
     *  DataTypeConversion: '<S1>/Data Type Conversion10'
     *  Math: '<S101>/Math Function1'
     */
    DynModel_Y.Lat_Lon_Alt[1] = (real32_T)(rt_modd(rtb_Sum_e + 180.0, 360.0) +
      -180.0);
  } else {
    /* Outport: '<Root>/Lat_Lon_Alt' incorporates:
     *  DataTypeConversion: '<S1>/Data Type Conversion10'
     */
    DynModel_Y.Lat_Lon_Alt[1] = (real32_T)rtb_Sum_e;
  }

  /* End of Switch: '<S101>/Switch' */

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
   *  Trigonometry: '<S18>/Trigonometric Function1'
   */
  DynModel_Y.RPY[0] = (real32_T)(rtb_ixj + DynModel_B.Output);
  DynModel_Y.RPY[1] = (real32_T)(rtb_VectorConcatenate_idx_1 + DynModel_B.Output);
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

  /* Outport: '<Root>/Forces' incorporates:
   *  DataTypeConversion: '<S1>/Data Type Conversion2'
   */
  DynModel_Y.Forces[0] = (real32_T)rtb_Add7_idx_0;
  DynModel_Y.Forces[1] = (real32_T)rtb_Add7_idx_1;
  DynModel_Y.Forces[2] = (real32_T)rtb_Add7;

  /* Outport: '<Root>/Torques' incorporates:
   *  DataTypeConversion: '<S1>/Data Type Conversion19'
   */
  DynModel_Y.Torques[0] = (real32_T)rtb_Product3_nw;
  DynModel_Y.Torques[1] = (real32_T)rtb_Add6_idx_1;
  DynModel_Y.Torques[2] = (real32_T)rtb_kxj;

  /* Outport: '<Root>/Thrusts' incorporates:
   *  DataTypeConversion: '<S1>/Data Type Conversion20'
   */
  DynModel_Y.Thrusts[0] = (real32_T)rtb_thrust_idx_0;
  DynModel_Y.Thrusts[1] = (real32_T)rtb_thrust_idx_1;
  DynModel_Y.Thrusts[2] = (real32_T)rtb_thrust_idx_2;
  DynModel_Y.Thrusts[3] = (real32_T)rtb_thrust;
  if (rtmIsMajorTimeStep(DynModel_M)) {
    /* Outport: '<Root>/Rotor_Speed' */
    DynModel_Y.Rotor_Speed[0] = rtb_DataTypeConversion21_idx_0;
    DynModel_Y.Rotor_Speed[1] = rtb_DataTypeConversion21_idx_1;
    DynModel_Y.Rotor_Speed[2] = rtb_DataTypeConversion21_idx_2;
    DynModel_Y.Rotor_Speed[3] = rtb_DataTypeConversion21_idx_3;
  }

  /* Outport: '<Root>/Xe' incorporates:
   *  DataTypeConversion: '<S1>/Data Type Conversion22'
   */
  DynModel_Y.Xe[0] = (real32_T)DynModel_B.xeyeze[0];
  DynModel_Y.Xe[1] = (real32_T)DynModel_B.xeyeze[1];
  DynModel_Y.Xe[2] = (real32_T)DynModel_B.xeyeze[2];
  if (rtmIsMajorTimeStep(DynModel_M)) {
    if (rtmIsMajorTimeStep(DynModel_M)) {
      /* Update for RandomNumber: '<S70>/Random Number' */
      DynModel_DW.NextOutput = rt_nrand_Upu32_Yd_f_pw(&DynModel_DW.RandSeed) *
        0.0031622776601683794;

      /* Update for RandomNumber: '<S68>/Random Number' */
      DynModel_DW.NextOutput_a = rt_nrand_Upu32_Yd_f_pw(&DynModel_DW.RandSeed_f)
        * 0.001;

      /* Update for RandomNumber: '<S66>/Random Number' */
      DynModel_DW.NextOutput_l = rt_nrand_Upu32_Yd_f_pw(&DynModel_DW.RandSeed_fw)
        * 0.001;

      /* Update for RandomNumber: '<S58>/Random Number' */
      DynModel_DW.NextOutput_n = rt_nrand_Upu32_Yd_f_pw(&DynModel_DW.RandSeed_fm)
        * 0.0031622776601683794;

      /* Update for RandomNumber: '<S62>/Random Number' */
      DynModel_DW.NextOutput_o[0] = rt_nrand_Upu32_Yd_f_pw
        (&DynModel_DW.RandSeed_e[0]) * 0.1;
      DynModel_DW.NextOutput_o[1] = rt_nrand_Upu32_Yd_f_pw
        (&DynModel_DW.RandSeed_e[1]) * 0.1;
      DynModel_DW.NextOutput_o[2] = rt_nrand_Upu32_Yd_f_pw
        (&DynModel_DW.RandSeed_e[2]) * 0.1;

      /* Update for RandomNumber: '<S62>/Random Number1' */
      DynModel_DW.NextOutput_h[0] = rt_nrand_Upu32_Yd_f_pw
        (&DynModel_DW.RandSeed_i[0]) * 0.00031622776601683794;
      DynModel_DW.NextOutput_h[1] = rt_nrand_Upu32_Yd_f_pw
        (&DynModel_DW.RandSeed_i[1]) * 0.00031622776601683794;
      DynModel_DW.NextOutput_h[2] = rt_nrand_Upu32_Yd_f_pw
        (&DynModel_DW.RandSeed_i[2]) * 0.00031622776601683794;
    }

    /* Update for Integrator: '<S10>/q0 q1 q2 q3' */
    DynModel_DW.q0q1q2q3_IWORK.IcNeedsLoading = 0;
    if (rtmIsMajorTimeStep(DynModel_M)) {
      /* Update for RandomNumber: '<S65>/Random Number' */
      DynModel_DW.NextOutput_am = rt_nrand_Upu32_Yd_f_pw(&DynModel_DW.RandSeed_p)
        * 0.001;

      /* Update for RandomNumber: '<S59>/White Noise' */
      DynModel_DW.NextOutput_lh = rt_nrand_Upu32_Yd_f_pw(&DynModel_DW.RandSeed_l);

      /* Update for UnitDelay: '<S44>/UD' */
      DynModel_DW.UD_DSTATE = rtb_TSamp;

      /* Update for Memory: '<S2>/Memory2' */
      DynModel_DW.Memory2_PreviousInput = DynModel_B.Product3;

      /* Update for DiscreteTransferFcn: '<S45>/Discrete Transfer Fcn' */
      DynModel_DW.DiscreteTransferFcn_states =
        DynModel_DW.DiscreteTransferFcn_tmp;

      /* Update for DiscreteTransferFcn: '<S46>/Discrete Transfer Fcn' */
      DynModel_DW.DiscreteTransferFcn_states_m =
        DynModel_DW.DiscreteTransferFcn_tmp_i;

      /* Update for DiscreteTransferFcn: '<S47>/Discrete Transfer Fcn' */
      DynModel_DW.DiscreteTransferFcn_states_c =
        DynModel_DW.DiscreteTransferFcn_tmp_b;

      /* Update for DiscreteTransferFcn: '<S48>/Discrete Transfer Fcn' */
      DynModel_DW.DiscreteTransferFcn_states_a =
        DynModel_DW.DiscreteTransferFcn_tmp_e;

      /* Update for RandomNumber: '<S137>/White Noise' */
      DynModel_DW.NextOutput_k[0] = rt_nrand_Upu32_Yd_f_pw
        (&DynModel_DW.RandSeed_ls[0]);
      DynModel_DW.NextOutput_k[1] = rt_nrand_Upu32_Yd_f_pw
        (&DynModel_DW.RandSeed_ls[1]);
      DynModel_DW.NextOutput_k[2] = rt_nrand_Upu32_Yd_f_pw
        (&DynModel_DW.RandSeed_ls[2]);

      /* Update for RandomNumber: '<S154>/White Noise' */
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

  /* Derivatives for Integrator: '<S10>/q0 q1 q2 q3' */
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

  /* Start for If: '<S60>/If' */
  DynModel_DW.If_ActiveSubsystem = -1;

  /* Start for IfAction SubSystem: '<S60>/Negative Trace' */
  /* Start for If: '<S71>/Find Maximum Diagonal Value' */
  DynModel_DW.FindMaximumDiagonalValue_ActiveSubsystem = -1;

  /* End of Start for SubSystem: '<S60>/Negative Trace' */

  /* Start for Merge: '<S60>/Merge' */
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

    /* InitializeConditions for RandomNumber: '<S70>/Random Number' */
    DynModel_DW.RandSeed = 1144108930UL;
    DynModel_DW.NextOutput = rt_nrand_Upu32_Yd_f_pw(&DynModel_DW.RandSeed) *
      0.0031622776601683794;

    /* InitializeConditions for RandomNumber: '<S68>/Random Number' */
    DynModel_DW.RandSeed_f = 1144108930UL;
    DynModel_DW.NextOutput_a = rt_nrand_Upu32_Yd_f_pw(&DynModel_DW.RandSeed_f) *
      0.001;

    /* InitializeConditions for Integrator: '<S4>/ub,vb,wb' */
    DynModel_X.ubvbwb_CSTATE[0] = 0.0;
    DynModel_X.ubvbwb_CSTATE[1] = 0.0;
    DynModel_X.ubvbwb_CSTATE[2] = 0.0;

    /* InitializeConditions for RandomNumber: '<S66>/Random Number' */
    DynModel_DW.RandSeed_fw = 1144108930UL;
    DynModel_DW.NextOutput_l = rt_nrand_Upu32_Yd_f_pw(&DynModel_DW.RandSeed_fw) *
      0.001;

    /* InitializeConditions for RandomNumber: '<S58>/Random Number' */
    DynModel_DW.RandSeed_fm = 1144108930UL;
    DynModel_DW.NextOutput_n = rt_nrand_Upu32_Yd_f_pw(&DynModel_DW.RandSeed_fm) *
      0.0031622776601683794;

    /* InitializeConditions for RandomNumber: '<S62>/Random Number' */
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

    /* InitializeConditions for RandomNumber: '<S62>/Random Number1' */
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

    /* InitializeConditions for Integrator: '<S10>/q0 q1 q2 q3' */
    if (rtmIsFirstInitCond(DynModel_M)) {
      DynModel_X.q0q1q2q3_CSTATE[0] = 0.0;
      DynModel_X.q0q1q2q3_CSTATE[1] = 0.0;
      DynModel_X.q0q1q2q3_CSTATE[2] = 0.0;
      DynModel_X.q0q1q2q3_CSTATE[3] = 0.0;
    }

    DynModel_DW.q0q1q2q3_IWORK.IcNeedsLoading = 1;

    /* InitializeConditions for RandomNumber: '<S65>/Random Number' */
    DynModel_DW.RandSeed_p = 1144108930UL;
    DynModel_DW.NextOutput_am = rt_nrand_Upu32_Yd_f_pw(&DynModel_DW.RandSeed_p) *
      0.001;

    /* InitializeConditions for RandomNumber: '<S59>/White Noise' */
    DynModel_DW.RandSeed_l = 931168259UL;
    DynModel_DW.NextOutput_lh = rt_nrand_Upu32_Yd_f_pw(&DynModel_DW.RandSeed_l);

    /* InitializeConditions for UnitDelay: '<S44>/UD' */
    DynModel_DW.UD_DSTATE = 0.0;

    /* InitializeConditions for Memory: '<S2>/Memory2' */
    DynModel_DW.Memory2_PreviousInput = 0.0;

    /* InitializeConditions for DiscreteTransferFcn: '<S45>/Discrete Transfer Fcn' */
    DynModel_DW.DiscreteTransferFcn_states = 0.0;

    /* InitializeConditions for DiscreteTransferFcn: '<S46>/Discrete Transfer Fcn' */
    DynModel_DW.DiscreteTransferFcn_states_m = 0.0;

    /* InitializeConditions for DiscreteTransferFcn: '<S47>/Discrete Transfer Fcn' */
    DynModel_DW.DiscreteTransferFcn_states_c = 0.0;

    /* InitializeConditions for DiscreteTransferFcn: '<S48>/Discrete Transfer Fcn' */
    DynModel_DW.DiscreteTransferFcn_states_a = 0.0;

    /* InitializeConditions for Integrator: '<S4>/p,q,r ' */
    DynModel_X.pqr_CSTATE[0] = 0.0;
    DynModel_X.pqr_CSTATE[1] = 0.0;
    DynModel_X.pqr_CSTATE[2] = 0.0;

    /* InitializeConditions for RandomNumber: '<S137>/White Noise' */
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

    /* InitializeConditions for RandomNumber: '<S154>/White Noise' */
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
