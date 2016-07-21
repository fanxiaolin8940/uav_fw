/*
 * SYSTEM0.c
 *
 * Code generation for model "SYSTEM0".
 *
 * Model version              : 1.311
 * Simulink Coder version : 8.8 (R2015a) 09-Feb-2015
 * C source code generated on : Wed Jul 20 23:23:46 2016
 *
 * Target selection: grt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64
 * Emulation hardware selection:
 *    Differs from embedded hardware (Generic->16-bit Embedded Processor)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "SYSTEM0.h"
#include "SYSTEM0_private.h"

/* Block signals (auto storage) */
B_SYSTEM0_T SYSTEM0_B;

/* Continuous states */
X_SYSTEM0_T SYSTEM0_X;

/* Block states (auto storage) */
DW_SYSTEM0_T SYSTEM0_DW;

/* External inputs (root inport signals with auto storage) */
ExtU_SYSTEM0_T SYSTEM0_U;

/* External outputs (root outports fed by signals with auto storage) */
ExtY_SYSTEM0_T SYSTEM0_Y;

/* Real-time model */
RT_MODEL_SYSTEM0_T SYSTEM0_M_;
RT_MODEL_SYSTEM0_T *const SYSTEM0_M = &SYSTEM0_M_;

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
  SYSTEM0_derivatives();

  /* f1 = f(t + (h/2), y + (h/2)*f0) */
  temp = 0.5 * h;
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (temp*f0[i]);
  }

  rtsiSetT(si, t + temp);
  rtsiSetdX(si, f1);
  SYSTEM0_step0();
  SYSTEM0_derivatives();

  /* f2 = f(t + (h/2), y + (h/2)*f1) */
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (temp*f1[i]);
  }

  rtsiSetdX(si, f2);
  SYSTEM0_step0();
  SYSTEM0_derivatives();

  /* f3 = f(t + h, y + h*f2) */
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (h*f2[i]);
  }

  rtsiSetT(si, tnew);
  rtsiSetdX(si, f3);
  SYSTEM0_step0();
  SYSTEM0_derivatives();

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

/* Model step function for TID0 */
void SYSTEM0_step0(void)               /* Sample time: [0.0s, 0.0s] */
{
  /* local block i/o variables */
  real_T rtb_q0q1q2q3[4];
  real32_T rtb_DataTypeConversion;
  real_T rtb_IntegratorSecondOrder_o2_c1;
  real_T rtb_IntegratorSecondOrder_o1_c;
  real_T rtb_IntegratorSecondOrder_o2_c;
  real_T rtb_IntegratorSecondOrder_o1_l;
  real_T rtb_VectorConcatenate[9];
  real_T rtb_forces[3];
  real_T rtb_Product_k4[3];
  real_T rtb_Saturation_p[3];
  real_T rtb_Product2_pl;
  real32_T rtb_DataTypeConversion1;
  real32_T rtb_DataTypeConversion10;
  real32_T rtb_DataTypeConversion11;
  real32_T rtb_DataTypeConversion12;
  real32_T rtb_DataTypeConversion13;
  real32_T rtb_DataTypeConversion14;
  real32_T rtb_DataTypeConversion15;
  real32_T rtb_DataTypeConversion16;
  real32_T rtb_DataTypeConversion18;
  real32_T rtb_DataTypeConversion19;
  real32_T rtb_DataTypeConversion2;
  real32_T rtb_DataTypeConversion3;
  real32_T rtb_DataTypeConversion4;
  real32_T rtb_DataTypeConversion5;
  real32_T rtb_DataTypeConversion6;
  real_T tmp[3];
  int16_T i;
  real_T rtb_RPM2RADS_idx_0;
  real_T rtb_RPM2RADS_idx_1;
  real_T rtb_RPM2RADS_idx_2;
  real_T rtb_RPM2RADS_idx_3;
  real_T rtb_thrust_idx_0;
  real32_T rtb_DataTypeConversion17_idx_0;
  real32_T rtb_DataTypeConversion17_idx_1;
  real32_T rtb_DataTypeConversion17_idx_2;
  real_T rtb_thrust_idx_1;
  real_T rtb_thrust_idx_2;
  real_T rtb_Sum_h_idx_0;
  real_T rtb_Sum_h_idx_1;
  real_T rtb_Sum_h_idx_2;
  real_T rtb_ZeroOrderHold2_idx_0;
  real_T rtb_ZeroOrderHold2_idx_1;
  real_T rtb_ZeroOrderHold2_idx_2;
  if (rtmIsMajorTimeStep(SYSTEM0_M)) {
    /* set solver stop time */
    if (!(SYSTEM0_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&SYSTEM0_M->solverInfo,
                            ((SYSTEM0_M->Timing.clockTickH0 + 1) *
        SYSTEM0_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&SYSTEM0_M->solverInfo,
                            ((SYSTEM0_M->Timing.clockTick0 + 1) *
        SYSTEM0_M->Timing.stepSize0 + SYSTEM0_M->Timing.clockTickH0 *
        SYSTEM0_M->Timing.stepSize0 * 4294967296.0));
    }

    /* Update the flag to indicate when data transfers from
     *  Sample time: [0.0s, 0.0s] to Sample time: [0.04s, 0.0s]  */
    (SYSTEM0_M->Timing.RateInteraction.TID0_2)++;
    if ((SYSTEM0_M->Timing.RateInteraction.TID0_2) > 9) {
      SYSTEM0_M->Timing.RateInteraction.TID0_2 = 0;
    }

    /* Update the flag to indicate when data transfers from
     *  Sample time: [0.0s, 0.0s] to Sample time: [1.0s, 0.0s]  */
    (SYSTEM0_M->Timing.RateInteraction.TID0_3)++;
    if ((SYSTEM0_M->Timing.RateInteraction.TID0_3) > 249) {
      SYSTEM0_M->Timing.RateInteraction.TID0_3 = 0;
    }

    /* Update the flag to indicate when data transfers from
     *  Sample time: [0.004s, 0.0s] to Sample time: [0.04s, 0.0s]  */
    (SYSTEM0_M->Timing.RateInteraction.TID1_2)++;
    if ((SYSTEM0_M->Timing.RateInteraction.TID1_2) > 9) {
      SYSTEM0_M->Timing.RateInteraction.TID1_2 = 0;
    }

    /* Update the flag to indicate when data transfers from
     *  Sample time: [0.004s, 0.0s] to Sample time: [1.0s, 0.0s]  */
    (SYSTEM0_M->Timing.RateInteraction.TID1_3)++;
    if ((SYSTEM0_M->Timing.RateInteraction.TID1_3) > 249) {
      SYSTEM0_M->Timing.RateInteraction.TID1_3 = 0;
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(SYSTEM0_M)) {
    SYSTEM0_M->Timing.t[0] = rtsiGetT(&SYSTEM0_M->solverInfo);
  }

  /* Clock: '<S3>/runtime' */
  SYSTEM0_B.runtime = SYSTEM0_M->Timing.t[0];
  if (rtmIsMajorTimeStep(SYSTEM0_M)) {
    /* DataTypeConversion: '<S2>/Data Type Conversion' incorporates:
     *  Fcn: '<S3>/Fcn'
     */
    rtb_DataTypeConversion = (real32_T)(SYSTEM0_B.runtime * 1000000.0);
  }

  /* Integrator: '<S10>/q0 q1 q2 q3' */
  if (SYSTEM0_DW.q0q1q2q3_IWORK.IcNeedsLoading) {
    SYSTEM0_X.q0q1q2q3_CSTATE[0] = SYSTEM0_ConstB.q0;
    SYSTEM0_X.q0q1q2q3_CSTATE[1] = SYSTEM0_ConstB.q1;
    SYSTEM0_X.q0q1q2q3_CSTATE[2] = SYSTEM0_ConstB.q2;
    SYSTEM0_X.q0q1q2q3_CSTATE[3] = SYSTEM0_ConstB.q3;
  }

  rtb_q0q1q2q3[0] = SYSTEM0_X.q0q1q2q3_CSTATE[0];
  rtb_q0q1q2q3[1] = SYSTEM0_X.q0q1q2q3_CSTATE[1];
  rtb_q0q1q2q3[2] = SYSTEM0_X.q0q1q2q3_CSTATE[2];
  rtb_q0q1q2q3[3] = SYSTEM0_X.q0q1q2q3_CSTATE[3];

  /* Sqrt: '<S32>/sqrt' incorporates:
   *  Product: '<S33>/Product'
   *  Product: '<S33>/Product1'
   *  Product: '<S33>/Product2'
   *  Product: '<S33>/Product3'
   *  Sum: '<S33>/Sum'
   */
  rtb_IntegratorSecondOrder_o2_c1 = sqrt(((rtb_q0q1q2q3[0] * rtb_q0q1q2q3[0] +
    rtb_q0q1q2q3[1] * rtb_q0q1q2q3[1]) + rtb_q0q1q2q3[2] * rtb_q0q1q2q3[2]) +
    rtb_q0q1q2q3[3] * rtb_q0q1q2q3[3]);

  /* Product: '<S31>/Product' */
  rtb_IntegratorSecondOrder_o1_c = rtb_q0q1q2q3[0] /
    rtb_IntegratorSecondOrder_o2_c1;

  /* Product: '<S31>/Product1' */
  rtb_IntegratorSecondOrder_o2_c = rtb_q0q1q2q3[1] /
    rtb_IntegratorSecondOrder_o2_c1;

  /* Product: '<S31>/Product2' */
  rtb_IntegratorSecondOrder_o1_l = rtb_q0q1q2q3[2] /
    rtb_IntegratorSecondOrder_o2_c1;

  /* Product: '<S31>/Product3' */
  rtb_IntegratorSecondOrder_o2_c1 = rtb_q0q1q2q3[3] /
    rtb_IntegratorSecondOrder_o2_c1;

  /* Sum: '<S21>/Sum' incorporates:
   *  Product: '<S21>/Product'
   *  Product: '<S21>/Product1'
   *  Product: '<S21>/Product2'
   *  Product: '<S21>/Product3'
   */
  rtb_VectorConcatenate[0] = ((rtb_IntegratorSecondOrder_o1_c *
    rtb_IntegratorSecondOrder_o1_c + rtb_IntegratorSecondOrder_o2_c *
    rtb_IntegratorSecondOrder_o2_c) - rtb_IntegratorSecondOrder_o1_l *
    rtb_IntegratorSecondOrder_o1_l) - rtb_IntegratorSecondOrder_o2_c1 *
    rtb_IntegratorSecondOrder_o2_c1;

  /* Gain: '<S24>/Gain' incorporates:
   *  Product: '<S24>/Product2'
   *  Product: '<S24>/Product3'
   *  Sum: '<S24>/Sum'
   */
  rtb_VectorConcatenate[1] = (rtb_IntegratorSecondOrder_o2_c *
    rtb_IntegratorSecondOrder_o1_l - rtb_IntegratorSecondOrder_o2_c1 *
    rtb_IntegratorSecondOrder_o1_c) * 2.0;

  /* Gain: '<S27>/Gain' incorporates:
   *  Product: '<S27>/Product1'
   *  Product: '<S27>/Product2'
   *  Sum: '<S27>/Sum'
   */
  rtb_VectorConcatenate[2] = (rtb_IntegratorSecondOrder_o1_c *
    rtb_IntegratorSecondOrder_o1_l + rtb_IntegratorSecondOrder_o2_c *
    rtb_IntegratorSecondOrder_o2_c1) * 2.0;

  /* Gain: '<S22>/Gain' incorporates:
   *  Product: '<S22>/Product2'
   *  Product: '<S22>/Product3'
   *  Sum: '<S22>/Sum'
   */
  rtb_VectorConcatenate[3] = (rtb_IntegratorSecondOrder_o2_c1 *
    rtb_IntegratorSecondOrder_o1_c + rtb_IntegratorSecondOrder_o2_c *
    rtb_IntegratorSecondOrder_o1_l) * 2.0;

  /* Sum: '<S25>/Sum' incorporates:
   *  Product: '<S25>/Product'
   *  Product: '<S25>/Product1'
   *  Product: '<S25>/Product2'
   *  Product: '<S25>/Product3'
   */
  rtb_VectorConcatenate[4] = ((rtb_IntegratorSecondOrder_o1_c *
    rtb_IntegratorSecondOrder_o1_c - rtb_IntegratorSecondOrder_o2_c *
    rtb_IntegratorSecondOrder_o2_c) + rtb_IntegratorSecondOrder_o1_l *
    rtb_IntegratorSecondOrder_o1_l) - rtb_IntegratorSecondOrder_o2_c1 *
    rtb_IntegratorSecondOrder_o2_c1;

  /* Gain: '<S28>/Gain' incorporates:
   *  Product: '<S28>/Product1'
   *  Product: '<S28>/Product2'
   *  Sum: '<S28>/Sum'
   */
  rtb_VectorConcatenate[5] = (rtb_IntegratorSecondOrder_o1_l *
    rtb_IntegratorSecondOrder_o2_c1 - rtb_IntegratorSecondOrder_o1_c *
    rtb_IntegratorSecondOrder_o2_c) * 2.0;

  /* Gain: '<S23>/Gain' incorporates:
   *  Product: '<S23>/Product1'
   *  Product: '<S23>/Product2'
   *  Sum: '<S23>/Sum'
   */
  rtb_VectorConcatenate[6] = (rtb_IntegratorSecondOrder_o2_c *
    rtb_IntegratorSecondOrder_o2_c1 - rtb_IntegratorSecondOrder_o1_c *
    rtb_IntegratorSecondOrder_o1_l) * 2.0;

  /* Gain: '<S26>/Gain' incorporates:
   *  Product: '<S26>/Product1'
   *  Product: '<S26>/Product2'
   *  Sum: '<S26>/Sum'
   */
  rtb_VectorConcatenate[7] = (rtb_IntegratorSecondOrder_o1_c *
    rtb_IntegratorSecondOrder_o2_c + rtb_IntegratorSecondOrder_o1_l *
    rtb_IntegratorSecondOrder_o2_c1) * 2.0;

  /* Sum: '<S29>/Sum' incorporates:
   *  Product: '<S29>/Product'
   *  Product: '<S29>/Product1'
   *  Product: '<S29>/Product2'
   *  Product: '<S29>/Product3'
   */
  rtb_VectorConcatenate[8] = ((rtb_IntegratorSecondOrder_o1_c *
    rtb_IntegratorSecondOrder_o1_c - rtb_IntegratorSecondOrder_o2_c *
    rtb_IntegratorSecondOrder_o2_c) - rtb_IntegratorSecondOrder_o1_l *
    rtb_IntegratorSecondOrder_o1_l) + rtb_IntegratorSecondOrder_o2_c1 *
    rtb_IntegratorSecondOrder_o2_c1;

  /* Sum: '<S56>/Sum1' incorporates:
   *  UnaryMinus: '<S56>/Ze2height'
   */
  SYSTEM0_B.Sum1 = -SYSTEM0_X.xeyeze_CSTATE[2];
  if (rtmIsMajorTimeStep(SYSTEM0_M)) {
    /* Sum: '<S4>/Add8' incorporates:
     *  Gain: '<S4>/Gain3'
     */
    SYSTEM0_B.Add8[0] = 0.0;
    SYSTEM0_B.Add8[1] = 0.0;
    SYSTEM0_B.Add8[2] = 5.0 * SYSTEM0_B.Sum1 + -11.772;

    /* Gain: '<S4>/Gain4' */
    SYSTEM0_B.Gain4 = -SYSTEM0_B.Sum1;
  }

  /* Switch: '<S4>/Switch' incorporates:
   *  Product: '<S4>/Matrix Multiply'
   */
  if (SYSTEM0_B.Gain4 >= 0.0) {
    for (i = 0; i < 3; i++) {
      SYSTEM0_B.Product[i] = 0.0;
      SYSTEM0_B.Product[i] += rtb_VectorConcatenate[i] * SYSTEM0_B.Add8[0];
      SYSTEM0_B.Product[i] += rtb_VectorConcatenate[i + 3] * SYSTEM0_B.Add8[1];
      SYSTEM0_B.Product[i] += rtb_VectorConcatenate[i + 6] * SYSTEM0_B.Add8[2];
    }
  } else {
    SYSTEM0_B.Product[0] = 0.0;
    SYSTEM0_B.Product[1] = 0.0;
    SYSTEM0_B.Product[2] = 0.0;
  }

  /* End of Switch: '<S4>/Switch' */
  if (rtmIsMajorTimeStep(SYSTEM0_M)) {
    /* Memory: '<S4>/Memory2' */
    SYSTEM0_B.Memory2 = SYSTEM0_DW.Memory2_PreviousInput;
  }

  /* Gain: '<S9>/RPM2RADS' incorporates:
   *  Gain: '<S9>/V2RPM'
   *  SecondOrderIntegrator: '<S49>/Integrator, Second-Order'
   *  SecondOrderIntegrator: '<S50>/Integrator, Second-Order'
   *  SecondOrderIntegrator: '<S51>/Integrator, Second-Order'
   *  SecondOrderIntegrator: '<S52>/Integrator, Second-Order'
   */
  rtb_RPM2RADS_idx_0 = 950.0 * SYSTEM0_X.IntegratorSecondOrder_CSTATE[0] *
    0.10471975511965977;
  rtb_RPM2RADS_idx_1 = 950.0 * SYSTEM0_X.IntegratorSecondOrder_CSTATE_h[0] *
    0.10471975511965977;
  rtb_RPM2RADS_idx_2 = 950.0 * SYSTEM0_X.IntegratorSecondOrder_CSTATE_n[0] *
    0.10471975511965977;
  rtb_RPM2RADS_idx_3 = 950.0 * SYSTEM0_X.IntegratorSecondOrder_CSTATE_d[0] *
    0.10471975511965977;

  /* MATLAB Function: '<S8>/multicopter' incorporates:
   *  Constant: '<S8>/h_ref10'
   *  Constant: '<S8>/h_ref12'
   *  Constant: '<S8>/h_ref8'
   */
  /* MATLAB Function 'SYSTEM/dynamic_model/Dynamics/Subsystem/multicopter': '<S46>:1' */
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
  /* '<S46>:1:36' */
  /*  drag force */
  /* --------------------------------Thrust Model------------------------------ */
  /* '<S46>:1:39' */
  rtb_IntegratorSecondOrder_o2_c1 = 1.2247084269789534E-5 * SYSTEM0_B.Memory2;
  rtb_thrust_idx_0 = rtb_RPM2RADS_idx_0 * rtb_RPM2RADS_idx_0 *
    rtb_IntegratorSecondOrder_o2_c1;
  rtb_thrust_idx_1 = rtb_RPM2RADS_idx_1 * rtb_RPM2RADS_idx_1 *
    rtb_IntegratorSecondOrder_o2_c1;
  rtb_thrust_idx_2 = rtb_RPM2RADS_idx_2 * rtb_RPM2RADS_idx_2 *
    rtb_IntegratorSecondOrder_o2_c1;
  rtb_IntegratorSecondOrder_o2_c1 *= rtb_RPM2RADS_idx_3 * rtb_RPM2RADS_idx_3;

  /*  rotor thrust */
  /* -------------------------------------------------------------------------- */
  /* '<S46>:1:42' */
  /* '<S46>:1:43' */
  rtb_forces[2] = (-SYSTEM0_X.ubvbwb_CSTATE[2] * 10.0 * SYSTEM0_B.Memory2 *
                   0.18845573684677208 + rtb_VectorConcatenate[8] * 9.81 * 1.2)
    - (((rtb_thrust_idx_0 + rtb_thrust_idx_1) + rtb_thrust_idx_2) +
       rtb_IntegratorSecondOrder_o2_c1);

  /* ==================================Moments================================= */
  /*  Thrusts contributions to momentum */
  /* '<S46>:1:48' */
  rtb_IntegratorSecondOrder_o1_c = rtb_thrust_idx_0 * 0.2 * 1.4142135623730951 /
    2.0;
  rtb_IntegratorSecondOrder_o2_c = -rtb_thrust_idx_1 * 0.2 * 1.4142135623730951 /
    2.0;
  rtb_IntegratorSecondOrder_o1_l = -rtb_thrust_idx_2 * 0.2 * 1.4142135623730951 /
    2.0;

  /*  x moment */
  /* '<S46>:1:49' */
  rtb_thrust_idx_0 = rtb_thrust_idx_0 * 0.2 * 1.4142135623730951 / 2.0;
  rtb_thrust_idx_1 = rtb_thrust_idx_1 * 0.2 * 1.4142135623730951 / 2.0;
  rtb_thrust_idx_2 = -rtb_thrust_idx_2 * 0.2 * 1.4142135623730951 / 2.0;

  /* Product: '<S6>/Product' incorporates:
   *  Constant: '<S12>/Constant'
   *  Constant: '<S8>/h_ref12'
   *  MATLAB Function: '<S8>/multicopter'
   *  Sum: '<S4>/Add'
   *  Sum: '<S4>/Add7'
   */
  /*  y moment */
  /*  Torques contributions */
  /* momentum_x = sum(abs(MIX(:, 3)) .* rotor_inertia .* rotors) * omega(1);     % x rotor momentum */
  /* momentum_y = sum(abs(MIX(:, 3)) .* rotor_inertia .* rotors) * omega(2);     % y rotor momentum */
  /* momentum_z = sum(MIX(:, 3) .* rotor_inertia .* rotors) * omega(3);          % z rotor momentum */
  /* --------------------------------Torque Model------------------------------ */
  /* '<S46>:1:57' */
  /*  rotor torque */
  /* -------------------------------------------------------------------------- */
  /* '<S46>:1:60' */
  /*  - [momentum_x; momentum_y; momentum_z]; */
  /* ========================================================================== */
  SYSTEM0_B.Product[0] = (SYSTEM0_B.Product[0] + (-SYSTEM0_X.ubvbwb_CSTATE[0] *
    10.0 * SYSTEM0_B.Memory2 * 0.016813708498984763 + rtb_VectorConcatenate[6] *
    9.81 * 1.2)) / 1.2;
  SYSTEM0_B.Product[1] = (SYSTEM0_B.Product[1] + (-SYSTEM0_X.ubvbwb_CSTATE[1] *
    10.0 * SYSTEM0_B.Memory2 * 0.018813708498984762 + rtb_VectorConcatenate[7] *
    9.81 * 1.2)) / 1.2;
  SYSTEM0_B.Product[2] = (SYSTEM0_B.Product[2] + rtb_forces[2]) / 1.2;
  if (rtmIsMajorTimeStep(SYSTEM0_M)) {
    /* ZeroOrderHold: '<S129>/Zero-Order Hold1' */
    rtb_Sum_h_idx_0 = SYSTEM0_B.Product[0];
    rtb_Sum_h_idx_1 = SYSTEM0_B.Product[1];
    rtb_Sum_h_idx_2 = SYSTEM0_B.Product[2];
  }

  /* Product: '<S5>/Matrix Multiply1' */
  for (i = 0; i < 3; i++) {
    SYSTEM0_B.MatrixMultiply1[i] = 0.0;
    SYSTEM0_B.MatrixMultiply1[i] += rtb_VectorConcatenate[i + 6] * 9.81;
  }

  /* End of Product: '<S5>/Matrix Multiply1' */
  if (rtmIsMajorTimeStep(SYSTEM0_M)) {
    /* ZeroOrderHold: '<S129>/Zero-Order Hold2' */
    rtb_ZeroOrderHold2_idx_0 = SYSTEM0_B.MatrixMultiply1[0];
    rtb_ZeroOrderHold2_idx_1 = SYSTEM0_B.MatrixMultiply1[1];
    rtb_ZeroOrderHold2_idx_2 = SYSTEM0_B.MatrixMultiply1[2];
  }

  /* Integrator: '<S6>/p,q,r ' */
  SYSTEM0_B.pqr[0] = SYSTEM0_X.pqr_CSTATE[0];
  SYSTEM0_B.pqr[1] = SYSTEM0_X.pqr_CSTATE[1];
  SYSTEM0_B.pqr[2] = SYSTEM0_X.pqr_CSTATE[2];
  if (rtmIsMajorTimeStep(SYSTEM0_M)) {
    /* Gain: '<S129>/Gain' incorporates:
     *  Constant: '<S129>/wl_ins'
     *  Constant: '<S5>/center of gravity'
     *  Sum: '<S129>/Sum7'
     */
    rtb_Product_k4[0] = 0.0;
    rtb_Product_k4[1] = 0.0;
    rtb_Product_k4[2] = 0.0;

    /* Sum: '<S140>/Sum' incorporates:
     *  Product: '<S142>/i x j'
     *  Product: '<S142>/j x k'
     *  Product: '<S142>/k x i'
     *  Product: '<S143>/i x k'
     *  Product: '<S143>/j x i'
     *  Product: '<S143>/k x j'
     */
    rtb_Saturation_p[0] = 0.0;
    rtb_Saturation_p[1] = 0.0;
    rtb_Saturation_p[2] = 0.0;
  }

  /* Product: '<S37>/Product' */
  for (i = 0; i < 3; i++) {
    rtb_forces[i] = SYSTEM0_ConstB.Selector[i + 6] * SYSTEM0_B.pqr[2] +
      (SYSTEM0_ConstB.Selector[i + 3] * SYSTEM0_B.pqr[1] +
       SYSTEM0_ConstB.Selector[i] * SYSTEM0_B.pqr[0]);
  }

  /* End of Product: '<S37>/Product' */

  /* Product: '<S41>/j x i' */
  rtb_Product2_pl = SYSTEM0_B.pqr[1] * rtb_forces[0];

  /* Product: '<S11>/Product2' incorporates:
   *  Constant: '<S8>/h_ref10'
   *  MATLAB Function: '<S8>/multicopter'
   *  Product: '<S40>/i x j'
   *  Product: '<S40>/j x k'
   *  Product: '<S40>/k x i'
   *  Product: '<S41>/i x k'
   *  Product: '<S41>/k x j'
   *  Sum: '<S11>/Sum2'
   *  Sum: '<S39>/Sum'
   */
  tmp[0] = (rtb_IntegratorSecondOrder_o2_c1 * 0.2 * 1.4142135623730951 / 2.0 +
            ((rtb_IntegratorSecondOrder_o1_c + rtb_IntegratorSecondOrder_o2_c) +
             rtb_IntegratorSecondOrder_o1_l)) - (SYSTEM0_B.pqr[1] * rtb_forces[2]
    - SYSTEM0_B.pqr[2] * rtb_forces[1]);
  tmp[1] = (-rtb_IntegratorSecondOrder_o2_c1 * 0.2 * 1.4142135623730951 / 2.0 +
            ((rtb_thrust_idx_0 + rtb_thrust_idx_1) + rtb_thrust_idx_2)) -
    (SYSTEM0_B.pqr[2] * rtb_forces[0] - SYSTEM0_B.pqr[0] * rtb_forces[2]);
  tmp[2] = (((-7.129366502583864E-8 * SYSTEM0_B.Memory2 * (rtb_RPM2RADS_idx_0 *
    rtb_RPM2RADS_idx_0) + 7.129366502583864E-8 * SYSTEM0_B.Memory2 *
              (rtb_RPM2RADS_idx_1 * rtb_RPM2RADS_idx_1)) + -7.129366502583864E-8
             * SYSTEM0_B.Memory2 * (rtb_RPM2RADS_idx_2 * rtb_RPM2RADS_idx_2)) +
            7.129366502583864E-8 * SYSTEM0_B.Memory2 * (rtb_RPM2RADS_idx_3 *
             rtb_RPM2RADS_idx_3)) - (SYSTEM0_B.pqr[0] * rtb_forces[1] -
    rtb_Product2_pl);
  rt_mrdivide_U1d1x3_U2d3x3_Yd1x3(tmp, SYSTEM0_ConstB.Selector2,
    SYSTEM0_B.Product2);
  if (rtmIsMajorTimeStep(SYSTEM0_M)) {
    /* Sum: '<S134>/Sum' */
    rtb_Product_k4[0] = 0.0;
    rtb_Product_k4[1] = 0.0;
    rtb_Product_k4[2] = 0.0;

    /* Sum: '<S129>/Sum' */
    rtb_Sum_h_idx_0 = (rtb_Sum_h_idx_0 - rtb_ZeroOrderHold2_idx_0) +
      rtb_Saturation_p[0];
    rtb_Sum_h_idx_1 = (rtb_Sum_h_idx_1 - rtb_ZeroOrderHold2_idx_1) +
      rtb_Saturation_p[1];
    rtb_IntegratorSecondOrder_o2_c1 = (rtb_Sum_h_idx_2 -
      rtb_ZeroOrderHold2_idx_2) + rtb_Saturation_p[2];

    /* Product: '<S129>/Product' incorporates:
     *  Constant: '<S129>/Scale factors & Cross-coupling  errors'
     */
    for (i = 0; i < 3; i++) {
      rtb_Saturation_p[i] = SYSTEM0_ConstP.pooled10[i + 6] *
        rtb_IntegratorSecondOrder_o2_c1 + (SYSTEM0_ConstP.pooled10[i + 3] *
        rtb_Sum_h_idx_1 + SYSTEM0_ConstP.pooled10[i] * rtb_Sum_h_idx_0);
    }

    /* End of Product: '<S129>/Product' */

    /* Sum: '<S129>/Sum4' */
    rtb_forces[2] = rtb_Saturation_p[2];

    /* Saturate: '<S129>/Saturation' incorporates:
     *  Gain: '<S132>/Output'
     *  RandomNumber: '<S132>/White Noise'
     *  Sum: '<S129>/Sum1'
     *  Sum: '<S129>/Sum4'
     */
    rtb_IntegratorSecondOrder_o1_c = 0.011180339887498949 *
      SYSTEM0_DW.NextOutput[0] + rtb_Saturation_p[0];
    rtb_IntegratorSecondOrder_o2_c = 0.011180339887498949 *
      SYSTEM0_DW.NextOutput[1] + rtb_Saturation_p[1];
    if (rtb_IntegratorSecondOrder_o2_c > 19.62) {
      rtb_Saturation_p[1] = 19.62;
    } else if (rtb_IntegratorSecondOrder_o2_c < -19.62) {
      rtb_Saturation_p[1] = -19.62;
    } else {
      rtb_Saturation_p[1] = rtb_IntegratorSecondOrder_o2_c;
    }

    rtb_IntegratorSecondOrder_o2_c = 0.011180339887498949 *
      SYSTEM0_DW.NextOutput[2] + rtb_forces[2];
    if (rtb_IntegratorSecondOrder_o2_c > 19.62) {
      rtb_Saturation_p[2] = 19.62;
    } else if (rtb_IntegratorSecondOrder_o2_c < -19.62) {
      rtb_Saturation_p[2] = -19.62;
    } else {
      rtb_Saturation_p[2] = rtb_IntegratorSecondOrder_o2_c;
    }

    if (rtb_IntegratorSecondOrder_o1_c > 19.62) {
      /* DataTypeConversion: '<S2>/Data Type Conversion1' */
      rtb_DataTypeConversion1 = 19.62F;
    } else if (rtb_IntegratorSecondOrder_o1_c < -19.62) {
      /* DataTypeConversion: '<S2>/Data Type Conversion1' */
      rtb_DataTypeConversion1 = -19.62F;
    } else {
      /* DataTypeConversion: '<S2>/Data Type Conversion1' */
      rtb_DataTypeConversion1 = (real32_T)rtb_IntegratorSecondOrder_o1_c;
    }

    /* End of Saturate: '<S129>/Saturation' */

    /* RateTransition: '<S3>/Rate Transition1' */
    if (SYSTEM0_M->Timing.RateInteraction.TID1_2 == 1) {
      SYSTEM0_B.RateTransition1 = SYSTEM0_DW.RateTransition1_Buffer0;
    }

    /* End of RateTransition: '<S3>/Rate Transition1' */

    /* DataTypeConversion: '<S2>/Data Type Conversion10' */
    rtb_DataTypeConversion10 = (real32_T)SYSTEM0_B.RateTransition1;

    /* RateTransition: '<S3>/Rate Transition3' */
    if (SYSTEM0_M->Timing.RateInteraction.TID1_2 == 1) {
      SYSTEM0_B.RateTransition3 = SYSTEM0_DW.RateTransition3_Buffer0;
    }

    /* End of RateTransition: '<S3>/Rate Transition3' */

    /* DataTypeConversion: '<S2>/Data Type Conversion11' */
    rtb_DataTypeConversion11 = (real32_T)SYSTEM0_B.RateTransition3;

    /* RandomNumber: '<S53>/Random Number' */
    SYSTEM0_B.RandomNumber = SYSTEM0_DW.NextOutput_n;
  }

  /* Sum: '<S53>/Add1' */
  rtb_Product2_pl = SYSTEM0_B.Sum1 + SYSTEM0_B.RandomNumber;

  /* Saturate: '<S53>/Saturation1' */
  if (rtb_Product2_pl >= 0.0) {
    SYSTEM0_B.Saturation1 = rtb_Product2_pl;
  } else {
    SYSTEM0_B.Saturation1 = 0.0;
  }

  /* End of Saturate: '<S53>/Saturation1' */
  if (rtmIsMajorTimeStep(SYSTEM0_M)) {
    /* DataTypeConversion: '<S2>/Data Type Conversion12' */
    rtb_DataTypeConversion12 = (real32_T)SYSTEM0_B.Saturation1;
  }

  /* Saturate: '<S58>/Limit  altitude  to troposhere' */
  if (SYSTEM0_B.Sum1 > 11000.0) {
    rtb_IntegratorSecondOrder_o2_c1 = 11000.0;
  } else if (SYSTEM0_B.Sum1 < 0.0) {
    rtb_IntegratorSecondOrder_o2_c1 = 0.0;
  } else {
    rtb_IntegratorSecondOrder_o2_c1 = SYSTEM0_B.Sum1;
  }

  /* Sum: '<S58>/Sum1' incorporates:
   *  Constant: '<S58>/Sea Level  Temperature'
   *  Gain: '<S58>/Lapse Rate'
   *  Saturate: '<S58>/Limit  altitude  to troposhere'
   */
  SYSTEM0_B.Sum1_e = 288.15 - 0.0065 * rtb_IntegratorSecondOrder_o2_c1;
  if (rtmIsMajorTimeStep(SYSTEM0_M)) {
    /* Sum: '<S65>/Sum2' incorporates:
     *  Constant: '<S65>/K2C'
     *  RandomNumber: '<S65>/Random Number'
     *  Sum: '<S65>/Add1'
     */
    rtb_IntegratorSecondOrder_o1_c = (273.15 + SYSTEM0_B.Sum1_e) +
      SYSTEM0_DW.NextOutput_h;

    /* Saturate: '<S65>/Saturation' */
    if (rtb_IntegratorSecondOrder_o1_c > 85.0) {
      /* DataTypeConversion: '<S2>/Data Type Conversion13' */
      rtb_DataTypeConversion13 = 85.0F;
    } else if (rtb_IntegratorSecondOrder_o1_c < 40.0) {
      /* DataTypeConversion: '<S2>/Data Type Conversion13' */
      rtb_DataTypeConversion13 = 40.0F;
    } else {
      /* DataTypeConversion: '<S2>/Data Type Conversion13' */
      rtb_DataTypeConversion13 = (real32_T)rtb_IntegratorSecondOrder_o1_c;
    }

    /* End of Saturate: '<S65>/Saturation' */

    /* RateTransition: '<S3>/Rate Transition6' */
    if (SYSTEM0_M->Timing.RateInteraction.TID1_3 == 1) {
      SYSTEM0_B.RateTransition6 = SYSTEM0_DW.RateTransition6_Buffer0;
    }

    /* End of RateTransition: '<S3>/Rate Transition6' */

    /* DataTypeConversion: '<S2>/Data Type Conversion14' */
    rtb_DataTypeConversion14 = (real32_T)SYSTEM0_B.RateTransition6;

    /* RateTransition: '<S3>/Rate Transition5' */
    if (SYSTEM0_M->Timing.RateInteraction.TID1_3 == 1) {
      SYSTEM0_B.RateTransition5 = SYSTEM0_DW.RateTransition5_Buffer0;
    }

    /* End of RateTransition: '<S3>/Rate Transition5' */

    /* DataTypeConversion: '<S2>/Data Type Conversion15' */
    rtb_DataTypeConversion15 = (real32_T)SYSTEM0_B.RateTransition5;

    /* RateTransition: '<S3>/Rate Transition7' */
    if (SYSTEM0_M->Timing.RateInteraction.TID1_3 == 1) {
      SYSTEM0_B.RateTransition7 = SYSTEM0_DW.RateTransition7_Buffer0;
    }

    /* End of RateTransition: '<S3>/Rate Transition7' */

    /* DataTypeConversion: '<S2>/Data Type Conversion16' */
    rtb_DataTypeConversion16 = (real32_T)SYSTEM0_B.RateTransition7;

    /* RateTransition: '<S3>/Rate Transition8' */
    if (SYSTEM0_M->Timing.RateInteraction.TID1_3 == 1) {
      SYSTEM0_B.RateTransition8[0] = SYSTEM0_DW.RateTransition8_Buffer0[0];
      SYSTEM0_B.RateTransition8[1] = SYSTEM0_DW.RateTransition8_Buffer0[1];
      SYSTEM0_B.RateTransition8[2] = SYSTEM0_DW.RateTransition8_Buffer0[2];
    }

    /* End of RateTransition: '<S3>/Rate Transition8' */

    /* DataTypeConversion: '<S2>/Data Type Conversion17' */
    rtb_DataTypeConversion17_idx_0 = (real32_T)SYSTEM0_B.RateTransition8[0];
    rtb_DataTypeConversion17_idx_1 = (real32_T)SYSTEM0_B.RateTransition8[1];
    rtb_DataTypeConversion17_idx_2 = (real32_T)SYSTEM0_B.RateTransition8[2];

    /* RateTransition: '<S3>/Rate Transition10' */
    if (SYSTEM0_M->Timing.RateInteraction.TID1_3 == 1) {
      SYSTEM0_B.RateTransition10 = SYSTEM0_DW.RateTransition10_Buffer0;
    }

    /* End of RateTransition: '<S3>/Rate Transition10' */

    /* DataTypeConversion: '<S2>/Data Type Conversion18' */
    rtb_DataTypeConversion18 = (real32_T)SYSTEM0_B.RateTransition10;

    /* DataTypeConversion: '<S2>/Data Type Conversion19' incorporates:
     *  DotProduct: '<S2>/Dot Product'
     *  Sqrt: '<S2>/Sqrt'
     */
    rtb_DataTypeConversion19 = (real32_T)sqrt((SYSTEM0_B.RateTransition8[0] *
      SYSTEM0_B.RateTransition8[0] + SYSTEM0_B.RateTransition8[1] *
      SYSTEM0_B.RateTransition8[1]) + SYSTEM0_B.RateTransition8[2] *
      SYSTEM0_B.RateTransition8[2]);

    /* DataTypeConversion: '<S2>/Data Type Conversion2' */
    rtb_DataTypeConversion2 = (real32_T)rtb_Saturation_p[1];

    /* DataTypeConversion: '<S2>/Data Type Conversion3' */
    rtb_DataTypeConversion3 = (real32_T)rtb_Saturation_p[2];

    /* Product: '<S130>/Product' incorporates:
     *  Constant: '<S130>/Scale factors & Cross-coupling  errors '
     *  ZeroOrderHold: '<S130>/Zero-Order Hold'
     */
    for (i = 0; i < 3; i++) {
      rtb_Product_k4[i] = SYSTEM0_ConstP.pooled10[i + 6] * SYSTEM0_B.pqr[2] +
        (SYSTEM0_ConstP.pooled10[i + 3] * SYSTEM0_B.pqr[1] +
         SYSTEM0_ConstP.pooled10[i] * SYSTEM0_B.pqr[0]);
    }

    /* End of Product: '<S130>/Product' */
  }

  /* Gain: '<S128>/Unit Conversion' */
  SYSTEM0_B.UnitConversion[0] = 0.10197162129779283 * SYSTEM0_B.Product[0];
  SYSTEM0_B.UnitConversion[1] = 0.10197162129779283 * SYSTEM0_B.Product[1];
  SYSTEM0_B.UnitConversion[2] = 0.10197162129779283 * SYSTEM0_B.Product[2];
  if (rtmIsMajorTimeStep(SYSTEM0_M)) {
    /* Saturate: '<S130>/Saturation' incorporates:
     *  Gain: '<S149>/Output'
     *  RandomNumber: '<S149>/White Noise'
     *  Sum: '<S130>/Sum1'
     *  Sum: '<S130>/Sum4'
     */
    rtb_IntegratorSecondOrder_o1_c = 0.00070710678118654762 *
      SYSTEM0_DW.NextOutput_p[0] + rtb_Product_k4[0];
    rtb_IntegratorSecondOrder_o2_c = 0.00070710678118654762 *
      SYSTEM0_DW.NextOutput_p[1] + rtb_Product_k4[1];
    rtb_IntegratorSecondOrder_o2_c1 = 0.00070710678118654762 *
      SYSTEM0_DW.NextOutput_p[2] + rtb_Product_k4[2];
    if (rtb_IntegratorSecondOrder_o1_c > 4.36) {
      /* DataTypeConversion: '<S2>/Data Type Conversion4' */
      rtb_DataTypeConversion4 = 4.36F;
    } else if (rtb_IntegratorSecondOrder_o1_c < -4.36) {
      /* DataTypeConversion: '<S2>/Data Type Conversion4' */
      rtb_DataTypeConversion4 = -4.36F;
    } else {
      /* DataTypeConversion: '<S2>/Data Type Conversion4' */
      rtb_DataTypeConversion4 = (real32_T)rtb_IntegratorSecondOrder_o1_c;
    }

    if (rtb_IntegratorSecondOrder_o2_c > 4.36) {
      /* DataTypeConversion: '<S2>/Data Type Conversion5' */
      rtb_DataTypeConversion5 = 4.36F;
    } else if (rtb_IntegratorSecondOrder_o2_c < -4.36) {
      /* DataTypeConversion: '<S2>/Data Type Conversion5' */
      rtb_DataTypeConversion5 = -4.36F;
    } else {
      /* DataTypeConversion: '<S2>/Data Type Conversion5' */
      rtb_DataTypeConversion5 = (real32_T)rtb_IntegratorSecondOrder_o2_c;
    }

    if (rtb_IntegratorSecondOrder_o2_c1 > 4.36) {
      /* DataTypeConversion: '<S2>/Data Type Conversion6' */
      rtb_DataTypeConversion6 = 4.36F;
    } else if (rtb_IntegratorSecondOrder_o2_c1 < -4.36) {
      /* DataTypeConversion: '<S2>/Data Type Conversion6' */
      rtb_DataTypeConversion6 = -4.36F;
    } else {
      /* DataTypeConversion: '<S2>/Data Type Conversion6' */
      rtb_DataTypeConversion6 = (real32_T)rtb_IntegratorSecondOrder_o2_c1;
    }

    /* End of Saturate: '<S130>/Saturation' */

    /* RandomNumber: '<S60>/Random Number' */
    SYSTEM0_B.RandomNumber_n = SYSTEM0_DW.NextOutput_a;
  }

  /* Sum: '<S60>/Sum2' incorporates:
   *  Product: '<S60>/Matrix Multiply2'
   */
  for (i = 0; i < 3; i++) {
    SYSTEM0_B.Sum2[i] = (rtb_VectorConcatenate[i + 3] * 0.49999999999999994 +
                         rtb_VectorConcatenate[i] * 0.86602540378443871) +
      SYSTEM0_B.RandomNumber_n;
  }

  /* End of Sum: '<S60>/Sum2' */
  if (rtmIsMajorTimeStep(SYSTEM0_M)) {
    /* Outport: '<Root>/signal1' */
    SYSTEM0_Y.signal1 = rtb_DataTypeConversion;

    /* Outport: '<Root>/signal2' */
    SYSTEM0_Y.signal2 = rtb_DataTypeConversion1;

    /* Outport: '<Root>/signal3' */
    SYSTEM0_Y.signal3 = rtb_DataTypeConversion2;

    /* Outport: '<Root>/signal4' */
    SYSTEM0_Y.signal4 = rtb_DataTypeConversion3;

    /* Outport: '<Root>/signal5' */
    SYSTEM0_Y.signal5 = rtb_DataTypeConversion4;

    /* Outport: '<Root>/signal6' */
    SYSTEM0_Y.signal6 = rtb_DataTypeConversion5;

    /* Outport: '<Root>/signal7' */
    SYSTEM0_Y.signal7 = rtb_DataTypeConversion6;

    /* Saturate: '<S60>/Saturation' incorporates:
     *  ZeroOrderHold: '<S60>/Zero-Order Hold1'
     */
    if (SYSTEM0_B.Sum2[0] > 2.0) {
      /* Outport: '<Root>/signal8' incorporates:
       *  DataTypeConversion: '<S2>/Data Type Conversion7'
       */
      SYSTEM0_Y.signal8 = 2.0F;
    } else if (SYSTEM0_B.Sum2[0] < -2.0) {
      /* Outport: '<Root>/signal8' incorporates:
       *  DataTypeConversion: '<S2>/Data Type Conversion7'
       */
      SYSTEM0_Y.signal8 = -2.0F;
    } else {
      /* Outport: '<Root>/signal8' incorporates:
       *  DataTypeConversion: '<S2>/Data Type Conversion7'
       */
      SYSTEM0_Y.signal8 = (real32_T)SYSTEM0_B.Sum2[0];
    }

    if (SYSTEM0_B.Sum2[1] > 2.0) {
      /* Outport: '<Root>/signal9' incorporates:
       *  DataTypeConversion: '<S2>/Data Type Conversion8'
       */
      SYSTEM0_Y.signal9 = 2.0F;
    } else if (SYSTEM0_B.Sum2[1] < -2.0) {
      /* Outport: '<Root>/signal9' incorporates:
       *  DataTypeConversion: '<S2>/Data Type Conversion8'
       */
      SYSTEM0_Y.signal9 = -2.0F;
    } else {
      /* Outport: '<Root>/signal9' incorporates:
       *  DataTypeConversion: '<S2>/Data Type Conversion8'
       */
      SYSTEM0_Y.signal9 = (real32_T)SYSTEM0_B.Sum2[1];
    }

    if (SYSTEM0_B.Sum2[2] > 2.0) {
      /* Outport: '<Root>/signal10' incorporates:
       *  DataTypeConversion: '<S2>/Data Type Conversion9'
       */
      SYSTEM0_Y.signal10 = 2.0F;
    } else if (SYSTEM0_B.Sum2[2] < -2.0) {
      /* Outport: '<Root>/signal10' incorporates:
       *  DataTypeConversion: '<S2>/Data Type Conversion9'
       */
      SYSTEM0_Y.signal10 = -2.0F;
    } else {
      /* Outport: '<Root>/signal10' incorporates:
       *  DataTypeConversion: '<S2>/Data Type Conversion9'
       */
      SYSTEM0_Y.signal10 = (real32_T)SYSTEM0_B.Sum2[2];
    }

    /* End of Saturate: '<S60>/Saturation' */

    /* Outport: '<Root>/signal11' */
    SYSTEM0_Y.signal11 = rtb_DataTypeConversion10;

    /* Outport: '<Root>/signal12' */
    SYSTEM0_Y.signal12 = rtb_DataTypeConversion11;

    /* Outport: '<Root>/signal13' */
    SYSTEM0_Y.signal13 = rtb_DataTypeConversion12;

    /* Outport: '<Root>/signal14' */
    SYSTEM0_Y.signal14 = rtb_DataTypeConversion13;

    /* Outport: '<Root>/signal15' */
    SYSTEM0_Y.signal15 = rtb_DataTypeConversion14;

    /* Outport: '<Root>/signal16' */
    SYSTEM0_Y.signal16 = rtb_DataTypeConversion15;

    /* Outport: '<Root>/signal17' */
    SYSTEM0_Y.signal17 = rtb_DataTypeConversion16;

    /* Outport: '<Root>/signal18' */
    SYSTEM0_Y.signal18 = rtb_DataTypeConversion19;

    /* Outport: '<Root>/signal19' */
    SYSTEM0_Y.signal19 = rtb_DataTypeConversion17_idx_0;

    /* Outport: '<Root>/signal20' */
    SYSTEM0_Y.signal20 = rtb_DataTypeConversion17_idx_1;

    /* Outport: '<Root>/signal21' */
    SYSTEM0_Y.signal21 = rtb_DataTypeConversion17_idx_2;

    /* Outport: '<Root>/signal22' */
    SYSTEM0_Y.signal22 = rtb_DataTypeConversion18;

    /* SignalConversion: '<S7>/TmpSignal ConversionAt SFunction Inport1' incorporates:
     *  Inport: '<Root>/In1'
     *  Inport: '<Root>/In2'
     *  Inport: '<Root>/In3'
     *  Inport: '<Root>/In4'
     *  MATLAB Function: '<S4>/LiPo Battery'
     */
    SYSTEM0_B.voltage[0] = SYSTEM0_U.In1;
    SYSTEM0_B.voltage[1] = SYSTEM0_U.In2;
    SYSTEM0_B.voltage[2] = SYSTEM0_U.In3;
    SYSTEM0_B.voltage[3] = SYSTEM0_U.In4;

    /* MATLAB Function: '<S4>/LiPo Battery' */
    /* MATLAB Function 'SYSTEM/dynamic_model/Dynamics/LiPo Battery': '<S7>:1' */
    /* ================================Constants================================= */
    /* '<S7>:1:6' */
    /*  maximum current per motor             (A) */
    /* '<S7>:1:7' */
    /*  current of auxillary components       (A)         */
    /*  LiPo current capacity                 (Ah)                                  */
    /*  LiPo cell count */
    /* ========================================================================== */
    /* '<S7>:1:19' */
    /* '<S7>:1:20' */
    SYSTEM0_DW.discharge += ((((SYSTEM0_B.voltage[0] * SYSTEM0_B.voltage[0] *
      5.5 + SYSTEM0_B.voltage[1] * SYSTEM0_B.voltage[1] * 5.5) +
      SYSTEM0_B.voltage[2] * SYSTEM0_B.voltage[2] * 5.5) + SYSTEM0_B.voltage[3] *
      SYSTEM0_B.voltage[3] * 5.5) + 0.1) * 1.111111111111111E-6;

    /* '<S7>:1:22' */
    if ((0.0 < SYSTEM0_DW.discharge) && (SYSTEM0_DW.discharge <= 0.2)) {
      /* '<S7>:1:24' */
      /* '<S7>:1:25' */
      rtb_IntegratorSecondOrder_o2_c1 = ((SYSTEM0_DW.discharge *
        SYSTEM0_DW.discharge * 16.975 + -14.029 * pow(SYSTEM0_DW.discharge, 3.0))
        - 5.3339 * SYSTEM0_DW.discharge) + 4.2;
    } else if ((0.2 < SYSTEM0_DW.discharge) && (SYSTEM0_DW.discharge < 0.7)) {
      /* '<S7>:1:26' */
      /* '<S7>:1:27' */
      rtb_IntegratorSecondOrder_o2_c1 = -0.2 * SYSTEM0_DW.discharge + 3.74;
    } else {
      /* '<S7>:1:29' */
      rtb_IntegratorSecondOrder_o2_c1 = ((SYSTEM0_DW.discharge *
        SYSTEM0_DW.discharge * 89.6 + -48.0 * pow(SYSTEM0_DW.discharge, 3.0)) -
        55.08 * SYSTEM0_DW.discharge) + 14.716;
    }

    if (rtb_IntegratorSecondOrder_o2_c1 < 2.5) {
      /* '<S7>:1:32' */
      /* '<S7>:1:33' */
      rtb_IntegratorSecondOrder_o2_c1 = 0.0;
    }

    /* '<S7>:1:36' */
    rtb_IntegratorSecondOrder_o2_c1 *= 2.0;

    /* '<S7>:1:37' */
    SYSTEM0_B.voltage[0] *= rtb_IntegratorSecondOrder_o2_c1;
    SYSTEM0_B.voltage[1] *= rtb_IntegratorSecondOrder_o2_c1;
    SYSTEM0_B.voltage[2] *= rtb_IntegratorSecondOrder_o2_c1;
    SYSTEM0_B.voltage[3] *= rtb_IntegratorSecondOrder_o2_c1;

    /* ========================================================================== */
  }

  /* Product: '<S42>/j x k' */
  rtb_Product2_pl = SYSTEM0_X.ubvbwb_CSTATE[1] * SYSTEM0_B.pqr[2];

  /* Sum: '<S6>/Sum' incorporates:
   *  Product: '<S42>/i x j'
   *  Product: '<S42>/k x i'
   *  Product: '<S43>/i x k'
   *  Product: '<S43>/j x i'
   *  Product: '<S43>/k x j'
   *  Sum: '<S13>/Sum'
   */
  SYSTEM0_B.Sum[0] = (rtb_Product2_pl - SYSTEM0_X.ubvbwb_CSTATE[2] *
                      SYSTEM0_B.pqr[1]) + SYSTEM0_B.Product[0];
  SYSTEM0_B.Sum[1] = (SYSTEM0_X.ubvbwb_CSTATE[2] * SYSTEM0_B.pqr[0] -
                      SYSTEM0_X.ubvbwb_CSTATE[0] * SYSTEM0_B.pqr[2]) +
    SYSTEM0_B.Product[1];
  SYSTEM0_B.Sum[2] = (SYSTEM0_X.ubvbwb_CSTATE[0] * SYSTEM0_B.pqr[1] -
                      SYSTEM0_X.ubvbwb_CSTATE[1] * SYSTEM0_B.pqr[0]) +
    SYSTEM0_B.Product[2];

  /* DotProduct: '<S20>/Dot Product' */
  rtb_Product2_pl = ((rtb_q0q1q2q3[0] * rtb_q0q1q2q3[0] + rtb_q0q1q2q3[1] *
                      rtb_q0q1q2q3[1]) + rtb_q0q1q2q3[2] * rtb_q0q1q2q3[2]) +
    rtb_q0q1q2q3[3] * rtb_q0q1q2q3[3];

  /* Sum: '<S20>/Sum' incorporates:
   *  Constant: '<S20>/Constant'
   */
  rtb_Product2_pl = 1.0 - rtb_Product2_pl;

  /* Fcn: '<S20>/q0dot' */
  SYSTEM0_B.q0dot = ((rtb_q0q1q2q3[1] * SYSTEM0_B.pqr[0] + rtb_q0q1q2q3[2] *
                      SYSTEM0_B.pqr[1]) + rtb_q0q1q2q3[3] * SYSTEM0_B.pqr[2]) *
    -0.5 + rtb_Product2_pl * rtb_q0q1q2q3[0];

  /* Fcn: '<S20>/q1dot' */
  SYSTEM0_B.q1dot = ((rtb_q0q1q2q3[0] * SYSTEM0_B.pqr[0] + rtb_q0q1q2q3[2] *
                      SYSTEM0_B.pqr[2]) - rtb_q0q1q2q3[3] * SYSTEM0_B.pqr[1]) *
    0.5 + rtb_Product2_pl * rtb_q0q1q2q3[1];

  /* Fcn: '<S20>/q2dot' */
  SYSTEM0_B.q2dot = ((rtb_q0q1q2q3[0] * SYSTEM0_B.pqr[1] + rtb_q0q1q2q3[3] *
                      SYSTEM0_B.pqr[0]) - rtb_q0q1q2q3[1] * SYSTEM0_B.pqr[2]) *
    0.5 + rtb_Product2_pl * rtb_q0q1q2q3[2];

  /* Fcn: '<S20>/q3dot' */
  SYSTEM0_B.q3dot = ((rtb_q0q1q2q3[0] * SYSTEM0_B.pqr[2] + rtb_q0q1q2q3[1] *
                      SYSTEM0_B.pqr[1]) - rtb_q0q1q2q3[2] * SYSTEM0_B.pqr[0]) *
    0.5 + rtb_Product2_pl * rtb_q0q1q2q3[3];

  /* Product: '<S16>/Product' incorporates:
   *  Math: '<S6>/Transpose'
   */
  for (i = 0; i < 3; i++) {
    SYSTEM0_B.Product_j[i] = 0.0;
    SYSTEM0_B.Product_j[i] += rtb_VectorConcatenate[3 * i] *
      SYSTEM0_X.ubvbwb_CSTATE[0];
    SYSTEM0_B.Product_j[i] += rtb_VectorConcatenate[3 * i + 1] *
      SYSTEM0_X.ubvbwb_CSTATE[1];
    SYSTEM0_B.Product_j[i] += rtb_VectorConcatenate[3 * i + 2] *
      SYSTEM0_X.ubvbwb_CSTATE[2];
  }

  /* End of Product: '<S16>/Product' */

  /* Sum: '<S49>/Sum2' incorporates:
   *  Gain: '<S49>/2*zeta * wn'
   *  Gain: '<S49>/wn^2'
   *  Inport: '<Root>/In1'
   *  Product: '<S9>/Product'
   *  SecondOrderIntegrator: '<S49>/Integrator, Second-Order'
   *  Sum: '<S49>/Sum3'
   */
  SYSTEM0_B.Sum2_n = (SYSTEM0_U.In1 * SYSTEM0_B.voltage[0] -
                      SYSTEM0_X.IntegratorSecondOrder_CSTATE[0]) * 4900.0 -
    140.0 * SYSTEM0_X.IntegratorSecondOrder_CSTATE[1];

  /* Sum: '<S50>/Sum2' incorporates:
   *  Gain: '<S50>/2*zeta * wn'
   *  Gain: '<S50>/wn^2'
   *  Inport: '<Root>/In2'
   *  Product: '<S9>/Product'
   *  SecondOrderIntegrator: '<S50>/Integrator, Second-Order'
   *  Sum: '<S50>/Sum3'
   */
  SYSTEM0_B.Sum2_j = (SYSTEM0_U.In2 * SYSTEM0_B.voltage[1] -
                      SYSTEM0_X.IntegratorSecondOrder_CSTATE_h[0]) * 4900.0 -
    140.0 * SYSTEM0_X.IntegratorSecondOrder_CSTATE_h[1];

  /* Sum: '<S51>/Sum2' incorporates:
   *  Gain: '<S51>/2*zeta * wn'
   *  Gain: '<S51>/wn^2'
   *  Inport: '<Root>/In3'
   *  Product: '<S9>/Product'
   *  SecondOrderIntegrator: '<S51>/Integrator, Second-Order'
   *  Sum: '<S51>/Sum3'
   */
  SYSTEM0_B.Sum2_c = (SYSTEM0_U.In3 * SYSTEM0_B.voltage[2] -
                      SYSTEM0_X.IntegratorSecondOrder_CSTATE_n[0]) * 4900.0 -
    140.0 * SYSTEM0_X.IntegratorSecondOrder_CSTATE_n[1];

  /* Sum: '<S52>/Sum2' incorporates:
   *  Gain: '<S52>/2*zeta * wn'
   *  Gain: '<S52>/wn^2'
   *  Inport: '<Root>/In4'
   *  Product: '<S9>/Product'
   *  SecondOrderIntegrator: '<S52>/Integrator, Second-Order'
   *  Sum: '<S52>/Sum3'
   */
  SYSTEM0_B.Sum2_p = (SYSTEM0_U.In4 * SYSTEM0_B.voltage[3] -
                      SYSTEM0_X.IntegratorSecondOrder_CSTATE_d[0]) * 4900.0 -
    140.0 * SYSTEM0_X.IntegratorSecondOrder_CSTATE_d[1];

  /* RateTransition: '<S57>/Rate Transition1' */
  if (rtmIsMajorTimeStep(SYSTEM0_M) && SYSTEM0_M->Timing.RateInteraction.TID0_3 ==
      1) {
    SYSTEM0_B.RateTransition1_a[0] = SYSTEM0_X.xeyeze_CSTATE[0];
    SYSTEM0_B.RateTransition1_a[1] = SYSTEM0_X.xeyeze_CSTATE[1];
    SYSTEM0_B.RateTransition1_a[2] = SYSTEM0_X.xeyeze_CSTATE[2];
  }

  /* End of RateTransition: '<S57>/Rate Transition1' */

  /* RateTransition: '<S57>/Rate Transition' incorporates:
   *  Constant: '<S5>/h_ref'
   */
  if (rtmIsMajorTimeStep(SYSTEM0_M) && (SYSTEM0_M->Timing.RateInteraction.TID1_3
       == 1)) {
    SYSTEM0_B.RateTransition = 0.0;
  }

  /* End of RateTransition: '<S57>/Rate Transition' */

  /* RateTransition: '<S57>/Rate Transition2' */
  if (rtmIsMajorTimeStep(SYSTEM0_M) && SYSTEM0_M->Timing.RateInteraction.TID0_3 ==
      1) {
    SYSTEM0_B.RateTransition2[0] = SYSTEM0_B.Product_j[0];
    SYSTEM0_B.RateTransition2[1] = SYSTEM0_B.Product_j[1];
    SYSTEM0_B.RateTransition2[2] = SYSTEM0_B.Product_j[2];
  }

  /* End of RateTransition: '<S57>/Rate Transition2' */

  /* Gain: '<S58>/1//T0' */
  rtb_Product2_pl = 0.00347041471455839 * SYSTEM0_B.Sum1_e;

  /* Math: '<S58>/(T//T0)^(g//LR) ' */
  if (rtb_Product2_pl < 0.0) {
    rtb_IntegratorSecondOrder_o1_c = -pow(-rtb_Product2_pl, 5.2558756014667134);
  } else {
    rtb_IntegratorSecondOrder_o1_c = pow(rtb_Product2_pl, 5.2558756014667134);
  }

  /* End of Math: '<S58>/(T//T0)^(g//LR) ' */

  /* Product: '<S58>/Product' */
  rtb_IntegratorSecondOrder_o2_c1 = rtb_IntegratorSecondOrder_o1_c /
    rtb_Product2_pl;

  /* Sum: '<S58>/Sum' incorporates:
   *  Constant: '<S58>/Altitude of Troposphere'
   */
  rtb_Product2_pl = 11000.0 - SYSTEM0_B.Sum1;

  /* Saturate: '<S58>/Limit  altitude  to Stratosphere' */
  if (rtb_Product2_pl > 0.0) {
    rtb_Product2_pl = 0.0;
  } else {
    if (rtb_Product2_pl < -9000.0) {
      rtb_Product2_pl = -9000.0;
    }
  }

  /* Math: '<S58>/Stratosphere Model' incorporates:
   *  Gain: '<S58>/g//R'
   *  Product: '<S58>/Product1'
   *  Saturate: '<S58>/Limit  altitude  to Stratosphere'
   *
   * About '<S58>/Stratosphere Model':
   *  Operator: exp
   */
  rtb_Product2_pl = exp(1.0 / SYSTEM0_B.Sum1_e * (0.034163191409533639 *
    rtb_Product2_pl));

  /* Product: '<S58>/Product3' incorporates:
   *  Gain: '<S58>/rho0'
   */
  SYSTEM0_B.Product3 = 1.225 * rtb_IntegratorSecondOrder_o2_c1 * rtb_Product2_pl;

  /* Product: '<S58>/Product2' incorporates:
   *  Gain: '<S58>/P0'
   */
  rtb_Product2_pl *= 101325.0 * rtb_IntegratorSecondOrder_o1_c;

  /* RateTransition: '<S61>/Rate Transition' */
  if (rtmIsMajorTimeStep(SYSTEM0_M) && SYSTEM0_M->Timing.RateInteraction.TID0_2 ==
      1) {
    SYSTEM0_B.RateTransition_b = SYSTEM0_B.Product3;

    /* RateTransition: '<S61>/Rate Transition1' */
    SYSTEM0_B.RateTransition1_p[0] = SYSTEM0_X.ubvbwb_CSTATE[0];
    SYSTEM0_B.RateTransition1_p[1] = SYSTEM0_X.ubvbwb_CSTATE[1];
    SYSTEM0_B.RateTransition1_p[2] = SYSTEM0_X.ubvbwb_CSTATE[2];

    /* RateTransition: '<S63>/Rate Transition' */
    SYSTEM0_B.RateTransition_c = rtb_Product2_pl;
  }

  /* End of RateTransition: '<S61>/Rate Transition' */

  /* If: '<S55>/If' incorporates:
   *  If: '<S66>/Find Maximum Diagonal Value'
   *  Sum: '<S68>/Add'
   */
  if (rtmIsMajorTimeStep(SYSTEM0_M)) {
    SYSTEM0_DW.If_ActiveSubsystem = (int8_T)!((rtb_VectorConcatenate[0] +
      rtb_VectorConcatenate[4]) + rtb_VectorConcatenate[8] > 0.0);
  }

  /* Outputs for IfAction SubSystem: '<S55>/Negative Trace' incorporates:
   *  ActionPort: '<S66>/Action Port'
   */
  if ((SYSTEM0_DW.If_ActiveSubsystem == 1) && rtmIsMajorTimeStep(SYSTEM0_M)) {
    if ((rtb_VectorConcatenate[4] > rtb_VectorConcatenate[0]) &&
        (rtb_VectorConcatenate[4] > rtb_VectorConcatenate[8])) {
      SYSTEM0_DW.FindMaximumDiagonalValue_ActiveSubsystem = 0;
    } else if (rtb_VectorConcatenate[8] > rtb_VectorConcatenate[0]) {
      SYSTEM0_DW.FindMaximumDiagonalValue_ActiveSubsystem = 1;
    } else {
      SYSTEM0_DW.FindMaximumDiagonalValue_ActiveSubsystem = 2;
    }
  }

  /* End of If: '<S55>/If' */
  /* End of Outputs for SubSystem: '<S55>/Negative Trace' */
  if (rtmIsMajorTimeStep(SYSTEM0_M)) {
    /* Update for Integrator: '<S10>/q0 q1 q2 q3' */
    SYSTEM0_DW.q0q1q2q3_IWORK.IcNeedsLoading = 0;
    if (rtmIsMajorTimeStep(SYSTEM0_M)) {
      /* Update for Memory: '<S4>/Memory2' */
      SYSTEM0_DW.Memory2_PreviousInput = SYSTEM0_B.Product3;

      /* Update for RandomNumber: '<S132>/White Noise' */
      SYSTEM0_DW.NextOutput[0] = rt_nrand_Upu32_Yd_f_pw(&SYSTEM0_DW.RandSeed[0]);
      SYSTEM0_DW.NextOutput[1] = rt_nrand_Upu32_Yd_f_pw(&SYSTEM0_DW.RandSeed[1]);
      SYSTEM0_DW.NextOutput[2] = rt_nrand_Upu32_Yd_f_pw(&SYSTEM0_DW.RandSeed[2]);

      /* Update for RandomNumber: '<S53>/Random Number' */
      SYSTEM0_DW.NextOutput_n = rt_nrand_Upu32_Yd_f_pw(&SYSTEM0_DW.RandSeed_f) *
        0.01;

      /* Update for RandomNumber: '<S65>/Random Number' */
      SYSTEM0_DW.NextOutput_h = rt_nrand_Upu32_Yd_f_pw(&SYSTEM0_DW.RandSeed_d) *
        0.01;

      /* Update for RandomNumber: '<S149>/White Noise' */
      SYSTEM0_DW.NextOutput_p[0] = rt_nrand_Upu32_Yd_f_pw
        (&SYSTEM0_DW.RandSeed_j[0]);
      SYSTEM0_DW.NextOutput_p[1] = rt_nrand_Upu32_Yd_f_pw
        (&SYSTEM0_DW.RandSeed_j[1]);
      SYSTEM0_DW.NextOutput_p[2] = rt_nrand_Upu32_Yd_f_pw
        (&SYSTEM0_DW.RandSeed_j[2]);

      /* Update for RandomNumber: '<S60>/Random Number' */
      SYSTEM0_DW.NextOutput_a = rt_nrand_Upu32_Yd_f_pw(&SYSTEM0_DW.RandSeed_p) *
        0.0031622776601683794;
    }
  }                                    /* end MajorTimeStep */

  if (rtmIsMajorTimeStep(SYSTEM0_M)) {
    rt_ertODEUpdateContinuousStates(&SYSTEM0_M->solverInfo);

    /* Update absolute time */
    /* The "clockTick0" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick0"
     * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick0 and the high bits
     * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++SYSTEM0_M->Timing.clockTick0)) {
      ++SYSTEM0_M->Timing.clockTickH0;
    }

    SYSTEM0_M->Timing.t[0] = rtsiGetSolverStopTime(&SYSTEM0_M->solverInfo);

    /* Update absolute time */
    /* The "clockTick1" counts the number of times the code of this task has
     * been executed. The resolution of this integer timer is 0.004, which is the step size
     * of the task. Size of "clockTick1" ensures timer will not overflow during the
     * application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick1 and the high bits
     * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
     */
    SYSTEM0_M->Timing.clockTick1++;
    if (!SYSTEM0_M->Timing.clockTick1) {
      SYSTEM0_M->Timing.clockTickH1++;
    }
  }                                    /* end MajorTimeStep */
}

/* Derivatives for root system: '<Root>' */
void SYSTEM0_derivatives(void)
{
  XDot_SYSTEM0_T *_rtXdot;
  _rtXdot = ((XDot_SYSTEM0_T *) SYSTEM0_M->ModelData.derivs);

  /* Derivatives for Integrator: '<S10>/q0 q1 q2 q3' */
  {
    ((XDot_SYSTEM0_T *) SYSTEM0_M->ModelData.derivs)->q0q1q2q3_CSTATE[0] =
      SYSTEM0_B.q0dot;
    ((XDot_SYSTEM0_T *) SYSTEM0_M->ModelData.derivs)->q0q1q2q3_CSTATE[1] =
      SYSTEM0_B.q1dot;
    ((XDot_SYSTEM0_T *) SYSTEM0_M->ModelData.derivs)->q0q1q2q3_CSTATE[2] =
      SYSTEM0_B.q2dot;
    ((XDot_SYSTEM0_T *) SYSTEM0_M->ModelData.derivs)->q0q1q2q3_CSTATE[3] =
      SYSTEM0_B.q3dot;
  }

  /* Derivatives for Integrator: '<S6>/xe,ye,ze' */
  _rtXdot->xeyeze_CSTATE[0] = SYSTEM0_B.Product_j[0];
  _rtXdot->xeyeze_CSTATE[1] = SYSTEM0_B.Product_j[1];
  _rtXdot->xeyeze_CSTATE[2] = SYSTEM0_B.Product_j[2];

  /* Derivatives for Integrator: '<S6>/ub,vb,wb' */
  _rtXdot->ubvbwb_CSTATE[0] = SYSTEM0_B.Sum[0];
  _rtXdot->ubvbwb_CSTATE[1] = SYSTEM0_B.Sum[1];
  _rtXdot->ubvbwb_CSTATE[2] = SYSTEM0_B.Sum[2];

  /* Derivatives for SecondOrderIntegrator: '<S49>/Integrator, Second-Order' */
  if (SYSTEM0_DW.IntegratorSecondOrder_MODE == 0) {
    _rtXdot->IntegratorSecondOrder_CSTATE[0] =
      SYSTEM0_X.IntegratorSecondOrder_CSTATE[1];
    _rtXdot->IntegratorSecondOrder_CSTATE[1] = SYSTEM0_B.Sum2_n;
  }

  /* End of Derivatives for SecondOrderIntegrator: '<S49>/Integrator, Second-Order' */

  /* Derivatives for SecondOrderIntegrator: '<S50>/Integrator, Second-Order' */
  if (SYSTEM0_DW.IntegratorSecondOrder_MODE_p == 0) {
    _rtXdot->IntegratorSecondOrder_CSTATE_h[0] =
      SYSTEM0_X.IntegratorSecondOrder_CSTATE_h[1];
    _rtXdot->IntegratorSecondOrder_CSTATE_h[1] = SYSTEM0_B.Sum2_j;
  }

  /* End of Derivatives for SecondOrderIntegrator: '<S50>/Integrator, Second-Order' */

  /* Derivatives for SecondOrderIntegrator: '<S51>/Integrator, Second-Order' */
  if (SYSTEM0_DW.IntegratorSecondOrder_MODE_a == 0) {
    _rtXdot->IntegratorSecondOrder_CSTATE_n[0] =
      SYSTEM0_X.IntegratorSecondOrder_CSTATE_n[1];
    _rtXdot->IntegratorSecondOrder_CSTATE_n[1] = SYSTEM0_B.Sum2_c;
  }

  /* End of Derivatives for SecondOrderIntegrator: '<S51>/Integrator, Second-Order' */

  /* Derivatives for SecondOrderIntegrator: '<S52>/Integrator, Second-Order' */
  if (SYSTEM0_DW.IntegratorSecondOrder_MODE_pu == 0) {
    _rtXdot->IntegratorSecondOrder_CSTATE_d[0] =
      SYSTEM0_X.IntegratorSecondOrder_CSTATE_d[1];
    _rtXdot->IntegratorSecondOrder_CSTATE_d[1] = SYSTEM0_B.Sum2_p;
  }

  /* End of Derivatives for SecondOrderIntegrator: '<S52>/Integrator, Second-Order' */

  /* Derivatives for Integrator: '<S6>/p,q,r ' */
  _rtXdot->pqr_CSTATE[0] = SYSTEM0_B.Product2[0];
  _rtXdot->pqr_CSTATE[1] = SYSTEM0_B.Product2[1];
  _rtXdot->pqr_CSTATE[2] = SYSTEM0_B.Product2[2];
}

/* Model step function for TID2 */
void SYSTEM0_step2(void)               /* Sample time: [0.04s, 0.0s] */
{
  /* local block i/o variables */
  real_T rtb_Saturation;
  real_T rtb_Saturation_h;

  /* Sum: '<S61>/Sum2' incorporates:
   *  DotProduct: '<S61>/Dot Product'
   *  Gain: '<S61>/Bar2mBar'
   *  Gain: '<S61>/Gain'
   *  Gain: '<S61>/Pa2Bar'
   *  Product: '<S61>/Product'
   *  RandomNumber: '<S61>/Random Number'
   */
  rtb_Saturation = ((SYSTEM0_B.RateTransition1_p[0] *
                     SYSTEM0_B.RateTransition1_p[0] +
                     SYSTEM0_B.RateTransition1_p[1] *
                     SYSTEM0_B.RateTransition1_p[1]) +
                    SYSTEM0_B.RateTransition1_p[2] *
                    SYSTEM0_B.RateTransition1_p[2]) * SYSTEM0_B.RateTransition_b
    * 0.5 * 1.0E-5 * 1000.0 + SYSTEM0_DW.NextOutput_l;

  /* Saturate: '<S61>/Saturation' */
  if (rtb_Saturation > 1000.0) {
    /* Sum: '<S61>/Sum2' */
    rtb_Saturation = 1000.0;
  } else {
    if (rtb_Saturation < 0.0) {
      /* Sum: '<S61>/Sum2' */
      rtb_Saturation = 0.0;
    }
  }

  /* End of Saturate: '<S61>/Saturation' */

  /* Sum: '<S63>/Sum2' incorporates:
   *  Gain: '<S63>/Bar2mBar'
   *  Gain: '<S63>/Pa2Bar'
   *  RandomNumber: '<S63>/Random Number'
   */
  rtb_Saturation_h = 1.0E-5 * SYSTEM0_B.RateTransition_c * 1000.0 +
    SYSTEM0_DW.NextOutput_lc;

  /* Saturate: '<S63>/Saturation' */
  if (rtb_Saturation_h > 1200.0) {
    /* Sum: '<S63>/Sum2' */
    rtb_Saturation_h = 1200.0;
  } else {
    if (rtb_Saturation_h < 10.0) {
      /* Sum: '<S63>/Sum2' */
      rtb_Saturation_h = 10.0;
    }
  }

  /* End of Saturate: '<S63>/Saturation' */

  /* Update for RateTransition: '<S3>/Rate Transition1' */
  SYSTEM0_DW.RateTransition1_Buffer0 = rtb_Saturation_h;

  /* Update for RateTransition: '<S3>/Rate Transition3' */
  SYSTEM0_DW.RateTransition3_Buffer0 = rtb_Saturation;

  /* Update for RandomNumber: '<S61>/Random Number' */
  SYSTEM0_DW.NextOutput_l = rt_nrand_Upu32_Yd_f_pw(&SYSTEM0_DW.RandSeed_fw) *
    0.0031622776601683794;

  /* Update for RandomNumber: '<S63>/Random Number' */
  SYSTEM0_DW.NextOutput_lc = rt_nrand_Upu32_Yd_f_pw(&SYSTEM0_DW.RandSeed_h) *
    0.0031622776601683794;
}

/* Model step function for TID3 */
void SYSTEM0_step3(void)               /* Sample time: [1.0s, 0.0s] */
{
  /* local block i/o variables */
  real_T rtb_Switch_b;
  real_T rtb_Switch_m;
  real_T rtb_Sum1_g;
  real_T rtb_Add2[3];
  real_T rtb_TrigonometricFunction;
  real_T rtb_Abs1;
  int16_T rtb_Compare_0;

  /* RandomNumber: '<S57>/Random Number' */
  rtb_Add2[0] = SYSTEM0_DW.NextOutput_o[0];
  rtb_Add2[1] = SYSTEM0_DW.NextOutput_o[1];
  rtb_Add2[2] = SYSTEM0_DW.NextOutput_o[2];

  /* Sum: '<S57>/Add1' */
  rtb_Add2[0] += SYSTEM0_B.RateTransition1_a[0];
  rtb_Add2[1] += SYSTEM0_B.RateTransition1_a[1];
  rtb_Add2[2] += SYSTEM0_B.RateTransition1_a[2];

  /* Sum: '<S109>/Sum' incorporates:
   *  Gain: '<S113>/Unit Conversion'
   *  Product: '<S112>/rad lat'
   *  Product: '<S112>/x*cos'
   */
  rtb_Switch_b = rtb_Add2[0] * 1.5784225029068334E-7 * 57.295779513082323;

  /* Switch: '<S117>/Switch' incorporates:
   *  Abs: '<S117>/Abs'
   */
  if (fabs(rtb_Switch_b) > 180.0) {
    /* Sum: '<S109>/Sum' incorporates:
     *  Bias: '<S117>/Bias'
     *  Bias: '<S117>/Bias1'
     *  Constant: '<S117>/Constant2'
     *  Math: '<S117>/Math Function1'
     */
    rtb_Switch_b = rt_modd(rtb_Switch_b + 180.0, 360.0) + -180.0;
  }

  /* End of Switch: '<S117>/Switch' */

  /* Abs: '<S114>/Abs1' */
  rtb_Abs1 = fabs(rtb_Switch_b);

  /* Switch: '<S114>/Switch' incorporates:
   *  Constant: '<S110>/Constant'
   *  Constant: '<S110>/Constant1'
   *  Constant: '<S116>/Constant'
   *  RelationalOperator: '<S116>/Compare'
   *  Switch: '<S110>/Switch1'
   */
  if ((rtb_Abs1 > 90.0) > 0) {
    /* Signum: '<S114>/Sign1' */
    if (rtb_Switch_b < 0.0) {
      /* Sum: '<S109>/Sum' */
      rtb_Switch_b = -1.0;
    } else {
      if (rtb_Switch_b > 0.0) {
        /* Sum: '<S109>/Sum' */
        rtb_Switch_b = 1.0;
      }
    }

    /* Sum: '<S109>/Sum' incorporates:
     *  Bias: '<S114>/Bias'
     *  Bias: '<S114>/Bias1'
     *  Gain: '<S114>/Gain'
     *  Product: '<S114>/Divide1'
     *  Signum: '<S114>/Sign1'
     */
    rtb_Switch_b *= -(rtb_Abs1 + -90.0) + 90.0;
    rtb_Compare_0 = 180;
  } else {
    rtb_Compare_0 = 0;
  }

  /* End of Switch: '<S114>/Switch' */

  /* Sum: '<S110>/Sum' incorporates:
   *  Gain: '<S113>/Unit Conversion'
   *  Product: '<S112>/rad long '
   *  Product: '<S112>/y*cos'
   *  Sum: '<S109>/Sum'
   */
  rtb_Switch_m = (1.5678559428873849E-7 * rtb_Add2[1] * 57.295779513082323 +
                  SYSTEM0_ConstB.Switch_d) + (real_T)rtb_Compare_0;

  /* Switch: '<S115>/Switch' incorporates:
   *  Abs: '<S115>/Abs'
   *  Bias: '<S115>/Bias'
   *  Bias: '<S115>/Bias1'
   *  Constant: '<S115>/Constant2'
   *  Math: '<S115>/Math Function1'
   */
  if (fabs(rtb_Switch_m) > 180.0) {
    rtb_Switch_m = rt_modd(rtb_Switch_m + 180.0, 360.0) + -180.0;
  }

  /* End of Switch: '<S115>/Switch' */

  /* Sum: '<S109>/Sum1' incorporates:
   *  UnaryMinus: '<S109>/Ze2height'
   */
  rtb_Sum1_g = -rtb_Add2[2] - SYSTEM0_B.RateTransition;

  /* RandomNumber: '<S57>/Random Number1' */
  rtb_Add2[0] = SYSTEM0_DW.NextOutput_hg[0];
  rtb_Add2[1] = SYSTEM0_DW.NextOutput_hg[1];
  rtb_Add2[2] = SYSTEM0_DW.NextOutput_hg[2];

  /* Sum: '<S57>/Add2' */
  rtb_Add2[0] += SYSTEM0_B.RateTransition2[0];
  rtb_Add2[1] += SYSTEM0_B.RateTransition2[1];
  rtb_Add2[2] += SYSTEM0_B.RateTransition2[2];

  /* Trigonometry: '<S57>/Trigonometric Function' */
  rtb_TrigonometricFunction = atan2(rtb_Add2[1], rtb_Add2[0]);

  /* Update for RateTransition: '<S3>/Rate Transition6' */
  SYSTEM0_DW.RateTransition6_Buffer0 = rtb_Switch_m;

  /* Update for RateTransition: '<S3>/Rate Transition5' */
  SYSTEM0_DW.RateTransition5_Buffer0 = rtb_Switch_b;

  /* Update for RateTransition: '<S3>/Rate Transition7' */
  SYSTEM0_DW.RateTransition7_Buffer0 = rtb_Sum1_g;

  /* Update for RateTransition: '<S3>/Rate Transition8' */
  SYSTEM0_DW.RateTransition8_Buffer0[0] = rtb_Add2[0];
  SYSTEM0_DW.RateTransition8_Buffer0[1] = rtb_Add2[1];
  SYSTEM0_DW.RateTransition8_Buffer0[2] = rtb_Add2[2];

  /* Update for RateTransition: '<S3>/Rate Transition10' */
  SYSTEM0_DW.RateTransition10_Buffer0 = rtb_TrigonometricFunction;

  /* Update for RandomNumber: '<S57>/Random Number' */
  SYSTEM0_DW.NextOutput_o[0] = rt_nrand_Upu32_Yd_f_pw(&SYSTEM0_DW.RandSeed_e[0])
    * 15.0;
  SYSTEM0_DW.NextOutput_o[1] = rt_nrand_Upu32_Yd_f_pw(&SYSTEM0_DW.RandSeed_e[1])
    * 15.0;
  SYSTEM0_DW.NextOutput_o[2] = rt_nrand_Upu32_Yd_f_pw(&SYSTEM0_DW.RandSeed_e[2])
    * 20.0;

  /* Update for RandomNumber: '<S57>/Random Number1' */
  SYSTEM0_DW.NextOutput_hg[0] = rt_nrand_Upu32_Yd_f_pw(&SYSTEM0_DW.RandSeed_i[0])
    * 0.031622776601683791;
  SYSTEM0_DW.NextOutput_hg[1] = rt_nrand_Upu32_Yd_f_pw(&SYSTEM0_DW.RandSeed_i[1])
    * 0.031622776601683791;
  SYSTEM0_DW.NextOutput_hg[2] = rt_nrand_Upu32_Yd_f_pw(&SYSTEM0_DW.RandSeed_i[2])
    * 0.031622776601683791;
}

/* Model step wrapper function for compatibility with a static main program */
void SYSTEM0_step(int_T tid)
{
  switch (tid) {
   case 0 :
    SYSTEM0_step0();
    break;

   case 2 :
    SYSTEM0_step2();
    break;

   case 3 :
    SYSTEM0_step3();
    break;

   default :
    break;
  }
}

/* Model initialize function */
void SYSTEM0_initialize(void)
{
  /* Registration code */

  /* initialize real-time model */
  (void) memset((void *)SYSTEM0_M, 0,
                sizeof(RT_MODEL_SYSTEM0_T));
  (SYSTEM0_M)->Timing.TaskCounters.cLimit[0] = 1;
  (SYSTEM0_M)->Timing.TaskCounters.cLimit[1] = 1;
  (SYSTEM0_M)->Timing.TaskCounters.cLimit[2] = 10;
  (SYSTEM0_M)->Timing.TaskCounters.cLimit[3] = 250;

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&SYSTEM0_M->solverInfo, &SYSTEM0_M->Timing.simTimeStep);
    rtsiSetTPtr(&SYSTEM0_M->solverInfo, &rtmGetTPtr(SYSTEM0_M));
    rtsiSetStepSizePtr(&SYSTEM0_M->solverInfo, &SYSTEM0_M->Timing.stepSize0);
    rtsiSetdXPtr(&SYSTEM0_M->solverInfo, &SYSTEM0_M->ModelData.derivs);
    rtsiSetContStatesPtr(&SYSTEM0_M->solverInfo, (real_T **)
                         &SYSTEM0_M->ModelData.contStates);
    rtsiSetNumContStatesPtr(&SYSTEM0_M->solverInfo,
      &SYSTEM0_M->Sizes.numContStates);
    rtsiSetErrorStatusPtr(&SYSTEM0_M->solverInfo, (&rtmGetErrorStatus(SYSTEM0_M)));
    rtsiSetRTModelPtr(&SYSTEM0_M->solverInfo, SYSTEM0_M);
  }

  rtsiSetSimTimeStep(&SYSTEM0_M->solverInfo, MAJOR_TIME_STEP);
  SYSTEM0_M->ModelData.intgData.y = SYSTEM0_M->ModelData.odeY;
  SYSTEM0_M->ModelData.intgData.f[0] = SYSTEM0_M->ModelData.odeF[0];
  SYSTEM0_M->ModelData.intgData.f[1] = SYSTEM0_M->ModelData.odeF[1];
  SYSTEM0_M->ModelData.intgData.f[2] = SYSTEM0_M->ModelData.odeF[2];
  SYSTEM0_M->ModelData.intgData.f[3] = SYSTEM0_M->ModelData.odeF[3];
  SYSTEM0_M->ModelData.contStates = ((X_SYSTEM0_T *) &SYSTEM0_X);
  rtsiSetSolverData(&SYSTEM0_M->solverInfo, (void *)
                    &SYSTEM0_M->ModelData.intgData);
  rtsiSetSolverName(&SYSTEM0_M->solverInfo,"ode4");
  rtmSetTPtr(SYSTEM0_M, &SYSTEM0_M->Timing.tArray[0]);
  SYSTEM0_M->Timing.stepSize0 = 0.004;
  rtmSetFirstInitCond(SYSTEM0_M, 1);

  /* block I/O */
  (void) memset(((void *) &SYSTEM0_B), 0,
                sizeof(B_SYSTEM0_T));

  /* states (continuous) */
  {
    (void) memset((void *)&SYSTEM0_X, 0,
                  sizeof(X_SYSTEM0_T));
  }

  /* states (dwork) */
  (void) memset((void *)&SYSTEM0_DW, 0,
                sizeof(DW_SYSTEM0_T));

  /* external inputs */
  (void) memset((void *)&SYSTEM0_U, 0,
                sizeof(ExtU_SYSTEM0_T));

  /* external outputs */
  (void) memset((void *)&SYSTEM0_Y, 0,
                sizeof(ExtY_SYSTEM0_T));

  /* Start for RateTransition: '<S3>/Rate Transition1' */
  SYSTEM0_B.RateTransition1 = 0.0;

  /* Start for RateTransition: '<S3>/Rate Transition3' */
  SYSTEM0_B.RateTransition3 = 0.0;

  /* Start for RateTransition: '<S3>/Rate Transition6' */
  SYSTEM0_B.RateTransition6 = 0.0;

  /* Start for RateTransition: '<S3>/Rate Transition5' */
  SYSTEM0_B.RateTransition5 = 0.0;

  /* Start for RateTransition: '<S3>/Rate Transition7' */
  SYSTEM0_B.RateTransition7 = 0.0;

  /* Start for RateTransition: '<S3>/Rate Transition8' */
  SYSTEM0_B.RateTransition8[0] = 0.0;
  SYSTEM0_B.RateTransition8[1] = 0.0;
  SYSTEM0_B.RateTransition8[2] = 0.0;

  /* Start for RateTransition: '<S3>/Rate Transition10' */
  SYSTEM0_B.RateTransition10 = 0.0;

  /* Start for If: '<S55>/If' */
  SYSTEM0_DW.If_ActiveSubsystem = -1;

  /* Start for IfAction SubSystem: '<S55>/Negative Trace' */
  /* Start for If: '<S66>/Find Maximum Diagonal Value' */
  SYSTEM0_DW.FindMaximumDiagonalValue_ActiveSubsystem = -1;

  /* End of Start for SubSystem: '<S55>/Negative Trace' */
  {
    uint32_T y;
    real_T y1;

    /* InitializeConditions for Integrator: '<S10>/q0 q1 q2 q3' */
    if (rtmIsFirstInitCond(SYSTEM0_M)) {
      SYSTEM0_X.q0q1q2q3_CSTATE[0] = 0.0;
      SYSTEM0_X.q0q1q2q3_CSTATE[1] = 0.0;
      SYSTEM0_X.q0q1q2q3_CSTATE[2] = 0.0;
      SYSTEM0_X.q0q1q2q3_CSTATE[3] = 0.0;
    }

    SYSTEM0_DW.q0q1q2q3_IWORK.IcNeedsLoading = 1;

    /* InitializeConditions for Integrator: '<S6>/xe,ye,ze' */
    SYSTEM0_X.xeyeze_CSTATE[0] = 0.0;
    SYSTEM0_X.xeyeze_CSTATE[1] = 0.0;
    SYSTEM0_X.xeyeze_CSTATE[2] = 0.0;

    /* InitializeConditions for Integrator: '<S6>/ub,vb,wb' */
    SYSTEM0_X.ubvbwb_CSTATE[0] = 0.0;
    SYSTEM0_X.ubvbwb_CSTATE[1] = 0.0;
    SYSTEM0_X.ubvbwb_CSTATE[2] = 0.0;

    /* InitializeConditions for Memory: '<S4>/Memory2' */
    SYSTEM0_DW.Memory2_PreviousInput = 0.0;

    /* InitializeConditions for SecondOrderIntegrator: '<S49>/Integrator, Second-Order' */
    SYSTEM0_X.IntegratorSecondOrder_CSTATE[0] = 0.0;
    SYSTEM0_X.IntegratorSecondOrder_CSTATE[1] = 0.0;
    SYSTEM0_DW.IntegratorSecondOrder_MODE = 0;

    /* InitializeConditions for SecondOrderIntegrator: '<S50>/Integrator, Second-Order' */
    SYSTEM0_X.IntegratorSecondOrder_CSTATE_h[0] = 0.0;
    SYSTEM0_X.IntegratorSecondOrder_CSTATE_h[1] = 0.0;
    SYSTEM0_DW.IntegratorSecondOrder_MODE_p = 0;

    /* InitializeConditions for SecondOrderIntegrator: '<S51>/Integrator, Second-Order' */
    SYSTEM0_X.IntegratorSecondOrder_CSTATE_n[0] = 0.0;
    SYSTEM0_X.IntegratorSecondOrder_CSTATE_n[1] = 0.0;
    SYSTEM0_DW.IntegratorSecondOrder_MODE_a = 0;

    /* InitializeConditions for SecondOrderIntegrator: '<S52>/Integrator, Second-Order' */
    SYSTEM0_X.IntegratorSecondOrder_CSTATE_d[0] = 0.0;
    SYSTEM0_X.IntegratorSecondOrder_CSTATE_d[1] = 0.0;
    SYSTEM0_DW.IntegratorSecondOrder_MODE_pu = 0;

    /* InitializeConditions for Integrator: '<S6>/p,q,r ' */
    SYSTEM0_X.pqr_CSTATE[0] = 0.0;
    SYSTEM0_X.pqr_CSTATE[1] = 0.0;
    SYSTEM0_X.pqr_CSTATE[2] = 0.0;

    /* InitializeConditions for RandomNumber: '<S132>/White Noise' */
    y = 1373044741UL;
    y1 = rt_nrand_Upu32_Yd_f_pw(&y);
    SYSTEM0_DW.NextOutput[0] = y1;
    SYSTEM0_DW.RandSeed[0] = y;
    y = 411009029UL;
    y1 = rt_nrand_Upu32_Yd_f_pw(&y);
    SYSTEM0_DW.NextOutput[1] = y1;
    SYSTEM0_DW.RandSeed[1] = y;
    y = 1845526542UL;
    y1 = rt_nrand_Upu32_Yd_f_pw(&y);
    SYSTEM0_DW.NextOutput[2] = y1;
    SYSTEM0_DW.RandSeed[2] = y;

    /* InitializeConditions for RateTransition: '<S3>/Rate Transition1' */
    SYSTEM0_DW.RateTransition1_Buffer0 = 0.0;

    /* InitializeConditions for RateTransition: '<S3>/Rate Transition3' */
    SYSTEM0_DW.RateTransition3_Buffer0 = 0.0;

    /* InitializeConditions for RandomNumber: '<S53>/Random Number' */
    SYSTEM0_DW.RandSeed_f = 1144108930UL;
    SYSTEM0_DW.NextOutput_n = rt_nrand_Upu32_Yd_f_pw(&SYSTEM0_DW.RandSeed_f) *
      0.01;

    /* InitializeConditions for RandomNumber: '<S65>/Random Number' */
    SYSTEM0_DW.RandSeed_d = 1144108930UL;
    SYSTEM0_DW.NextOutput_h = rt_nrand_Upu32_Yd_f_pw(&SYSTEM0_DW.RandSeed_d) *
      0.01;

    /* InitializeConditions for RateTransition: '<S3>/Rate Transition6' */
    SYSTEM0_DW.RateTransition6_Buffer0 = 0.0;

    /* InitializeConditions for RateTransition: '<S3>/Rate Transition5' */
    SYSTEM0_DW.RateTransition5_Buffer0 = 0.0;

    /* InitializeConditions for RateTransition: '<S3>/Rate Transition7' */
    SYSTEM0_DW.RateTransition7_Buffer0 = 0.0;

    /* InitializeConditions for RateTransition: '<S3>/Rate Transition8' */
    SYSTEM0_DW.RateTransition8_Buffer0[0] = 0.0;
    SYSTEM0_DW.RateTransition8_Buffer0[1] = 0.0;
    SYSTEM0_DW.RateTransition8_Buffer0[2] = 0.0;

    /* InitializeConditions for RateTransition: '<S3>/Rate Transition10' */
    SYSTEM0_DW.RateTransition10_Buffer0 = 0.0;

    /* InitializeConditions for RandomNumber: '<S149>/White Noise' */
    y = 1689616386UL;
    y1 = rt_nrand_Upu32_Yd_f_pw(&y);
    SYSTEM0_DW.NextOutput_p[0] = y1;
    SYSTEM0_DW.RandSeed_j[0] = y;
    y = 1998225409UL;
    y1 = rt_nrand_Upu32_Yd_f_pw(&y);
    SYSTEM0_DW.NextOutput_p[1] = y1;
    SYSTEM0_DW.RandSeed_j[1] = y;
    y = 1181220867UL;
    y1 = rt_nrand_Upu32_Yd_f_pw(&y);
    SYSTEM0_DW.NextOutput_p[2] = y1;
    SYSTEM0_DW.RandSeed_j[2] = y;

    /* InitializeConditions for RandomNumber: '<S60>/Random Number' */
    SYSTEM0_DW.RandSeed_p = 1144108930UL;
    SYSTEM0_DW.NextOutput_a = rt_nrand_Upu32_Yd_f_pw(&SYSTEM0_DW.RandSeed_p) *
      0.0031622776601683794;

    /* InitializeConditions for MATLAB Function: '<S4>/LiPo Battery' */
    SYSTEM0_DW.discharge = 0.0;

    /* InitializeConditions for RandomNumber: '<S57>/Random Number' */
    y = 1144108930UL;
    y1 = rt_nrand_Upu32_Yd_f_pw(&y) * 15.0;
    SYSTEM0_DW.NextOutput_o[0] = y1;
    SYSTEM0_DW.RandSeed_e[0] = y;
    y = 1144108930UL;
    y1 = rt_nrand_Upu32_Yd_f_pw(&y) * 15.0;
    SYSTEM0_DW.NextOutput_o[1] = y1;
    SYSTEM0_DW.RandSeed_e[1] = y;
    y = 1144108930UL;
    y1 = rt_nrand_Upu32_Yd_f_pw(&y) * 20.0;
    SYSTEM0_DW.NextOutput_o[2] = y1;
    SYSTEM0_DW.RandSeed_e[2] = y;

    /* InitializeConditions for RandomNumber: '<S57>/Random Number1' */
    y = 1144108930UL;
    y1 = rt_nrand_Upu32_Yd_f_pw(&y) * 0.031622776601683791;
    SYSTEM0_DW.NextOutput_hg[0] = y1;
    SYSTEM0_DW.RandSeed_i[0] = y;
    y = 1144108930UL;
    y1 = rt_nrand_Upu32_Yd_f_pw(&y) * 0.031622776601683791;
    SYSTEM0_DW.NextOutput_hg[1] = y1;
    SYSTEM0_DW.RandSeed_i[1] = y;
    y = 1144108930UL;
    y1 = rt_nrand_Upu32_Yd_f_pw(&y) * 0.031622776601683791;
    SYSTEM0_DW.NextOutput_hg[2] = y1;
    SYSTEM0_DW.RandSeed_i[2] = y;

    /* InitializeConditions for RandomNumber: '<S61>/Random Number' */
    SYSTEM0_DW.RandSeed_fw = 1144108930UL;
    SYSTEM0_DW.NextOutput_l = rt_nrand_Upu32_Yd_f_pw(&SYSTEM0_DW.RandSeed_fw) *
      0.0031622776601683794;

    /* InitializeConditions for RandomNumber: '<S63>/Random Number' */
    SYSTEM0_DW.RandSeed_h = 1144108930UL;
    SYSTEM0_DW.NextOutput_lc = rt_nrand_Upu32_Yd_f_pw(&SYSTEM0_DW.RandSeed_h) *
      0.0031622776601683794;

    /* set "at time zero" to false */
    if (rtmIsFirstInitCond(SYSTEM0_M)) {
      rtmSetFirstInitCond(SYSTEM0_M, 0);
    }
  }
}

/* Model terminate function */
void SYSTEM0_terminate(void)
{
  /* (no terminate code required) */
}
