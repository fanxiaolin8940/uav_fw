/*
 * Dynamics.c
 *
 * Code generation for model "Dynamics".
 *
 * Model version              : 1.306
 * Simulink Coder version : 8.8 (R2015a) 09-Feb-2015
 * C source code generated on : Wed Jul 20 21:36:49 2016
 *
 * Target selection: grt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64
 * Emulation hardware selection:
 *    Differs from embedded hardware (Generic->16-bit Embedded Processor)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "Dynamics.h"
#include "Dynamics_private.h"

/* Block signals (auto storage) */
B_Dynamics_T Dynamics_B;

/* Continuous states */
X_Dynamics_T Dynamics_X;

/* Block states (auto storage) */
DW_Dynamics_T Dynamics_DW;

/* External inputs (root inport signals with auto storage) */
ExtU_Dynamics_T Dynamics_U;

/* External outputs (root outports fed by signals with auto storage) */
ExtY_Dynamics_T Dynamics_Y;

/* Real-time model */
RT_MODEL_Dynamics_T Dynamics_M_;
RT_MODEL_Dynamics_T *const Dynamics_M = &Dynamics_M_;

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
  Dynamics_derivatives();

  /* f1 = f(t + (h/2), y + (h/2)*f0) */
  temp = 0.5 * h;
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (temp*f0[i]);
  }

  rtsiSetT(si, t + temp);
  rtsiSetdX(si, f1);
  Dynamics_step();
  Dynamics_derivatives();

  /* f2 = f(t + (h/2), y + (h/2)*f1) */
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (temp*f1[i]);
  }

  rtsiSetdX(si, f2);
  Dynamics_step();
  Dynamics_derivatives();

  /* f3 = f(t + h, y + h*f2) */
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (h*f2[i]);
  }

  rtsiSetT(si, tnew);
  rtsiSetdX(si, f3);
  Dynamics_step();
  Dynamics_derivatives();

  /* tnew = t + h
     ynew = y + (h/6)*(f0 + 2*f1 + 2*f2 + 2*f3) */
  temp = h / 6.0;
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + temp*(f0[i] + 2.0*f1[i] + 2.0*f2[i] + f3[i]);
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
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
void Dynamics_step(void)
{
  /* local block i/o variables */
  real_T rtb_q0q1q2q3[4];
  real_T rtb_IntegratorSecondOrder_o2_d;
  real_T rtb_IntegratorSecondOrder_o1_f;
  real_T rtb_IntegratorSecondOrder_o2_a;
  real_T rtb_IntegratorSecondOrder_o1_n;
  real_T rtb_VectorConcatenate[9];
  real_T rtb_IntegratorSecondOrder_o2_c;
  real_T rtb_IntegratorSecondOrder_o1_p;
  real_T tmp[3];
  int16_T i;
  real_T rtb_TmpSignalConversionAtSFunctionInport1_idx_0;
  real_T rtb_TmpSignalConversionAtSFunctionInport1_idx_1;
  real_T rtb_TmpSignalConversionAtSFunctionInport1_idx_2;
  real_T rtb_TmpSignalConversionAtSFunctionInport1_idx_3;
  real_T u0;
  real_T u1;
  real_T u0_0;
  if (rtmIsMajorTimeStep(Dynamics_M)) {
    /* set solver stop time */
    if (!(Dynamics_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&Dynamics_M->solverInfo,
                            ((Dynamics_M->Timing.clockTickH0 + 1) *
        Dynamics_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&Dynamics_M->solverInfo,
                            ((Dynamics_M->Timing.clockTick0 + 1) *
        Dynamics_M->Timing.stepSize0 + Dynamics_M->Timing.clockTickH0 *
        Dynamics_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(Dynamics_M)) {
    Dynamics_M->Timing.t[0] = rtsiGetT(&Dynamics_M->solverInfo);
  }

  /* Integrator: '<S6>/q0 q1 q2 q3' */
  if (Dynamics_DW.q0q1q2q3_IWORK.IcNeedsLoading) {
    Dynamics_X.q0q1q2q3_CSTATE[0] = Dynamics_ConstB.q0;
    Dynamics_X.q0q1q2q3_CSTATE[1] = Dynamics_ConstB.q1;
    Dynamics_X.q0q1q2q3_CSTATE[2] = Dynamics_ConstB.q2;
    Dynamics_X.q0q1q2q3_CSTATE[3] = Dynamics_ConstB.q3;
  }

  rtb_q0q1q2q3[0] = Dynamics_X.q0q1q2q3_CSTATE[0];
  rtb_q0q1q2q3[1] = Dynamics_X.q0q1q2q3_CSTATE[1];
  rtb_q0q1q2q3[2] = Dynamics_X.q0q1q2q3_CSTATE[2];
  rtb_q0q1q2q3[3] = Dynamics_X.q0q1q2q3_CSTATE[3];

  /* Sqrt: '<S28>/sqrt' incorporates:
   *  Product: '<S29>/Product'
   *  Product: '<S29>/Product1'
   *  Product: '<S29>/Product2'
   *  Product: '<S29>/Product3'
   *  Sum: '<S29>/Sum'
   */
  rtb_IntegratorSecondOrder_o2_d = sqrt(((rtb_q0q1q2q3[0] * rtb_q0q1q2q3[0] +
    rtb_q0q1q2q3[1] * rtb_q0q1q2q3[1]) + rtb_q0q1q2q3[2] * rtb_q0q1q2q3[2]) +
    rtb_q0q1q2q3[3] * rtb_q0q1q2q3[3]);

  /* Product: '<S27>/Product' */
  rtb_IntegratorSecondOrder_o1_f = rtb_q0q1q2q3[0] /
    rtb_IntegratorSecondOrder_o2_d;

  /* Product: '<S27>/Product1' */
  rtb_IntegratorSecondOrder_o2_a = rtb_q0q1q2q3[1] /
    rtb_IntegratorSecondOrder_o2_d;

  /* Product: '<S27>/Product2' */
  rtb_IntegratorSecondOrder_o1_n = rtb_q0q1q2q3[2] /
    rtb_IntegratorSecondOrder_o2_d;

  /* Product: '<S27>/Product3' */
  rtb_IntegratorSecondOrder_o2_d = rtb_q0q1q2q3[3] /
    rtb_IntegratorSecondOrder_o2_d;

  /* Sum: '<S17>/Sum' incorporates:
   *  Product: '<S17>/Product'
   *  Product: '<S17>/Product1'
   *  Product: '<S17>/Product2'
   *  Product: '<S17>/Product3'
   */
  rtb_VectorConcatenate[0] = ((rtb_IntegratorSecondOrder_o1_f *
    rtb_IntegratorSecondOrder_o1_f + rtb_IntegratorSecondOrder_o2_a *
    rtb_IntegratorSecondOrder_o2_a) - rtb_IntegratorSecondOrder_o1_n *
    rtb_IntegratorSecondOrder_o1_n) - rtb_IntegratorSecondOrder_o2_d *
    rtb_IntegratorSecondOrder_o2_d;

  /* Gain: '<S20>/Gain' incorporates:
   *  Product: '<S20>/Product2'
   *  Product: '<S20>/Product3'
   *  Sum: '<S20>/Sum'
   */
  rtb_VectorConcatenate[1] = (rtb_IntegratorSecondOrder_o2_a *
    rtb_IntegratorSecondOrder_o1_n - rtb_IntegratorSecondOrder_o2_d *
    rtb_IntegratorSecondOrder_o1_f) * 2.0;

  /* Gain: '<S23>/Gain' incorporates:
   *  Product: '<S23>/Product1'
   *  Product: '<S23>/Product2'
   *  Sum: '<S23>/Sum'
   */
  rtb_VectorConcatenate[2] = (rtb_IntegratorSecondOrder_o1_f *
    rtb_IntegratorSecondOrder_o1_n + rtb_IntegratorSecondOrder_o2_a *
    rtb_IntegratorSecondOrder_o2_d) * 2.0;

  /* Gain: '<S18>/Gain' incorporates:
   *  Product: '<S18>/Product2'
   *  Product: '<S18>/Product3'
   *  Sum: '<S18>/Sum'
   */
  rtb_VectorConcatenate[3] = (rtb_IntegratorSecondOrder_o2_d *
    rtb_IntegratorSecondOrder_o1_f + rtb_IntegratorSecondOrder_o2_a *
    rtb_IntegratorSecondOrder_o1_n) * 2.0;

  /* Sum: '<S21>/Sum' incorporates:
   *  Product: '<S21>/Product'
   *  Product: '<S21>/Product1'
   *  Product: '<S21>/Product2'
   *  Product: '<S21>/Product3'
   */
  rtb_VectorConcatenate[4] = ((rtb_IntegratorSecondOrder_o1_f *
    rtb_IntegratorSecondOrder_o1_f - rtb_IntegratorSecondOrder_o2_a *
    rtb_IntegratorSecondOrder_o2_a) + rtb_IntegratorSecondOrder_o1_n *
    rtb_IntegratorSecondOrder_o1_n) - rtb_IntegratorSecondOrder_o2_d *
    rtb_IntegratorSecondOrder_o2_d;

  /* Gain: '<S24>/Gain' incorporates:
   *  Product: '<S24>/Product1'
   *  Product: '<S24>/Product2'
   *  Sum: '<S24>/Sum'
   */
  rtb_VectorConcatenate[5] = (rtb_IntegratorSecondOrder_o1_n *
    rtb_IntegratorSecondOrder_o2_d - rtb_IntegratorSecondOrder_o1_f *
    rtb_IntegratorSecondOrder_o2_a) * 2.0;

  /* Gain: '<S19>/Gain' incorporates:
   *  Product: '<S19>/Product1'
   *  Product: '<S19>/Product2'
   *  Sum: '<S19>/Sum'
   */
  rtb_VectorConcatenate[6] = (rtb_IntegratorSecondOrder_o2_a *
    rtb_IntegratorSecondOrder_o2_d - rtb_IntegratorSecondOrder_o1_f *
    rtb_IntegratorSecondOrder_o1_n) * 2.0;

  /* Gain: '<S22>/Gain' incorporates:
   *  Product: '<S22>/Product1'
   *  Product: '<S22>/Product2'
   *  Sum: '<S22>/Sum'
   */
  rtb_VectorConcatenate[7] = (rtb_IntegratorSecondOrder_o1_f *
    rtb_IntegratorSecondOrder_o2_a + rtb_IntegratorSecondOrder_o1_n *
    rtb_IntegratorSecondOrder_o2_d) * 2.0;

  /* Sum: '<S25>/Sum' incorporates:
   *  Product: '<S25>/Product'
   *  Product: '<S25>/Product1'
   *  Product: '<S25>/Product2'
   *  Product: '<S25>/Product3'
   */
  rtb_VectorConcatenate[8] = ((rtb_IntegratorSecondOrder_o1_f *
    rtb_IntegratorSecondOrder_o1_f - rtb_IntegratorSecondOrder_o2_a *
    rtb_IntegratorSecondOrder_o2_a) - rtb_IntegratorSecondOrder_o1_n *
    rtb_IntegratorSecondOrder_o1_n) + rtb_IntegratorSecondOrder_o2_d *
    rtb_IntegratorSecondOrder_o2_d;

  /* Product: '<S12>/Product' incorporates:
   *  Integrator: '<S2>/ub,vb,wb'
   *  Math: '<S2>/Transpose'
   */
  for (i = 0; i < 3; i++) {
    Dynamics_B.Product[i] = 0.0;
    Dynamics_B.Product[i] += rtb_VectorConcatenate[3 * i] *
      Dynamics_X.ubvbwb_CSTATE[0];
    Dynamics_B.Product[i] += rtb_VectorConcatenate[3 * i + 1] *
      Dynamics_X.ubvbwb_CSTATE[1];
    Dynamics_B.Product[i] += rtb_VectorConcatenate[3 * i + 2] *
      Dynamics_X.ubvbwb_CSTATE[2];
  }

  /* End of Product: '<S12>/Product' */

  /* Sqrt: '<S31>/sqrt' incorporates:
   *  Product: '<S32>/Product'
   *  Product: '<S32>/Product1'
   *  Product: '<S32>/Product2'
   *  Product: '<S32>/Product3'
   *  Sum: '<S32>/Sum'
   */
  rtb_IntegratorSecondOrder_o1_n = sqrt(((rtb_q0q1q2q3[0] * rtb_q0q1q2q3[0] +
    rtb_q0q1q2q3[1] * rtb_q0q1q2q3[1]) + rtb_q0q1q2q3[2] * rtb_q0q1q2q3[2]) +
    rtb_q0q1q2q3[3] * rtb_q0q1q2q3[3]);

  /* Product: '<S30>/Product' */
  rtb_IntegratorSecondOrder_o2_a = rtb_q0q1q2q3[0] /
    rtb_IntegratorSecondOrder_o1_n;

  /* Product: '<S30>/Product1' */
  rtb_IntegratorSecondOrder_o1_f = rtb_q0q1q2q3[1] /
    rtb_IntegratorSecondOrder_o1_n;

  /* Product: '<S30>/Product2' */
  rtb_IntegratorSecondOrder_o2_d = rtb_q0q1q2q3[2] /
    rtb_IntegratorSecondOrder_o1_n;

  /* Product: '<S30>/Product3' */
  rtb_IntegratorSecondOrder_o1_n = rtb_q0q1q2q3[3] /
    rtb_IntegratorSecondOrder_o1_n;

  /* Fcn: '<S14>/fcn1' */
  u0 = (rtb_IntegratorSecondOrder_o1_f * rtb_IntegratorSecondOrder_o2_d +
        rtb_IntegratorSecondOrder_o2_a * rtb_IntegratorSecondOrder_o1_n) * 2.0;

  /* Fcn: '<S14>/fcn2' */
  u1 = ((rtb_IntegratorSecondOrder_o2_a * rtb_IntegratorSecondOrder_o2_a +
         rtb_IntegratorSecondOrder_o1_f * rtb_IntegratorSecondOrder_o1_f) -
        rtb_IntegratorSecondOrder_o2_d * rtb_IntegratorSecondOrder_o2_d) -
    rtb_IntegratorSecondOrder_o1_n * rtb_IntegratorSecondOrder_o1_n;

  /* Fcn: '<S14>/fcn3' */
  u0_0 = (rtb_IntegratorSecondOrder_o1_f * rtb_IntegratorSecondOrder_o1_n -
          rtb_IntegratorSecondOrder_o2_a * rtb_IntegratorSecondOrder_o2_d) *
    -2.0;

  /* Fcn: '<S14>/fcn4' */
  rtb_IntegratorSecondOrder_o1_p = (rtb_IntegratorSecondOrder_o2_d *
    rtb_IntegratorSecondOrder_o1_n + rtb_IntegratorSecondOrder_o2_a *
    rtb_IntegratorSecondOrder_o1_f) * 2.0;

  /* Fcn: '<S14>/fcn5' */
  rtb_IntegratorSecondOrder_o2_a = ((rtb_IntegratorSecondOrder_o2_a *
    rtb_IntegratorSecondOrder_o2_a - rtb_IntegratorSecondOrder_o1_f *
    rtb_IntegratorSecondOrder_o1_f) - rtb_IntegratorSecondOrder_o2_d *
    rtb_IntegratorSecondOrder_o2_d) + rtb_IntegratorSecondOrder_o1_n *
    rtb_IntegratorSecondOrder_o1_n;

  /* Product: '<S33>/Product' incorporates:
   *  Integrator: '<S2>/p,q,r '
   */
  for (i = 0; i < 3; i++) {
    Dynamics_Y.A_b[i] = 0.0;
    Dynamics_Y.A_b[i] += Dynamics_ConstB.Selector[i] * Dynamics_X.pqr_CSTATE[0];
    Dynamics_Y.A_b[i] += Dynamics_ConstB.Selector[i + 3] *
      Dynamics_X.pqr_CSTATE[1];
    Dynamics_Y.A_b[i] += Dynamics_ConstB.Selector[i + 6] *
      Dynamics_X.pqr_CSTATE[2];
  }

  /* End of Product: '<S33>/Product' */

  /* Product: '<S36>/k x i' incorporates:
   *  Integrator: '<S2>/p,q,r '
   */
  rtb_IntegratorSecondOrder_o2_c = Dynamics_X.pqr_CSTATE[2] * Dynamics_Y.A_b[0];

  /* Product: '<S36>/i x j' incorporates:
   *  Integrator: '<S2>/p,q,r '
   */
  rtb_IntegratorSecondOrder_o1_n = Dynamics_X.pqr_CSTATE[0] * Dynamics_Y.A_b[1];

  /* Product: '<S37>/i x k' incorporates:
   *  Integrator: '<S2>/p,q,r '
   */
  rtb_IntegratorSecondOrder_o1_f = Dynamics_X.pqr_CSTATE[0] * Dynamics_Y.A_b[2];

  /* Product: '<S37>/j x i' incorporates:
   *  Integrator: '<S2>/p,q,r '
   */
  rtb_IntegratorSecondOrder_o2_d = Dynamics_X.pqr_CSTATE[1] * Dynamics_Y.A_b[0];

  /* Sum: '<S35>/Sum' incorporates:
   *  Integrator: '<S2>/p,q,r '
   *  Product: '<S36>/j x k'
   *  Product: '<S37>/k x j'
   */
  Dynamics_Y.A_b[0] = Dynamics_X.pqr_CSTATE[1] * Dynamics_Y.A_b[2] -
    Dynamics_X.pqr_CSTATE[2] * Dynamics_Y.A_b[1];
  Dynamics_Y.A_b[1] = rtb_IntegratorSecondOrder_o2_c -
    rtb_IntegratorSecondOrder_o1_f;
  Dynamics_Y.A_b[2] = rtb_IntegratorSecondOrder_o1_n -
    rtb_IntegratorSecondOrder_o2_d;
  if (rtmIsMajorTimeStep(Dynamics_M)) {
    /* Memory: '<S1>/Memory2' */
    Dynamics_B.Memory2 = Dynamics_DW.Memory2_PreviousInput;
  }

  /* Gain: '<S5>/RPM2RADS' incorporates:
   *  Gain: '<S5>/V2RPM'
   *  SecondOrderIntegrator: '<S45>/Integrator, Second-Order'
   *  SecondOrderIntegrator: '<S46>/Integrator, Second-Order'
   *  SecondOrderIntegrator: '<S47>/Integrator, Second-Order'
   *  SecondOrderIntegrator: '<S48>/Integrator, Second-Order'
   */
  Dynamics_Y.Rotor_Speeds[0] = 950.0 * Dynamics_X.IntegratorSecondOrder_CSTATE[0]
    * 0.10471975511965977;
  Dynamics_Y.Rotor_Speeds[1] = 950.0 *
    Dynamics_X.IntegratorSecondOrder_CSTATE_h[0] * 0.10471975511965977;
  Dynamics_Y.Rotor_Speeds[2] = 950.0 *
    Dynamics_X.IntegratorSecondOrder_CSTATE_n[0] * 0.10471975511965977;
  Dynamics_Y.Rotor_Speeds[3] = 950.0 *
    Dynamics_X.IntegratorSecondOrder_CSTATE_d[0] * 0.10471975511965977;

  /* MATLAB Function: '<S4>/multicopter' incorporates:
   *  Constant: '<S4>/h_ref12'
   *  Constant: '<S4>/h_ref8'
   *  Integrator: '<S2>/ub,vb,wb'
   */
  /* MATLAB Function 'Dynamics/Subsystem/multicopter': '<S42>:1' */
  /* ===============================Parameters================================= */
  /*  */
  /*  a = [ax;ay;az];     Vector with the cross-sectional areas */
  /*  M = [M];            Frame Mass */
  /*  l = [l];            Lenght of the Quadcopter arm */
  /*  Kthr = [Kthr];      Coeff. for the computation of thrust */
  /*  Ktrq = [Ktrq];      Coeff. for the computation of the torque */
  /* ================================Constants================================= */
  /* '<S42>:1:16' */
  /*  coefficient of drag                                       */
  /* ==============================Actuator Mixer============================== */
  /*  [x(roll), y(pitch), z(yaw)] */
  /*  QUAD [X] */
  /* '<S42>:1:23' */
  /*  MIX = [  0   1  -1   ;      % QUAD [+] */
  /*          -1   0   1   ; */
  /*           0  -1  -1   ; */
  /*           1   0   1   ]; */
  /* ==================================Forces==================================  */
  /*  We are evaluating in Body frame     */
  /* '<S42>:1:35' */
  /*  gravitational force */
  /* '<S42>:1:36' */
  /*  drag force */
  /* --------------------------------Thrust Model------------------------------ */
  /* '<S42>:1:39' */
  Dynamics_Y.Thursts[0] = Dynamics_Y.Rotor_Speeds[0] * Dynamics_Y.Rotor_Speeds[0];
  Dynamics_Y.Thursts[1] = Dynamics_Y.Rotor_Speeds[1] * Dynamics_Y.Rotor_Speeds[1];
  Dynamics_Y.Thursts[2] = Dynamics_Y.Rotor_Speeds[2] * Dynamics_Y.Rotor_Speeds[2];
  Dynamics_Y.Thursts[3] = Dynamics_Y.Rotor_Speeds[3] * Dynamics_Y.Rotor_Speeds[3];
  rtb_IntegratorSecondOrder_o2_d = 1.2247084269789534E-5 * Dynamics_B.Memory2;
  Dynamics_Y.Thursts[0] *= rtb_IntegratorSecondOrder_o2_d;
  Dynamics_Y.Thursts[1] *= rtb_IntegratorSecondOrder_o2_d;
  Dynamics_Y.Thursts[2] *= rtb_IntegratorSecondOrder_o2_d;
  Dynamics_Y.Thursts[3] *= rtb_IntegratorSecondOrder_o2_d;

  /*  rotor thrust */
  /* -------------------------------------------------------------------------- */
  /* '<S42>:1:42' */
  Dynamics_Y.Forces[0] = -Dynamics_X.ubvbwb_CSTATE[0] * 10.0 *
    Dynamics_B.Memory2 * 0.016813708498984763 + rtb_VectorConcatenate[6] * 9.81 *
    1.2;
  Dynamics_Y.Forces[1] = -Dynamics_X.ubvbwb_CSTATE[1] * 10.0 *
    Dynamics_B.Memory2 * 0.018813708498984762 + rtb_VectorConcatenate[7] * 9.81 *
    1.2;
  Dynamics_Y.Forces[2] = -Dynamics_X.ubvbwb_CSTATE[2] * 10.0 *
    Dynamics_B.Memory2 * 0.18845573684677208 + rtb_VectorConcatenate[8] * 9.81 *
    1.2;

  /* '<S42>:1:43' */
  Dynamics_Y.Forces[2] -= ((Dynamics_Y.Thursts[0] + Dynamics_Y.Thursts[1]) +
    Dynamics_Y.Thursts[2]) + Dynamics_Y.Thursts[3];

  /* Sum: '<S1>/Add6' incorporates:
   *  Constant: '<S4>/h_ref10'
   *  MATLAB Function: '<S4>/multicopter'
   */
  /* ==================================Moments================================= */
  /*  Thrusts contributions to momentum */
  /* '<S42>:1:48' */
  /*  x moment */
  /* '<S42>:1:49' */
  /*  y moment */
  /*  Torques contributions */
  /* momentum_x = sum(abs(MIX(:, 3)) .* rotor_inertia .* rotors) * omega(1);     % x rotor momentum */
  /* momentum_y = sum(abs(MIX(:, 3)) .* rotor_inertia .* rotors) * omega(2);     % y rotor momentum */
  /* momentum_z = sum(MIX(:, 3) .* rotor_inertia .* rotors) * omega(3);          % z rotor momentum */
  /* --------------------------------Torque Model------------------------------ */
  /* '<S42>:1:57' */
  /*  rotor torque */
  /* -------------------------------------------------------------------------- */
  /* '<S42>:1:60' */
  /*  - [momentum_x; momentum_y; momentum_z]; */
  /* ========================================================================== */
  Dynamics_Y.Torques[0] = ((Dynamics_Y.Thursts[0] * 0.2 * 1.4142135623730951 /
    2.0 + -Dynamics_Y.Thursts[1] * 0.2 * 1.4142135623730951 / 2.0) +
    -Dynamics_Y.Thursts[2] * 0.2 * 1.4142135623730951 / 2.0) +
    Dynamics_Y.Thursts[3] * 0.2 * 1.4142135623730951 / 2.0;
  Dynamics_Y.Torques[1] = ((Dynamics_Y.Thursts[0] * 0.2 * 1.4142135623730951 /
    2.0 + Dynamics_Y.Thursts[1] * 0.2 * 1.4142135623730951 / 2.0) +
    -Dynamics_Y.Thursts[2] * 0.2 * 1.4142135623730951 / 2.0) +
    -Dynamics_Y.Thursts[3] * 0.2 * 1.4142135623730951 / 2.0;
  Dynamics_Y.Torques[2] = ((-7.129366502583864E-8 * Dynamics_B.Memory2 *
    (Dynamics_Y.Rotor_Speeds[0] * Dynamics_Y.Rotor_Speeds[0]) +
    7.129366502583864E-8 * Dynamics_B.Memory2 * (Dynamics_Y.Rotor_Speeds[1] *
    Dynamics_Y.Rotor_Speeds[1])) + -7.129366502583864E-8 * Dynamics_B.Memory2 *
    (Dynamics_Y.Rotor_Speeds[2] * Dynamics_Y.Rotor_Speeds[2])) +
    7.129366502583864E-8 * Dynamics_B.Memory2 * (Dynamics_Y.Rotor_Speeds[3] *
    Dynamics_Y.Rotor_Speeds[3]);

  /* Product: '<S7>/Product2' incorporates:
   *  Sum: '<S7>/Sum2'
   */
  tmp[0] = Dynamics_Y.Torques[0] - Dynamics_Y.A_b[0];
  tmp[1] = Dynamics_Y.Torques[1] - Dynamics_Y.A_b[1];
  tmp[2] = Dynamics_Y.Torques[2] - Dynamics_Y.A_b[2];
  rt_mrdivide_U1d1x3_U2d3x3_Yd1x3(tmp, Dynamics_ConstB.Selector2,
    Dynamics_B.Product2);
  if (rtmIsMajorTimeStep(Dynamics_M)) {
    /* Sum: '<S1>/Add8' incorporates:
     *  Gain: '<S1>/Gain3'
     *  Inport: '<Root>/h'
     */
    Dynamics_B.Add8[0] = 0.0;
    Dynamics_B.Add8[1] = 0.0;
    Dynamics_B.Add8[2] = 5.0 * Dynamics_U.h + -11.772;

    /* Gain: '<S1>/Gain4' incorporates:
     *  Inport: '<Root>/h'
     */
    Dynamics_B.Gain4 = -Dynamics_U.h;
  }

  /* Switch: '<S1>/Switch' incorporates:
   *  Product: '<S1>/Matrix Multiply'
   */
  if (Dynamics_B.Gain4 >= 0.0) {
    for (i = 0; i < 3; i++) {
      Dynamics_Y.A_b[i] = 0.0;
      Dynamics_Y.A_b[i] += rtb_VectorConcatenate[i] * Dynamics_B.Add8[0];
      Dynamics_Y.A_b[i] += rtb_VectorConcatenate[i + 3] * Dynamics_B.Add8[1];
      Dynamics_Y.A_b[i] += rtb_VectorConcatenate[i + 6] * Dynamics_B.Add8[2];
    }
  } else {
    Dynamics_Y.A_b[0] = 0.0;
    Dynamics_Y.A_b[1] = 0.0;
    Dynamics_Y.A_b[2] = 0.0;
  }

  /* End of Switch: '<S1>/Switch' */

  /* Product: '<S2>/Product' incorporates:
   *  Constant: '<S8>/Constant'
   *  Sum: '<S1>/Add'
   */
  Dynamics_Y.A_b[0] = (Dynamics_Y.A_b[0] + Dynamics_Y.Forces[0]) / 1.2;
  Dynamics_Y.A_b[1] = (Dynamics_Y.A_b[1] + Dynamics_Y.Forces[1]) / 1.2;
  Dynamics_Y.A_b[2] = (Dynamics_Y.A_b[2] + Dynamics_Y.Forces[2]) / 1.2;

  /* Sum: '<S2>/Sum' incorporates:
   *  Integrator: '<S2>/p,q,r '
   *  Integrator: '<S2>/ub,vb,wb'
   *  Product: '<S38>/i x j'
   *  Product: '<S38>/j x k'
   *  Product: '<S38>/k x i'
   *  Product: '<S39>/i x k'
   *  Product: '<S39>/j x i'
   *  Product: '<S39>/k x j'
   *  Sum: '<S9>/Sum'
   */
  Dynamics_B.Sum[0] = (Dynamics_X.ubvbwb_CSTATE[1] * Dynamics_X.pqr_CSTATE[2] -
                       Dynamics_X.ubvbwb_CSTATE[2] * Dynamics_X.pqr_CSTATE[1]) +
    Dynamics_Y.A_b[0];
  Dynamics_B.Sum[1] = (Dynamics_X.ubvbwb_CSTATE[2] * Dynamics_X.pqr_CSTATE[0] -
                       Dynamics_X.ubvbwb_CSTATE[0] * Dynamics_X.pqr_CSTATE[2]) +
    Dynamics_Y.A_b[1];
  Dynamics_B.Sum[2] = (Dynamics_X.ubvbwb_CSTATE[0] * Dynamics_X.pqr_CSTATE[1] -
                       Dynamics_X.ubvbwb_CSTATE[1] * Dynamics_X.pqr_CSTATE[0]) +
    Dynamics_Y.A_b[2];

  /* Outport: '<Root>/V_e' */
  Dynamics_Y.V_e[0] = Dynamics_B.Product[0];
  Dynamics_Y.V_e[1] = Dynamics_B.Product[1];
  Dynamics_Y.V_e[2] = Dynamics_B.Product[2];

  /* Outport: '<Root>/X_e' incorporates:
   *  Integrator: '<S2>/xe,ye,ze'
   */
  Dynamics_Y.X_e[0] = Dynamics_X.xeyeze_CSTATE[0];
  Dynamics_Y.X_e[1] = Dynamics_X.xeyeze_CSTATE[1];
  Dynamics_Y.X_e[2] = Dynamics_X.xeyeze_CSTATE[2];

  /* Outport: '<Root>/rpy' incorporates:
   *  Trigonometry: '<S14>/Trigonometric Function3'
   */
  Dynamics_Y.rpy[0] = atan2(rtb_IntegratorSecondOrder_o1_p,
    rtb_IntegratorSecondOrder_o2_a);

  /* Trigonometry: '<S14>/trigFcn' */
  if (u0_0 > 1.0) {
    u0_0 = 1.0;
  } else {
    if (u0_0 < -1.0) {
      u0_0 = -1.0;
    }
  }

  /* Outport: '<Root>/rpy' incorporates:
   *  Trigonometry: '<S14>/Trigonometric Function1'
   *  Trigonometry: '<S14>/trigFcn'
   */
  Dynamics_Y.rpy[1] = asin(u0_0);
  Dynamics_Y.rpy[2] = atan2(u0, u1);

  /* Outport: '<Root>/DCM_be' */
  memcpy(&Dynamics_Y.DCM_be[0L], &rtb_VectorConcatenate[0L], 9U * sizeof(real_T));

  /* Outport: '<Root>/V_b' incorporates:
   *  Integrator: '<S2>/ub,vb,wb'
   */
  Dynamics_Y.V_b[0] = Dynamics_X.ubvbwb_CSTATE[0];
  Dynamics_Y.V_b[1] = Dynamics_X.ubvbwb_CSTATE[1];
  Dynamics_Y.V_b[2] = Dynamics_X.ubvbwb_CSTATE[2];

  /* Outport: '<Root>/Omega' incorporates:
   *  Integrator: '<S2>/p,q,r '
   */
  Dynamics_Y.Omega[0] = Dynamics_X.pqr_CSTATE[0];
  Dynamics_Y.Omega[1] = Dynamics_X.pqr_CSTATE[1];
  Dynamics_Y.Omega[2] = Dynamics_X.pqr_CSTATE[2];
  if (rtmIsMajorTimeStep(Dynamics_M)) {
    /* MATLAB Function: '<S1>/LiPo Battery' incorporates:
     *  Inport: '<Root>/pwm1'
     *  Inport: '<Root>/pwm2'
     *  Inport: '<Root>/pwm3'
     *  Inport: '<Root>/pwm4'
     */
    /* MATLAB Function 'Dynamics/LiPo Battery': '<S3>:1' */
    /* ================================Constants================================= */
    /* '<S3>:1:6' */
    /*  maximum current per motor             (A) */
    /* '<S3>:1:7' */
    /*  current of auxillary components       (A)         */
    /*  LiPo current capacity                 (Ah)                                  */
    /*  LiPo cell count */
    /* ========================================================================== */
    /* '<S3>:1:19' */
    /* '<S3>:1:20' */
    Dynamics_DW.discharge += ((((Dynamics_U.pwm1 * Dynamics_U.pwm1 * 5.5 +
      Dynamics_U.pwm2 * Dynamics_U.pwm2 * 5.5) + Dynamics_U.pwm3 *
      Dynamics_U.pwm3 * 5.5) + Dynamics_U.pwm4 * Dynamics_U.pwm4 * 5.5) + 0.1) *
      1.111111111111111E-6;

    /* '<S3>:1:22' */
    if ((0.0 < Dynamics_DW.discharge) && (Dynamics_DW.discharge <= 0.2)) {
      /* '<S3>:1:24' */
      /* '<S3>:1:25' */
      rtb_IntegratorSecondOrder_o2_d = ((Dynamics_DW.discharge *
        Dynamics_DW.discharge * 16.975 + -14.029 * pow(Dynamics_DW.discharge,
        3.0)) - 5.3339 * Dynamics_DW.discharge) + 4.2;
    } else if ((0.2 < Dynamics_DW.discharge) && (Dynamics_DW.discharge < 0.7)) {
      /* '<S3>:1:26' */
      /* '<S3>:1:27' */
      rtb_IntegratorSecondOrder_o2_d = -0.2 * Dynamics_DW.discharge + 3.74;
    } else {
      /* '<S3>:1:29' */
      rtb_IntegratorSecondOrder_o2_d = ((Dynamics_DW.discharge *
        Dynamics_DW.discharge * 89.6 + -48.0 * pow(Dynamics_DW.discharge, 3.0))
        - 55.08 * Dynamics_DW.discharge) + 14.716;
    }

    if (rtb_IntegratorSecondOrder_o2_d < 2.5) {
      /* '<S3>:1:32' */
      /* '<S3>:1:33' */
      rtb_IntegratorSecondOrder_o2_d = 0.0;
    }

    /* '<S3>:1:36' */
    rtb_IntegratorSecondOrder_o2_d *= 2.0;

    /* '<S3>:1:37' */
    rtb_TmpSignalConversionAtSFunctionInport1_idx_0 =
      rtb_IntegratorSecondOrder_o2_d * Dynamics_U.pwm1;
    rtb_TmpSignalConversionAtSFunctionInport1_idx_1 =
      rtb_IntegratorSecondOrder_o2_d * Dynamics_U.pwm2;
    rtb_TmpSignalConversionAtSFunctionInport1_idx_2 =
      rtb_IntegratorSecondOrder_o2_d * Dynamics_U.pwm3;
    rtb_TmpSignalConversionAtSFunctionInport1_idx_3 =
      rtb_IntegratorSecondOrder_o2_d * Dynamics_U.pwm4;

    /* End of MATLAB Function: '<S1>/LiPo Battery' */
    /* ========================================================================== */
  }

  /* Outport: '<Root>/Omega_dot' */
  Dynamics_Y.Omega_dot[0] = Dynamics_B.Product2[0];
  Dynamics_Y.Omega_dot[1] = Dynamics_B.Product2[1];
  Dynamics_Y.Omega_dot[2] = Dynamics_B.Product2[2];

  /* DotProduct: '<S16>/Dot Product' */
  rtb_IntegratorSecondOrder_o2_d = ((rtb_q0q1q2q3[0] * rtb_q0q1q2q3[0] +
    rtb_q0q1q2q3[1] * rtb_q0q1q2q3[1]) + rtb_q0q1q2q3[2] * rtb_q0q1q2q3[2]) +
    rtb_q0q1q2q3[3] * rtb_q0q1q2q3[3];

  /* Fcn: '<S16>/q0dot' incorporates:
   *  Constant: '<S16>/Constant'
   *  DotProduct: '<S16>/Dot Product'
   *  Integrator: '<S2>/p,q,r '
   *  Sum: '<S16>/Sum'
   */
  Dynamics_B.q0dot = ((rtb_q0q1q2q3[1] * Dynamics_X.pqr_CSTATE[0] +
                       rtb_q0q1q2q3[2] * Dynamics_X.pqr_CSTATE[1]) +
                      rtb_q0q1q2q3[3] * Dynamics_X.pqr_CSTATE[2]) * -0.5 + (1.0
    - rtb_IntegratorSecondOrder_o2_d) * rtb_q0q1q2q3[0];

  /* Fcn: '<S16>/q1dot' incorporates:
   *  Constant: '<S16>/Constant'
   *  DotProduct: '<S16>/Dot Product'
   *  Integrator: '<S2>/p,q,r '
   *  Sum: '<S16>/Sum'
   */
  Dynamics_B.q1dot = ((rtb_q0q1q2q3[0] * Dynamics_X.pqr_CSTATE[0] +
                       rtb_q0q1q2q3[2] * Dynamics_X.pqr_CSTATE[2]) -
                      rtb_q0q1q2q3[3] * Dynamics_X.pqr_CSTATE[1]) * 0.5 + (1.0 -
    rtb_IntegratorSecondOrder_o2_d) * rtb_q0q1q2q3[1];

  /* Fcn: '<S16>/q2dot' incorporates:
   *  Constant: '<S16>/Constant'
   *  DotProduct: '<S16>/Dot Product'
   *  Integrator: '<S2>/p,q,r '
   *  Sum: '<S16>/Sum'
   */
  Dynamics_B.q2dot = ((rtb_q0q1q2q3[0] * Dynamics_X.pqr_CSTATE[1] +
                       rtb_q0q1q2q3[3] * Dynamics_X.pqr_CSTATE[0]) -
                      rtb_q0q1q2q3[1] * Dynamics_X.pqr_CSTATE[2]) * 0.5 + (1.0 -
    rtb_IntegratorSecondOrder_o2_d) * rtb_q0q1q2q3[2];

  /* Fcn: '<S16>/q3dot' incorporates:
   *  Constant: '<S16>/Constant'
   *  DotProduct: '<S16>/Dot Product'
   *  Integrator: '<S2>/p,q,r '
   *  Sum: '<S16>/Sum'
   */
  Dynamics_B.q3dot = ((rtb_q0q1q2q3[0] * Dynamics_X.pqr_CSTATE[2] +
                       rtb_q0q1q2q3[1] * Dynamics_X.pqr_CSTATE[1]) -
                      rtb_q0q1q2q3[2] * Dynamics_X.pqr_CSTATE[0]) * 0.5 + (1.0 -
    rtb_IntegratorSecondOrder_o2_d) * rtb_q0q1q2q3[3];
  if (rtmIsMajorTimeStep(Dynamics_M)) {
    /* Product: '<S5>/Product' incorporates:
     *  Inport: '<Root>/pwm1'
     *  Inport: '<Root>/pwm2'
     *  Inport: '<Root>/pwm3'
     *  Inport: '<Root>/pwm4'
     */
    Dynamics_B.Product_o[0] = Dynamics_U.pwm1 *
      rtb_TmpSignalConversionAtSFunctionInport1_idx_0;
    Dynamics_B.Product_o[1] = Dynamics_U.pwm2 *
      rtb_TmpSignalConversionAtSFunctionInport1_idx_1;
    Dynamics_B.Product_o[2] = Dynamics_U.pwm3 *
      rtb_TmpSignalConversionAtSFunctionInport1_idx_2;
    Dynamics_B.Product_o[3] = Dynamics_U.pwm4 *
      rtb_TmpSignalConversionAtSFunctionInport1_idx_3;
  }

  /* Sum: '<S45>/Sum2' incorporates:
   *  Gain: '<S45>/2*zeta * wn'
   *  Gain: '<S45>/wn^2'
   *  SecondOrderIntegrator: '<S45>/Integrator, Second-Order'
   *  Sum: '<S45>/Sum3'
   */
  Dynamics_B.Sum2 = (Dynamics_B.Product_o[0] -
                     Dynamics_X.IntegratorSecondOrder_CSTATE[0]) * 4900.0 -
    140.0 * Dynamics_X.IntegratorSecondOrder_CSTATE[1];

  /* Sum: '<S46>/Sum2' incorporates:
   *  Gain: '<S46>/2*zeta * wn'
   *  Gain: '<S46>/wn^2'
   *  SecondOrderIntegrator: '<S46>/Integrator, Second-Order'
   *  Sum: '<S46>/Sum3'
   */
  Dynamics_B.Sum2_j = (Dynamics_B.Product_o[1] -
                       Dynamics_X.IntegratorSecondOrder_CSTATE_h[0]) * 4900.0 -
    140.0 * Dynamics_X.IntegratorSecondOrder_CSTATE_h[1];

  /* Sum: '<S47>/Sum2' incorporates:
   *  Gain: '<S47>/2*zeta * wn'
   *  Gain: '<S47>/wn^2'
   *  SecondOrderIntegrator: '<S47>/Integrator, Second-Order'
   *  Sum: '<S47>/Sum3'
   */
  Dynamics_B.Sum2_c = (Dynamics_B.Product_o[2] -
                       Dynamics_X.IntegratorSecondOrder_CSTATE_n[0]) * 4900.0 -
    140.0 * Dynamics_X.IntegratorSecondOrder_CSTATE_n[1];

  /* Sum: '<S48>/Sum2' incorporates:
   *  Gain: '<S48>/2*zeta * wn'
   *  Gain: '<S48>/wn^2'
   *  SecondOrderIntegrator: '<S48>/Integrator, Second-Order'
   *  Sum: '<S48>/Sum3'
   */
  Dynamics_B.Sum2_p = (Dynamics_B.Product_o[3] -
                       Dynamics_X.IntegratorSecondOrder_CSTATE_d[0]) * 4900.0 -
    140.0 * Dynamics_X.IntegratorSecondOrder_CSTATE_d[1];
  if (rtmIsMajorTimeStep(Dynamics_M)) {
    /* Update for Integrator: '<S6>/q0 q1 q2 q3' */
    Dynamics_DW.q0q1q2q3_IWORK.IcNeedsLoading = 0;
    if (rtmIsMajorTimeStep(Dynamics_M)) {
      /* Update for Memory: '<S1>/Memory2' incorporates:
       *  Update for Inport: '<Root>/rho'
       */
      Dynamics_DW.Memory2_PreviousInput = Dynamics_U.rho;
    }
  }                                    /* end MajorTimeStep */

  if (rtmIsMajorTimeStep(Dynamics_M)) {
    rt_ertODEUpdateContinuousStates(&Dynamics_M->solverInfo);

    /* Update absolute time for base rate */
    /* The "clockTick0" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick0"
     * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick0 and the high bits
     * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++Dynamics_M->Timing.clockTick0)) {
      ++Dynamics_M->Timing.clockTickH0;
    }

    Dynamics_M->Timing.t[0] = rtsiGetSolverStopTime(&Dynamics_M->solverInfo);

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
      Dynamics_M->Timing.clockTick1++;
      if (!Dynamics_M->Timing.clockTick1) {
        Dynamics_M->Timing.clockTickH1++;
      }
    }
  }                                    /* end MajorTimeStep */
}

/* Derivatives for root system: '<Root>' */
void Dynamics_derivatives(void)
{
  XDot_Dynamics_T *_rtXdot;
  _rtXdot = ((XDot_Dynamics_T *) Dynamics_M->ModelData.derivs);

  /* Derivatives for Integrator: '<S6>/q0 q1 q2 q3' */
  {
    ((XDot_Dynamics_T *) Dynamics_M->ModelData.derivs)->q0q1q2q3_CSTATE[0] =
      Dynamics_B.q0dot;
    ((XDot_Dynamics_T *) Dynamics_M->ModelData.derivs)->q0q1q2q3_CSTATE[1] =
      Dynamics_B.q1dot;
    ((XDot_Dynamics_T *) Dynamics_M->ModelData.derivs)->q0q1q2q3_CSTATE[2] =
      Dynamics_B.q2dot;
    ((XDot_Dynamics_T *) Dynamics_M->ModelData.derivs)->q0q1q2q3_CSTATE[3] =
      Dynamics_B.q3dot;
  }

  /* Derivatives for Integrator: '<S2>/ub,vb,wb' */
  _rtXdot->ubvbwb_CSTATE[0] = Dynamics_B.Sum[0];
  _rtXdot->ubvbwb_CSTATE[1] = Dynamics_B.Sum[1];
  _rtXdot->ubvbwb_CSTATE[2] = Dynamics_B.Sum[2];

  /* Derivatives for Integrator: '<S2>/xe,ye,ze' */
  _rtXdot->xeyeze_CSTATE[0] = Dynamics_B.Product[0];
  _rtXdot->xeyeze_CSTATE[1] = Dynamics_B.Product[1];
  _rtXdot->xeyeze_CSTATE[2] = Dynamics_B.Product[2];

  /* Derivatives for Integrator: '<S2>/p,q,r ' */
  _rtXdot->pqr_CSTATE[0] = Dynamics_B.Product2[0];
  _rtXdot->pqr_CSTATE[1] = Dynamics_B.Product2[1];
  _rtXdot->pqr_CSTATE[2] = Dynamics_B.Product2[2];

  /* Derivatives for SecondOrderIntegrator: '<S45>/Integrator, Second-Order' */
  if (Dynamics_DW.IntegratorSecondOrder_MODE == 0) {
    _rtXdot->IntegratorSecondOrder_CSTATE[0] =
      Dynamics_X.IntegratorSecondOrder_CSTATE[1];
    _rtXdot->IntegratorSecondOrder_CSTATE[1] = Dynamics_B.Sum2;
  }

  /* End of Derivatives for SecondOrderIntegrator: '<S45>/Integrator, Second-Order' */

  /* Derivatives for SecondOrderIntegrator: '<S46>/Integrator, Second-Order' */
  if (Dynamics_DW.IntegratorSecondOrder_MODE_p == 0) {
    _rtXdot->IntegratorSecondOrder_CSTATE_h[0] =
      Dynamics_X.IntegratorSecondOrder_CSTATE_h[1];
    _rtXdot->IntegratorSecondOrder_CSTATE_h[1] = Dynamics_B.Sum2_j;
  }

  /* End of Derivatives for SecondOrderIntegrator: '<S46>/Integrator, Second-Order' */

  /* Derivatives for SecondOrderIntegrator: '<S47>/Integrator, Second-Order' */
  if (Dynamics_DW.IntegratorSecondOrder_MODE_a == 0) {
    _rtXdot->IntegratorSecondOrder_CSTATE_n[0] =
      Dynamics_X.IntegratorSecondOrder_CSTATE_n[1];
    _rtXdot->IntegratorSecondOrder_CSTATE_n[1] = Dynamics_B.Sum2_c;
  }

  /* End of Derivatives for SecondOrderIntegrator: '<S47>/Integrator, Second-Order' */

  /* Derivatives for SecondOrderIntegrator: '<S48>/Integrator, Second-Order' */
  if (Dynamics_DW.IntegratorSecondOrder_MODE_pu == 0) {
    _rtXdot->IntegratorSecondOrder_CSTATE_d[0] =
      Dynamics_X.IntegratorSecondOrder_CSTATE_d[1];
    _rtXdot->IntegratorSecondOrder_CSTATE_d[1] = Dynamics_B.Sum2_p;
  }

  /* End of Derivatives for SecondOrderIntegrator: '<S48>/Integrator, Second-Order' */
}

/* Model initialize function */
void Dynamics_initialize(void)
{
  /* Registration code */

  /* initialize real-time model */
  (void) memset((void *)Dynamics_M, 0,
                sizeof(RT_MODEL_Dynamics_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&Dynamics_M->solverInfo,
                          &Dynamics_M->Timing.simTimeStep);
    rtsiSetTPtr(&Dynamics_M->solverInfo, &rtmGetTPtr(Dynamics_M));
    rtsiSetStepSizePtr(&Dynamics_M->solverInfo, &Dynamics_M->Timing.stepSize0);
    rtsiSetdXPtr(&Dynamics_M->solverInfo, &Dynamics_M->ModelData.derivs);
    rtsiSetContStatesPtr(&Dynamics_M->solverInfo, (real_T **)
                         &Dynamics_M->ModelData.contStates);
    rtsiSetNumContStatesPtr(&Dynamics_M->solverInfo,
      &Dynamics_M->Sizes.numContStates);
    rtsiSetErrorStatusPtr(&Dynamics_M->solverInfo, (&rtmGetErrorStatus
      (Dynamics_M)));
    rtsiSetRTModelPtr(&Dynamics_M->solverInfo, Dynamics_M);
  }

  rtsiSetSimTimeStep(&Dynamics_M->solverInfo, MAJOR_TIME_STEP);
  Dynamics_M->ModelData.intgData.y = Dynamics_M->ModelData.odeY;
  Dynamics_M->ModelData.intgData.f[0] = Dynamics_M->ModelData.odeF[0];
  Dynamics_M->ModelData.intgData.f[1] = Dynamics_M->ModelData.odeF[1];
  Dynamics_M->ModelData.intgData.f[2] = Dynamics_M->ModelData.odeF[2];
  Dynamics_M->ModelData.intgData.f[3] = Dynamics_M->ModelData.odeF[3];
  Dynamics_M->ModelData.contStates = ((X_Dynamics_T *) &Dynamics_X);
  rtsiSetSolverData(&Dynamics_M->solverInfo, (void *)
                    &Dynamics_M->ModelData.intgData);
  rtsiSetSolverName(&Dynamics_M->solverInfo,"ode4");
  rtmSetTPtr(Dynamics_M, &Dynamics_M->Timing.tArray[0]);
  Dynamics_M->Timing.stepSize0 = 0.004;
  rtmSetFirstInitCond(Dynamics_M, 1);

  /* block I/O */
  (void) memset(((void *) &Dynamics_B), 0,
                sizeof(B_Dynamics_T));

  /* states (continuous) */
  {
    (void) memset((void *)&Dynamics_X, 0,
                  sizeof(X_Dynamics_T));
  }

  /* states (dwork) */
  (void) memset((void *)&Dynamics_DW, 0,
                sizeof(DW_Dynamics_T));

  /* external inputs */
  (void) memset((void *)&Dynamics_U, 0,
                sizeof(ExtU_Dynamics_T));

  /* external outputs */
  (void) memset((void *)&Dynamics_Y, 0,
                sizeof(ExtY_Dynamics_T));

  /* InitializeConditions for Integrator: '<S6>/q0 q1 q2 q3' */
  if (rtmIsFirstInitCond(Dynamics_M)) {
    Dynamics_X.q0q1q2q3_CSTATE[0] = 0.0;
    Dynamics_X.q0q1q2q3_CSTATE[1] = 0.0;
    Dynamics_X.q0q1q2q3_CSTATE[2] = 0.0;
    Dynamics_X.q0q1q2q3_CSTATE[3] = 0.0;
  }

  Dynamics_DW.q0q1q2q3_IWORK.IcNeedsLoading = 1;

  /* InitializeConditions for Integrator: '<S2>/ub,vb,wb' */
  Dynamics_X.ubvbwb_CSTATE[0] = 0.0;
  Dynamics_X.ubvbwb_CSTATE[1] = 0.0;
  Dynamics_X.ubvbwb_CSTATE[2] = 0.0;

  /* InitializeConditions for Integrator: '<S2>/xe,ye,ze' */
  Dynamics_X.xeyeze_CSTATE[0] = 0.0;
  Dynamics_X.xeyeze_CSTATE[1] = 0.0;
  Dynamics_X.xeyeze_CSTATE[2] = 0.0;

  /* InitializeConditions for Integrator: '<S2>/p,q,r ' */
  Dynamics_X.pqr_CSTATE[0] = 0.0;
  Dynamics_X.pqr_CSTATE[1] = 0.0;
  Dynamics_X.pqr_CSTATE[2] = 0.0;

  /* InitializeConditions for Memory: '<S1>/Memory2' */
  Dynamics_DW.Memory2_PreviousInput = 0.0;

  /* InitializeConditions for SecondOrderIntegrator: '<S45>/Integrator, Second-Order' */
  Dynamics_X.IntegratorSecondOrder_CSTATE[0] = 0.0;
  Dynamics_X.IntegratorSecondOrder_CSTATE[1] = 0.0;
  Dynamics_DW.IntegratorSecondOrder_MODE = 0;

  /* InitializeConditions for SecondOrderIntegrator: '<S46>/Integrator, Second-Order' */
  Dynamics_X.IntegratorSecondOrder_CSTATE_h[0] = 0.0;
  Dynamics_X.IntegratorSecondOrder_CSTATE_h[1] = 0.0;
  Dynamics_DW.IntegratorSecondOrder_MODE_p = 0;

  /* InitializeConditions for SecondOrderIntegrator: '<S47>/Integrator, Second-Order' */
  Dynamics_X.IntegratorSecondOrder_CSTATE_n[0] = 0.0;
  Dynamics_X.IntegratorSecondOrder_CSTATE_n[1] = 0.0;
  Dynamics_DW.IntegratorSecondOrder_MODE_a = 0;

  /* InitializeConditions for SecondOrderIntegrator: '<S48>/Integrator, Second-Order' */
  Dynamics_X.IntegratorSecondOrder_CSTATE_d[0] = 0.0;
  Dynamics_X.IntegratorSecondOrder_CSTATE_d[1] = 0.0;
  Dynamics_DW.IntegratorSecondOrder_MODE_pu = 0;

  /* InitializeConditions for MATLAB Function: '<S1>/LiPo Battery' */
  Dynamics_DW.discharge = 0.0;

  /* set "at time zero" to false */
  if (rtmIsFirstInitCond(Dynamics_M)) {
    rtmSetFirstInitCond(Dynamics_M, 0);
  }
}

/* Model terminate function */
void Dynamics_terminate(void)
{
  /* (no terminate code required) */
}
