/*
 * Sensors.c
 *
 * Code generation for model "Sensors".
 *
 * Model version              : 1.309
 * Simulink Coder version : 8.8 (R2015a) 09-Feb-2015
 * C source code generated on : Wed Jul 20 22:40:23 2016
 *
 * Target selection: grt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64
 * Emulation hardware selection:
 *    Differs from embedded hardware (Generic->16-bit Embedded Processor)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "Sensors.h"
#include "Sensors_private.h"

/* Block signals (auto storage) */
B_Sensors_T Sensors_B;

/* Block states (auto storage) */
DW_Sensors_T Sensors_DW;

/* External inputs (root inport signals with auto storage) */
ExtU_Sensors_T Sensors_U;

/* External outputs (root outports fed by signals with auto storage) */
ExtY_Sensors_T Sensors_Y;

/* Real-time model */
RT_MODEL_Sensors_T Sensors_M_;
RT_MODEL_Sensors_T *const Sensors_M = &Sensors_M_;
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

/* Model step function for TID0 */
void Sensors_step0(void)               /* Sample time: [0.004s, 0.0s] */
{
  real_T rtb_jxi;
  real_T rtb_Switch_i;
  real_T rtb_Switch_bz;
  real_T rtb_Product;
  real_T rtb_Output;
  real_T rtb_Sum1_f[3];
  real_T rtb_Product_jl;
  real_T tmp[3];
  int16_T i;
  real_T rtb_Sum1_e;
  real_T rtb_Sum_n_idx_0;
  real_T rtb_Sum4_idx_0;
  real_T rtb_Sum4_idx_1;
  real_T rtb_Sum4_idx_2;

  /* Update the flag to indicate when data transfers from
   *  Sample time: [0.004s, 0.0s] to Sample time: [0.04s, 0.0s]  */
  (Sensors_M->Timing.RateInteraction.TID0_1)++;
  if ((Sensors_M->Timing.RateInteraction.TID0_1) > 9) {
    Sensors_M->Timing.RateInteraction.TID0_1 = 0;
  }

  /* Update the flag to indicate when data transfers from
   *  Sample time: [0.004s, 0.0s] to Sample time: [1.0s, 0.0s]  */
  (Sensors_M->Timing.RateInteraction.TID0_2)++;
  if ((Sensors_M->Timing.RateInteraction.TID0_2) > 249) {
    Sensors_M->Timing.RateInteraction.TID0_2 = 0;
  }

  /* Saturate: '<S7>/Limit  altitude  to troposhere' incorporates:
   *  Inport: '<Root>/X_e'
   *  UnaryMinus: '<S5>/Ze2height'
   */
  if (-Sensors_U.X_e[2] > 11000.0) {
    rtb_Product = 11000.0;
  } else if (-Sensors_U.X_e[2] < 0.0) {
    rtb_Product = 0.0;
  } else {
    rtb_Product = -Sensors_U.X_e[2];
  }

  /* Sum: '<S7>/Sum1' incorporates:
   *  Constant: '<S7>/Sea Level  Temperature'
   *  Gain: '<S7>/Lapse Rate'
   *  Saturate: '<S7>/Limit  altitude  to troposhere'
   */
  rtb_jxi = 288.15 - 0.0065 * rtb_Product;

  /* Gain: '<S7>/1//T0' */
  rtb_Switch_i = 0.00347041471455839 * rtb_jxi;

  /* Math: '<S7>/(T//T0)^(g//LR) ' */
  rtb_Switch_bz = pow(rtb_Switch_i, 5.2558756014667134);

  /* Product: '<S7>/Product' */
  rtb_Product = rtb_Switch_bz / rtb_Switch_i;

  /* Sum: '<S7>/Sum' incorporates:
   *  Constant: '<S7>/Altitude of Troposphere'
   *  Inport: '<Root>/X_e'
   *  UnaryMinus: '<S5>/Ze2height'
   */
  rtb_Switch_i = 11000.0 - (-Sensors_U.X_e[2]);

  /* Saturate: '<S7>/Limit  altitude  to Stratosphere' */
  if (rtb_Switch_i > 0.0) {
    rtb_Switch_i = 0.0;
  } else {
    if (rtb_Switch_i < -9000.0) {
      rtb_Switch_i = -9000.0;
    }
  }

  /* Math: '<S7>/Stratosphere Model' incorporates:
   *  Gain: '<S7>/g//R'
   *  Product: '<S7>/Product1'
   *  Saturate: '<S7>/Limit  altitude  to Stratosphere'
   *
   * About '<S7>/Stratosphere Model':
   *  Operator: exp
   */
  rtb_Switch_i = exp(1.0 / rtb_jxi * (0.034163191409533639 * rtb_Switch_i));

  /* Product: '<S7>/Product3' incorporates:
   *  Gain: '<S7>/rho0'
   */
  Sensors_Y.rho = 1.225 * rtb_Product * rtb_Switch_i;

  /* Product: '<S7>/Product2' incorporates:
   *  Gain: '<S7>/P0'
   */
  rtb_Switch_i *= 101325.0 * rtb_Switch_bz;

  /* RateTransition: '<S12>/Rate Transition' */
  if (Sensors_M->Timing.RateInteraction.TID0_1 == 1) {
    Sensors_B.RateTransition = rtb_Switch_i;

    /* RateTransition: '<S10>/Rate Transition' */
    Sensors_B.RateTransition_b = Sensors_Y.rho;
  }

  /* End of RateTransition: '<S12>/Rate Transition' */

  /* RateTransition: '<S6>/Rate Transition1' incorporates:
   *  Inport: '<Root>/X_e'
   */
  if (Sensors_M->Timing.RateInteraction.TID0_2 == 1) {
    Sensors_B.RateTransition1[0] = Sensors_U.X_e[0];
    Sensors_B.RateTransition1[1] = Sensors_U.X_e[1];
    Sensors_B.RateTransition1[2] = Sensors_U.X_e[2];

    /* RateTransition: '<S6>/Rate Transition' incorporates:
     *  Constant: '<S1>/h_ref'
     *  Inport: '<Root>/X_e'
     */
    Sensors_B.RateTransition_e = 0.0;
  }

  /* End of RateTransition: '<S6>/Rate Transition1' */

  /* Product: '<S42>/rad long ' incorporates:
   *  Inport: '<Root>/X_e'
   *  Product: '<S42>/y*cos'
   */
  rtb_Switch_i = 1.5678559428873849E-7 * Sensors_U.X_e[1];

  /* Sum: '<S5>/Sum' incorporates:
   *  Gain: '<S43>/Unit Conversion'
   *  Inport: '<Root>/X_e'
   *  Product: '<S42>/rad lat'
   *  Product: '<S42>/x*cos'
   */
  rtb_Sum_n_idx_0 = Sensors_U.X_e[0] * 1.5784225029068334E-7 *
    57.295779513082323;
  rtb_Product_jl = 57.295779513082323 * rtb_Switch_i + Sensors_ConstB.Switch_b;

  /* Switch: '<S47>/Switch' incorporates:
   *  Abs: '<S47>/Abs'
   *  Bias: '<S47>/Bias'
   *  Bias: '<S47>/Bias1'
   *  Constant: '<S47>/Constant2'
   *  Math: '<S47>/Math Function1'
   */
  if (fabs(rtb_Sum_n_idx_0) > 180.0) {
    rtb_Sum_n_idx_0 = rt_modd(rtb_Sum_n_idx_0 + 180.0, 360.0) + -180.0;
  }

  /* End of Switch: '<S47>/Switch' */

  /* Abs: '<S44>/Abs1' */
  rtb_Product = fabs(rtb_Sum_n_idx_0);

  /* Switch: '<S44>/Switch' incorporates:
   *  Bias: '<S44>/Bias'
   *  Bias: '<S44>/Bias1'
   *  Constant: '<S40>/Constant'
   *  Constant: '<S40>/Constant1'
   *  Constant: '<S46>/Constant'
   *  Gain: '<S44>/Gain'
   *  Product: '<S44>/Divide1'
   *  RelationalOperator: '<S46>/Compare'
   *  Signum: '<S44>/Sign1'
   *  Switch: '<S40>/Switch1'
   */
  if ((rtb_Product > 90.0) > 0) {
    /* Signum: '<S44>/Sign1' */
    if (rtb_Sum_n_idx_0 < 0.0) {
      rtb_Sum_n_idx_0 = -1.0;
    } else {
      if (rtb_Sum_n_idx_0 > 0.0) {
        rtb_Sum_n_idx_0 = 1.0;
      }
    }

    rtb_Switch_bz = (-(rtb_Product + -90.0) + 90.0) * rtb_Sum_n_idx_0;
    rtb_Switch_i = 180.0;
  } else {
    rtb_Switch_bz = rtb_Sum_n_idx_0;
    rtb_Switch_i = 0.0;
  }

  /* End of Switch: '<S44>/Switch' */

  /* Sum: '<S40>/Sum' */
  rtb_Product = rtb_Switch_i + rtb_Product_jl;

  /* Abs: '<S45>/Abs' */
  rtb_Switch_i = fabs(rtb_Product);

  /* Switch: '<S45>/Switch' incorporates:
   *  Bias: '<S45>/Bias'
   *  Bias: '<S45>/Bias1'
   *  Constant: '<S45>/Constant2'
   *  Math: '<S45>/Math Function1'
   */
  if (rtb_Switch_i > 180.0) {
    rtb_Switch_i = rt_modd(rtb_Product + 180.0, 360.0) + -180.0;
  } else {
    rtb_Switch_i = rtb_Product;
  }

  /* End of Switch: '<S45>/Switch' */

  /* Gain: '<S3>/Output' incorporates:
   *  RandomNumber: '<S3>/White Noise'
   */
  rtb_Output = 0.00019364916731037085 * Sensors_DW.NextOutput_lh3;

  /* Sum: '<S89>/Sum' incorporates:
   *  Product: '<S91>/i x j'
   *  Product: '<S91>/j x k'
   *  Product: '<S91>/k x i'
   *  Product: '<S92>/i x k'
   *  Product: '<S92>/j x i'
   *  Product: '<S92>/k x j'
   */
  rtb_Sum1_f[0] = 0.0;
  rtb_Sum1_f[1] = 0.0;
  rtb_Sum1_f[2] = 0.0;

  /* Sum: '<S78>/Sum' incorporates:
   *  Inport: '<Root>/A_b'
   *  Inport: '<Root>/DCM'
   *  Product: '<S1>/Matrix Multiply1'
   *  Product: '<S78>/Product'
   */
  for (i = 0; i < 3; i++) {
    tmp[i] = (Sensors_U.A_b[i] - Sensors_U.DCM[i + 6] * 9.81) + rtb_Sum1_f[i];
  }

  /* End of Sum: '<S78>/Sum' */

  /* Product: '<S78>/Product' incorporates:
   *  Constant: '<S78>/Scale factors & Cross-coupling  errors'
   */
  for (i = 0; i < 3; i++) {
    rtb_Sum1_f[i] = Sensors_ConstP.pooled22[i + 6] * tmp[2] +
      (Sensors_ConstP.pooled22[i + 3] * tmp[1] + Sensors_ConstP.pooled22[i] *
       tmp[0]);
  }

  /* Sum: '<S78>/Sum4' */
  rtb_Sum4_idx_0 = rtb_Sum1_f[0];
  rtb_Sum4_idx_1 = rtb_Sum1_f[1];
  rtb_Sum4_idx_2 = rtb_Sum1_f[2];

  /* Product: '<S79>/Product' incorporates:
   *  Constant: '<S79>/Scale factors & Cross-coupling  errors '
   *  Inport: '<Root>/Omega '
   */
  for (i = 0; i < 3; i++) {
    rtb_Sum1_f[i] = Sensors_ConstP.pooled22[i + 6] * Sensors_U.Omega[2] +
      (Sensors_ConstP.pooled22[i + 3] * Sensors_U.Omega[1] +
       Sensors_ConstP.pooled22[i] * Sensors_U.Omega[0]);
  }

  /* End of Product: '<S79>/Product' */

  /* Sum: '<S79>/Sum4' */
  rtb_Product = rtb_Sum1_f[1];
  rtb_Sum_n_idx_0 = rtb_Sum1_f[2];

  /* Sum: '<S79>/Sum1' incorporates:
   *  Gain: '<S98>/Output'
   *  RandomNumber: '<S98>/White Noise'
   *  Sum: '<S79>/Sum4'
   */
  rtb_Sum1_f[0] += 0.00070710678118654762 * Sensors_DW.NextOutput_p[0];
  rtb_Sum1_f[1] = 0.00070710678118654762 * Sensors_DW.NextOutput_p[1] +
    rtb_Product;
  rtb_Sum1_e = 0.00070710678118654762 * Sensors_DW.NextOutput_p[2] +
    rtb_Sum_n_idx_0;

  /* Sum: '<S17>/Add' incorporates:
   *  Inport: '<Root>/DCM'
   */
  rtb_Product = (Sensors_U.DCM[0] + Sensors_U.DCM[4]) + Sensors_U.DCM[8];

  /* If: '<S4>/If' incorporates:
   *  If: '<S15>/Find Maximum Diagonal Value'
   *  Inport: '<Root>/DCM'
   *  Sum: '<S17>/Add'
   */
  if (rtb_Product > 0.0) {
    /* Outputs for IfAction SubSystem: '<S4>/Positive Trace' incorporates:
     *  ActionPort: '<S16>/Action Port'
     */
    /* Sqrt: '<S16>/sqrt' incorporates:
     *  Constant: '<S16>/Constant'
     *  Sum: '<S16>/Sum'
     */
    rtb_Product = sqrt(rtb_Product + 1.0);

    /* Gain: '<S16>/Gain' */
    Sensors_B.Merge[0] = 0.5 * rtb_Product;

    /* Gain: '<S16>/Gain1' */
    rtb_Product *= 2.0;

    /* Product: '<S16>/Product' incorporates:
     *  Inport: '<Root>/DCM'
     *  Sum: '<S37>/Add'
     *  Sum: '<S38>/Add'
     *  Sum: '<S39>/Add'
     */
    Sensors_B.Merge[1] = (Sensors_U.DCM[7] - Sensors_U.DCM[5]) / rtb_Product;
    Sensors_B.Merge[2] = (Sensors_U.DCM[2] - Sensors_U.DCM[6]) / rtb_Product;
    Sensors_B.Merge[3] = (Sensors_U.DCM[3] - Sensors_U.DCM[1]) / rtb_Product;

    /* End of Outputs for SubSystem: '<S4>/Positive Trace' */
  } else {
    /* Outputs for IfAction SubSystem: '<S4>/Negative Trace' incorporates:
     *  ActionPort: '<S15>/Action Port'
     */
    if ((Sensors_U.DCM[4] > Sensors_U.DCM[0]) && (Sensors_U.DCM[4] >
         Sensors_U.DCM[8])) {
      /* Outputs for IfAction SubSystem: '<S15>/Maximum Value at DCM(2,2)' incorporates:
       *  ActionPort: '<S19>/Action Port'
       */
      /* If: '<S15>/Find Maximum Diagonal Value' incorporates:
       *  Constant: '<S30>/Constant1'
       *  Constant: '<S31>/Constant'
       *  Gain: '<S19>/Gain'
       *  Gain: '<S19>/Gain1'
       *  Gain: '<S19>/Gain3'
       *  Gain: '<S19>/Gain4'
       *  Inport: '<Root>/DCM'
       *  Product: '<S19>/Product'
       *  Product: '<S30>/Product'
       *  Sqrt: '<S19>/sqrt'
       *  Sum: '<S27>/Add'
       *  Sum: '<S28>/Add'
       *  Sum: '<S29>/Add'
       *  Sum: '<S31>/Add'
       *  Switch: '<S30>/Switch'
       */
      rtb_Product_jl = sqrt(((Sensors_U.DCM[4] - Sensors_U.DCM[0]) -
        Sensors_U.DCM[8]) + 1.0);
      Sensors_B.Merge[2] = 0.5 * rtb_Product_jl;
      if (rtb_Product_jl != 0.0) {
        rtb_Sum_n_idx_0 = 0.5;
      } else {
        rtb_Sum_n_idx_0 = 0.0;
        rtb_Product_jl = 1.0;
      }

      rtb_Product_jl = rtb_Sum_n_idx_0 / rtb_Product_jl;
      Sensors_B.Merge[1] = (Sensors_U.DCM[1] + Sensors_U.DCM[3]) *
        rtb_Product_jl;
      Sensors_B.Merge[3] = (Sensors_U.DCM[5] + Sensors_U.DCM[7]) *
        rtb_Product_jl;
      Sensors_B.Merge[0] = (Sensors_U.DCM[2] - Sensors_U.DCM[6]) *
        rtb_Product_jl;

      /* End of Outputs for SubSystem: '<S15>/Maximum Value at DCM(2,2)' */
    } else if (Sensors_U.DCM[8] > Sensors_U.DCM[0]) {
      /* Outputs for IfAction SubSystem: '<S15>/Maximum Value at DCM(3,3)' incorporates:
       *  ActionPort: '<S20>/Action Port'
       */
      /* If: '<S15>/Find Maximum Diagonal Value' incorporates:
       *  Constant: '<S35>/Constant1'
       *  Constant: '<S36>/Constant'
       *  Gain: '<S20>/Gain'
       *  Gain: '<S20>/Gain1'
       *  Gain: '<S20>/Gain2'
       *  Gain: '<S20>/Gain3'
       *  Inport: '<Root>/DCM'
       *  Product: '<S20>/Product'
       *  Product: '<S35>/Product'
       *  Sqrt: '<S20>/sqrt'
       *  Sum: '<S32>/Add'
       *  Sum: '<S33>/Add'
       *  Sum: '<S34>/Add'
       *  Sum: '<S36>/Add'
       *  Switch: '<S35>/Switch'
       */
      rtb_Product_jl = sqrt(((Sensors_U.DCM[8] - Sensors_U.DCM[0]) -
        Sensors_U.DCM[4]) + 1.0);
      Sensors_B.Merge[3] = 0.5 * rtb_Product_jl;
      if (rtb_Product_jl != 0.0) {
        rtb_Sum_n_idx_0 = 0.5;
      } else {
        rtb_Sum_n_idx_0 = 0.0;
        rtb_Product_jl = 1.0;
      }

      rtb_Product_jl = rtb_Sum_n_idx_0 / rtb_Product_jl;
      Sensors_B.Merge[1] = (Sensors_U.DCM[2] + Sensors_U.DCM[6]) *
        rtb_Product_jl;
      Sensors_B.Merge[2] = (Sensors_U.DCM[5] + Sensors_U.DCM[7]) *
        rtb_Product_jl;
      Sensors_B.Merge[0] = (Sensors_U.DCM[3] - Sensors_U.DCM[1]) *
        rtb_Product_jl;

      /* End of Outputs for SubSystem: '<S15>/Maximum Value at DCM(3,3)' */
    } else {
      /* Outputs for IfAction SubSystem: '<S15>/Maximum Value at DCM(1,1)' incorporates:
       *  ActionPort: '<S18>/Action Port'
       */
      /* If: '<S15>/Find Maximum Diagonal Value' incorporates:
       *  Constant: '<S25>/Constant1'
       *  Constant: '<S26>/Constant'
       *  Gain: '<S18>/Gain'
       *  Gain: '<S18>/Gain1'
       *  Gain: '<S18>/Gain2'
       *  Gain: '<S18>/Gain3'
       *  Inport: '<Root>/DCM'
       *  Product: '<S18>/Product'
       *  Product: '<S25>/Product'
       *  Sqrt: '<S18>/sqrt'
       *  Sum: '<S22>/Add'
       *  Sum: '<S23>/Add'
       *  Sum: '<S24>/Add'
       *  Sum: '<S26>/Add'
       *  Switch: '<S25>/Switch'
       */
      rtb_Product_jl = sqrt(((Sensors_U.DCM[0] - Sensors_U.DCM[4]) -
        Sensors_U.DCM[8]) + 1.0);
      Sensors_B.Merge[1] = 0.5 * rtb_Product_jl;
      if (rtb_Product_jl != 0.0) {
        rtb_Sum_n_idx_0 = 0.5;
      } else {
        rtb_Sum_n_idx_0 = 0.0;
        rtb_Product_jl = 1.0;
      }

      rtb_Product_jl = rtb_Sum_n_idx_0 / rtb_Product_jl;
      Sensors_B.Merge[2] = (Sensors_U.DCM[1] + Sensors_U.DCM[3]) *
        rtb_Product_jl;
      Sensors_B.Merge[3] = (Sensors_U.DCM[2] + Sensors_U.DCM[6]) *
        rtb_Product_jl;
      Sensors_B.Merge[0] = (Sensors_U.DCM[7] - Sensors_U.DCM[5]) *
        rtb_Product_jl;

      /* End of Outputs for SubSystem: '<S15>/Maximum Value at DCM(1,1)' */
    }

    /* End of Outputs for SubSystem: '<S4>/Negative Trace' */
  }

  /* End of If: '<S4>/If' */

  /* Sum: '<S14>/Sum2' incorporates:
   *  Constant: '<S14>/K2C'
   *  RandomNumber: '<S14>/Random Number'
   *  Sum: '<S14>/Add1'
   */
  Sensors_Y.Temp = (273.15 + rtb_jxi) + Sensors_DW.NextOutput;

  /* Saturate: '<S14>/Saturation' */
  if (Sensors_Y.Temp > 85.0) {
    /* Sum: '<S14>/Sum2' incorporates:
     *  Outport: '<Root>/Temp'
     */
    Sensors_Y.Temp = 85.0;
  } else {
    if (Sensors_Y.Temp < 40.0) {
      /* Sum: '<S14>/Sum2' incorporates:
       *  Outport: '<Root>/Temp'
       */
      Sensors_Y.Temp = 40.0;
    }
  }

  /* End of Saturate: '<S14>/Saturation' */

  /* Sum: '<S2>/Add1' incorporates:
   *  Inport: '<Root>/X_e'
   *  RandomNumber: '<S2>/Random Number'
   *  UnaryMinus: '<S5>/Ze2height'
   */
  Sensors_Y.Baro_alt = -Sensors_U.X_e[2] + Sensors_DW.NextOutput_n;

  /* Saturate: '<S2>/Saturation1' */
  if (!(Sensors_Y.Baro_alt >= 0.0)) {
    /* Sum: '<S2>/Add1' incorporates:
     *  Outport: '<Root>/Baro_alt'
     */
    Sensors_Y.Baro_alt = 0.0;
  }

  /* End of Saturate: '<S2>/Saturation1' */

  /* Outport: '<Root>/Real_Lat_Lon' */
  Sensors_Y.Real_Lat_Lon[0] = rtb_Switch_bz;
  Sensors_Y.Real_Lat_Lon[1] = rtb_Switch_i;

  /* Outport: '<Root>/Real_Alt' incorporates:
   *  Inport: '<Root>/X_e'
   *  UnaryMinus: '<S5>/Ze2height'
   */
  Sensors_Y.Real_Alt = -Sensors_U.X_e[2];

  /* Sum: '<S9>/Sum2' incorporates:
   *  Inport: '<Root>/DCM'
   *  Product: '<S9>/Matrix Multiply2'
   *  RandomNumber: '<S9>/Random Number'
   *  Saturate: '<S9>/Saturation'
   */
  for (i = 0; i < 3; i++) {
    tmp[i] = (Sensors_U.DCM[i + 3] * 0.49999999999999994 + Sensors_U.DCM[i] *
              0.86602540378443871) + Sensors_DW.NextOutput_a;
  }

  /* End of Sum: '<S9>/Sum2' */

  /* Saturate: '<S9>/Saturation' incorporates:
   *  Product: '<S9>/Matrix Multiply2'
   */
  if (tmp[0] > 2.0) {
    /* Outport: '<Root>/Magn_meas' */
    Sensors_Y.Magn_meas[0] = 2.0;
  } else if (tmp[0] < -2.0) {
    /* Outport: '<Root>/Magn_meas' */
    Sensors_Y.Magn_meas[0] = -2.0;
  } else {
    /* Outport: '<Root>/Magn_meas' */
    Sensors_Y.Magn_meas[0] = tmp[0];
  }

  if (tmp[1] > 2.0) {
    /* Outport: '<Root>/Magn_meas' */
    Sensors_Y.Magn_meas[1] = 2.0;
  } else if (tmp[1] < -2.0) {
    /* Outport: '<Root>/Magn_meas' */
    Sensors_Y.Magn_meas[1] = -2.0;
  } else {
    /* Outport: '<Root>/Magn_meas' */
    Sensors_Y.Magn_meas[1] = tmp[1];
  }

  if (tmp[2] > 2.0) {
    /* Outport: '<Root>/Magn_meas' */
    Sensors_Y.Magn_meas[2] = 2.0;
  } else if (tmp[2] < -2.0) {
    /* Outport: '<Root>/Magn_meas' */
    Sensors_Y.Magn_meas[2] = -2.0;
  } else {
    /* Outport: '<Root>/Magn_meas' */
    Sensors_Y.Magn_meas[2] = tmp[2];
  }

  /* Outport: '<Root>/RPY' incorporates:
   *  Inport: '<Root>/rpy'
   *  Sum: '<S1>/Add5'
   */
  Sensors_Y.RPY[0] = Sensors_U.rpy[0] + rtb_Output;
  Sensors_Y.RPY[1] = Sensors_U.rpy[1] + rtb_Output;
  Sensors_Y.RPY[2] = Sensors_U.rpy[2] + rtb_Output;

  /* Sum: '<S78>/Sum1' incorporates:
   *  Gain: '<S81>/Output'
   *  RandomNumber: '<S81>/White Noise'
   */
  rtb_Product = 0.011180339887498949 * Sensors_DW.NextOutput_k[0] +
    rtb_Sum4_idx_0;

  /* Saturate: '<S78>/Saturation' */
  if (rtb_Product > 19.62) {
    /* Outport: '<Root>/A_meas' */
    Sensors_Y.A_meas[0] = 19.62;
  } else if (rtb_Product < -19.62) {
    /* Outport: '<Root>/A_meas' */
    Sensors_Y.A_meas[0] = -19.62;
  } else {
    /* Outport: '<Root>/A_meas' */
    Sensors_Y.A_meas[0] = rtb_Product;
  }

  /* Sum: '<S78>/Sum1' incorporates:
   *  Gain: '<S81>/Output'
   *  RandomNumber: '<S81>/White Noise'
   */
  rtb_Product = 0.011180339887498949 * Sensors_DW.NextOutput_k[1] +
    rtb_Sum4_idx_1;

  /* Saturate: '<S78>/Saturation' */
  if (rtb_Product > 19.62) {
    /* Outport: '<Root>/A_meas' */
    Sensors_Y.A_meas[1] = 19.62;
  } else if (rtb_Product < -19.62) {
    /* Outport: '<Root>/A_meas' */
    Sensors_Y.A_meas[1] = -19.62;
  } else {
    /* Outport: '<Root>/A_meas' */
    Sensors_Y.A_meas[1] = rtb_Product;
  }

  /* Sum: '<S78>/Sum1' incorporates:
   *  Gain: '<S81>/Output'
   *  RandomNumber: '<S81>/White Noise'
   */
  rtb_Product = 0.011180339887498949 * Sensors_DW.NextOutput_k[2] +
    rtb_Sum4_idx_2;

  /* Saturate: '<S78>/Saturation' */
  if (rtb_Product > 19.62) {
    /* Outport: '<Root>/A_meas' */
    Sensors_Y.A_meas[2] = 19.62;
  } else if (rtb_Product < -19.62) {
    /* Outport: '<Root>/A_meas' */
    Sensors_Y.A_meas[2] = -19.62;
  } else {
    /* Outport: '<Root>/A_meas' */
    Sensors_Y.A_meas[2] = rtb_Product;
  }

  /* Saturate: '<S79>/Saturation' */
  if (rtb_Sum1_f[0] > 4.36) {
    /* Outport: '<Root>/Omega_meas' */
    Sensors_Y.Omega_meas[0] = 4.36;
  } else if (rtb_Sum1_f[0] < -4.36) {
    /* Outport: '<Root>/Omega_meas' */
    Sensors_Y.Omega_meas[0] = -4.36;
  } else {
    /* Outport: '<Root>/Omega_meas' */
    Sensors_Y.Omega_meas[0] = rtb_Sum1_f[0];
  }

  if (rtb_Sum1_f[1] > 4.36) {
    /* Outport: '<Root>/Omega_meas' */
    Sensors_Y.Omega_meas[1] = 4.36;
  } else if (rtb_Sum1_f[1] < -4.36) {
    /* Outport: '<Root>/Omega_meas' */
    Sensors_Y.Omega_meas[1] = -4.36;
  } else {
    /* Outport: '<Root>/Omega_meas' */
    Sensors_Y.Omega_meas[1] = rtb_Sum1_f[1];
  }

  if (rtb_Sum1_e > 4.36) {
    /* Outport: '<Root>/Omega_meas' */
    Sensors_Y.Omega_meas[2] = 4.36;
  } else if (rtb_Sum1_e < -4.36) {
    /* Outport: '<Root>/Omega_meas' */
    Sensors_Y.Omega_meas[2] = -4.36;
  } else {
    /* Outport: '<Root>/Omega_meas' */
    Sensors_Y.Omega_meas[2] = rtb_Sum1_e;
  }

  /* End of Saturate: '<S79>/Saturation' */

  /* Outport: '<Root>/q' */
  Sensors_Y.q[0] = Sensors_B.Merge[0];
  Sensors_Y.q[1] = Sensors_B.Merge[1];
  Sensors_Y.q[2] = Sensors_B.Merge[2];
  Sensors_Y.q[3] = Sensors_B.Merge[3];

  /* Update for RandomNumber: '<S14>/Random Number' */
  Sensors_DW.NextOutput = rt_nrand_Upu32_Yd_f_pw(&Sensors_DW.RandSeed) * 0.01;

  /* Update for RandomNumber: '<S2>/Random Number' */
  Sensors_DW.NextOutput_n = rt_nrand_Upu32_Yd_f_pw(&Sensors_DW.RandSeed_fm) *
    0.01;

  /* Update for RandomNumber: '<S9>/Random Number' */
  Sensors_DW.NextOutput_a = rt_nrand_Upu32_Yd_f_pw(&Sensors_DW.RandSeed_p) *
    0.0031622776601683794;

  /* Update for RandomNumber: '<S3>/White Noise' */
  Sensors_DW.NextOutput_lh3 = rt_nrand_Upu32_Yd_f_pw(&Sensors_DW.RandSeed_l);

  /* Update for RandomNumber: '<S81>/White Noise' */
  Sensors_DW.NextOutput_k[0] = rt_nrand_Upu32_Yd_f_pw(&Sensors_DW.RandSeed_ls[0]);
  Sensors_DW.NextOutput_k[1] = rt_nrand_Upu32_Yd_f_pw(&Sensors_DW.RandSeed_ls[1]);
  Sensors_DW.NextOutput_k[2] = rt_nrand_Upu32_Yd_f_pw(&Sensors_DW.RandSeed_ls[2]);

  /* Update for RandomNumber: '<S98>/White Noise' */
  Sensors_DW.NextOutput_p[0] = rt_nrand_Upu32_Yd_f_pw(&Sensors_DW.RandSeed_j[0]);
  Sensors_DW.NextOutput_p[1] = rt_nrand_Upu32_Yd_f_pw(&Sensors_DW.RandSeed_j[1]);
  Sensors_DW.NextOutput_p[2] = rt_nrand_Upu32_Yd_f_pw(&Sensors_DW.RandSeed_j[2]);
}

/* Model step function for TID1 */
void Sensors_step1(void)               /* Sample time: [0.04s, 0.0s] */
{
  /* Sum: '<S12>/Sum2' incorporates:
   *  Gain: '<S12>/Bar2mBar'
   *  Gain: '<S12>/Pa2Bar'
   *  RandomNumber: '<S12>/Random Number'
   */
  Sensors_Y.Press = 1.0E-5 * Sensors_B.RateTransition * 1000.0 +
    Sensors_DW.NextOutput_l;

  /* Saturate: '<S12>/Saturation' */
  if (Sensors_Y.Press > 1200.0) {
    /* Sum: '<S12>/Sum2' incorporates:
     *  Outport: '<Root>/Press'
     */
    Sensors_Y.Press = 1200.0;
  } else {
    if (Sensors_Y.Press < 10.0) {
      /* Sum: '<S12>/Sum2' incorporates:
       *  Outport: '<Root>/Press'
       */
      Sensors_Y.Press = 10.0;
    }
  }

  /* End of Saturate: '<S12>/Saturation' */

  /* Sum: '<S10>/Sum2' incorporates:
   *  DotProduct: '<S10>/Dot Product'
   *  Gain: '<S10>/Bar2mBar'
   *  Gain: '<S10>/Gain'
   *  Gain: '<S10>/Pa2Bar'
   *  Inport: '<Root>/V_b'
   *  Product: '<S10>/Product'
   *  RandomNumber: '<S10>/Random Number'
   */
  Sensors_Y.diff_Press = ((Sensors_U.V_b[0] * Sensors_U.V_b[0] + Sensors_U.V_b[1]
    * Sensors_U.V_b[1]) + Sensors_U.V_b[2] * Sensors_U.V_b[2]) *
    Sensors_B.RateTransition_b * 0.5 * 1.0E-5 * 1000.0 +
    Sensors_DW.NextOutput_lh;

  /* Saturate: '<S10>/Saturation' */
  if (Sensors_Y.diff_Press > 1000.0) {
    /* Sum: '<S10>/Sum2' incorporates:
     *  Outport: '<Root>/diff_Press'
     */
    Sensors_Y.diff_Press = 1000.0;
  } else {
    if (Sensors_Y.diff_Press < 0.0) {
      /* Sum: '<S10>/Sum2' incorporates:
       *  Outport: '<Root>/diff_Press'
       */
      Sensors_Y.diff_Press = 0.0;
    }
  }

  /* End of Saturate: '<S10>/Saturation' */

  /* Update for RandomNumber: '<S12>/Random Number' */
  Sensors_DW.NextOutput_l = rt_nrand_Upu32_Yd_f_pw(&Sensors_DW.RandSeed_h) *
    0.0031622776601683794;

  /* Update for RandomNumber: '<S10>/Random Number' */
  Sensors_DW.NextOutput_lh = rt_nrand_Upu32_Yd_f_pw(&Sensors_DW.RandSeed_f) *
    0.0031622776601683794;
}

/* Model step function for TID2 */
void Sensors_step2(void)               /* Sample time: [1.0s, 0.0s] */
{
  real_T rtb_Abs1;
  int16_T rtb_Compare_0;

  /* RandomNumber: '<S6>/Random Number' */
  Sensors_Y.V_meas[0] = Sensors_DW.NextOutput_o[0];
  Sensors_Y.V_meas[1] = Sensors_DW.NextOutput_o[1];
  Sensors_Y.V_meas[2] = Sensors_DW.NextOutput_o[2];

  /* Sum: '<S6>/Add1' */
  Sensors_Y.V_meas[0] += Sensors_B.RateTransition1[0];
  Sensors_Y.V_meas[1] += Sensors_B.RateTransition1[1];
  Sensors_Y.V_meas[2] += Sensors_B.RateTransition1[2];

  /* Sum: '<S58>/Sum' incorporates:
   *  Gain: '<S62>/Unit Conversion'
   *  Product: '<S61>/rad lat'
   *  Product: '<S61>/x*cos'
   */
  Sensors_Y.Lat_meas = Sensors_Y.V_meas[0] * 1.5784225029068334E-7 *
    57.295779513082323;

  /* Switch: '<S66>/Switch' incorporates:
   *  Abs: '<S66>/Abs'
   */
  if (fabs(Sensors_Y.Lat_meas) > 180.0) {
    /* Sum: '<S58>/Sum' incorporates:
     *  Bias: '<S66>/Bias'
     *  Bias: '<S66>/Bias1'
     *  Constant: '<S66>/Constant2'
     *  Math: '<S66>/Math Function1'
     */
    Sensors_Y.Lat_meas = rt_modd(Sensors_Y.Lat_meas + 180.0, 360.0) + -180.0;
  }

  /* End of Switch: '<S66>/Switch' */

  /* Abs: '<S63>/Abs1' */
  rtb_Abs1 = fabs(Sensors_Y.Lat_meas);

  /* Switch: '<S63>/Switch' incorporates:
   *  Constant: '<S59>/Constant'
   *  Constant: '<S59>/Constant1'
   *  Constant: '<S65>/Constant'
   *  RelationalOperator: '<S65>/Compare'
   *  Switch: '<S59>/Switch1'
   */
  if ((rtb_Abs1 > 90.0) > 0) {
    /* Signum: '<S63>/Sign1' */
    if (Sensors_Y.Lat_meas < 0.0) {
      /* Sum: '<S58>/Sum' */
      Sensors_Y.Lat_meas = -1.0;
    } else {
      if (Sensors_Y.Lat_meas > 0.0) {
        /* Sum: '<S58>/Sum' */
        Sensors_Y.Lat_meas = 1.0;
      }
    }

    /* Sum: '<S58>/Sum' incorporates:
     *  Bias: '<S63>/Bias'
     *  Bias: '<S63>/Bias1'
     *  Gain: '<S63>/Gain'
     *  Outport: '<Root>/Lat_meas'
     *  Product: '<S63>/Divide1'
     *  Signum: '<S63>/Sign1'
     */
    Sensors_Y.Lat_meas *= -(rtb_Abs1 + -90.0) + 90.0;
    rtb_Compare_0 = 180;
  } else {
    rtb_Compare_0 = 0;
  }

  /* End of Switch: '<S63>/Switch' */

  /* Sum: '<S59>/Sum' incorporates:
   *  Gain: '<S62>/Unit Conversion'
   *  Product: '<S61>/rad long '
   *  Product: '<S61>/y*cos'
   *  Sum: '<S58>/Sum'
   */
  Sensors_Y.Lon_meas = (1.5678559428873849E-7 * Sensors_Y.V_meas[1] *
                        57.295779513082323 + Sensors_ConstB.Switch_d) + (real_T)
    rtb_Compare_0;

  /* Switch: '<S64>/Switch' incorporates:
   *  Abs: '<S64>/Abs'
   */
  if (fabs(Sensors_Y.Lon_meas) > 180.0) {
    /* Outport: '<Root>/Lon_meas' incorporates:
     *  Bias: '<S64>/Bias'
     *  Bias: '<S64>/Bias1'
     *  Constant: '<S64>/Constant2'
     *  Math: '<S64>/Math Function1'
     */
    Sensors_Y.Lon_meas = rt_modd(Sensors_Y.Lon_meas + 180.0, 360.0) + -180.0;
  }

  /* End of Switch: '<S64>/Switch' */

  /* UnaryMinus: '<S58>/Ze2height' */
  rtb_Abs1 = -Sensors_Y.V_meas[2];

  /* RandomNumber: '<S6>/Random Number1' */
  Sensors_Y.V_meas[0] = Sensors_DW.NextOutput_h[0];
  Sensors_Y.V_meas[1] = Sensors_DW.NextOutput_h[1];
  Sensors_Y.V_meas[2] = Sensors_DW.NextOutput_h[2];

  /* Sum: '<S6>/Add2' incorporates:
   *  Inport: '<Root>/V_e'
   */
  Sensors_Y.V_meas[0] += Sensors_U.V_e[0];
  Sensors_Y.V_meas[1] += Sensors_U.V_e[1];
  Sensors_Y.V_meas[2] += Sensors_U.V_e[2];

  /* Outport: '<Root>/Alt_meas' incorporates:
   *  Sum: '<S58>/Sum1'
   */
  Sensors_Y.Alt_meas = rtb_Abs1 - Sensors_B.RateTransition_e;

  /* Outport: '<Root>/COG' incorporates:
   *  Trigonometry: '<S6>/Trigonometric Function'
   */
  Sensors_Y.COG = atan2(Sensors_Y.V_meas[1], Sensors_Y.V_meas[0]);

  /* Update for RandomNumber: '<S6>/Random Number' */
  Sensors_DW.NextOutput_o[0] = rt_nrand_Upu32_Yd_f_pw(&Sensors_DW.RandSeed_e[0])
    * 15.0;
  Sensors_DW.NextOutput_o[1] = rt_nrand_Upu32_Yd_f_pw(&Sensors_DW.RandSeed_e[1])
    * 15.0;
  Sensors_DW.NextOutput_o[2] = rt_nrand_Upu32_Yd_f_pw(&Sensors_DW.RandSeed_e[2])
    * 20.0;

  /* Update for RandomNumber: '<S6>/Random Number1' */
  Sensors_DW.NextOutput_h[0] = rt_nrand_Upu32_Yd_f_pw(&Sensors_DW.RandSeed_i[0])
    * 0.031622776601683791;
  Sensors_DW.NextOutput_h[1] = rt_nrand_Upu32_Yd_f_pw(&Sensors_DW.RandSeed_i[1])
    * 0.031622776601683791;
  Sensors_DW.NextOutput_h[2] = rt_nrand_Upu32_Yd_f_pw(&Sensors_DW.RandSeed_i[2])
    * 0.031622776601683791;
}

/* Model step wrapper function for compatibility with a static main program */
void Sensors_step(int_T tid)
{
  switch (tid) {
   case 0 :
    Sensors_step0();
    break;

   case 1 :
    Sensors_step1();
    break;

   case 2 :
    Sensors_step2();
    break;

   default :
    break;
  }
}

/* Model initialize function */
void Sensors_initialize(void)
{
  /* Registration code */

  /* initialize real-time model */
  (void) memset((void *)Sensors_M, 0,
                sizeof(RT_MODEL_Sensors_T));
  (Sensors_M)->Timing.TaskCounters.cLimit[0] = 1;
  (Sensors_M)->Timing.TaskCounters.cLimit[1] = 10;
  (Sensors_M)->Timing.TaskCounters.cLimit[2] = 250;

  /* block I/O */
  (void) memset(((void *) &Sensors_B), 0,
                sizeof(B_Sensors_T));

  /* states (dwork) */
  (void) memset((void *)&Sensors_DW, 0,
                sizeof(DW_Sensors_T));

  /* external inputs */
  (void) memset((void *)&Sensors_U, 0,
                sizeof(ExtU_Sensors_T));

  /* external outputs */
  (void) memset((void *)&Sensors_Y, 0,
                sizeof(ExtY_Sensors_T));

  /* Start for Merge: '<S4>/Merge' */
  Sensors_B.Merge[0] = 1.0;
  Sensors_B.Merge[1] = 0.0;
  Sensors_B.Merge[2] = 0.0;
  Sensors_B.Merge[3] = 0.0;

  /* ConstCode for Outport: '<Root>/Sonar_alt' incorporates:
   *  Constant: '<S1>/Constant2'
   */
  Sensors_Y.Sonar_alt = 0.0;

  {
    uint32_T y;
    real_T y1;

    /* InitializeConditions for RandomNumber: '<S14>/Random Number' */
    Sensors_DW.RandSeed = 1144108930UL;
    Sensors_DW.NextOutput = rt_nrand_Upu32_Yd_f_pw(&Sensors_DW.RandSeed) * 0.01;

    /* InitializeConditions for RandomNumber: '<S12>/Random Number' */
    Sensors_DW.RandSeed_h = 1144108930UL;
    Sensors_DW.NextOutput_l = rt_nrand_Upu32_Yd_f_pw(&Sensors_DW.RandSeed_h) *
      0.0031622776601683794;

    /* InitializeConditions for RandomNumber: '<S10>/Random Number' */
    Sensors_DW.RandSeed_f = 1144108930UL;
    Sensors_DW.NextOutput_lh = rt_nrand_Upu32_Yd_f_pw(&Sensors_DW.RandSeed_f) *
      0.0031622776601683794;

    /* InitializeConditions for RandomNumber: '<S2>/Random Number' */
    Sensors_DW.RandSeed_fm = 1144108930UL;
    Sensors_DW.NextOutput_n = rt_nrand_Upu32_Yd_f_pw(&Sensors_DW.RandSeed_fm) *
      0.01;

    /* InitializeConditions for RandomNumber: '<S6>/Random Number' */
    y = 1144108930UL;
    y1 = rt_nrand_Upu32_Yd_f_pw(&y) * 15.0;
    Sensors_DW.NextOutput_o[0] = y1;
    Sensors_DW.RandSeed_e[0] = y;
    y = 1144108930UL;
    y1 = rt_nrand_Upu32_Yd_f_pw(&y) * 15.0;
    Sensors_DW.NextOutput_o[1] = y1;
    Sensors_DW.RandSeed_e[1] = y;
    y = 1144108930UL;
    y1 = rt_nrand_Upu32_Yd_f_pw(&y) * 20.0;
    Sensors_DW.NextOutput_o[2] = y1;
    Sensors_DW.RandSeed_e[2] = y;

    /* InitializeConditions for RandomNumber: '<S6>/Random Number1' */
    y = 1144108930UL;
    y1 = rt_nrand_Upu32_Yd_f_pw(&y) * 0.031622776601683791;
    Sensors_DW.NextOutput_h[0] = y1;
    Sensors_DW.RandSeed_i[0] = y;
    y = 1144108930UL;
    y1 = rt_nrand_Upu32_Yd_f_pw(&y) * 0.031622776601683791;
    Sensors_DW.NextOutput_h[1] = y1;
    Sensors_DW.RandSeed_i[1] = y;
    y = 1144108930UL;
    y1 = rt_nrand_Upu32_Yd_f_pw(&y) * 0.031622776601683791;
    Sensors_DW.NextOutput_h[2] = y1;
    Sensors_DW.RandSeed_i[2] = y;

    /* InitializeConditions for RandomNumber: '<S9>/Random Number' */
    Sensors_DW.RandSeed_p = 1144108930UL;
    Sensors_DW.NextOutput_a = rt_nrand_Upu32_Yd_f_pw(&Sensors_DW.RandSeed_p) *
      0.0031622776601683794;

    /* InitializeConditions for RandomNumber: '<S3>/White Noise' */
    Sensors_DW.RandSeed_l = 931168259UL;
    Sensors_DW.NextOutput_lh3 = rt_nrand_Upu32_Yd_f_pw(&Sensors_DW.RandSeed_l);

    /* InitializeConditions for RandomNumber: '<S81>/White Noise' */
    y = 1373044741UL;
    y1 = rt_nrand_Upu32_Yd_f_pw(&y);
    Sensors_DW.NextOutput_k[0] = y1;
    Sensors_DW.RandSeed_ls[0] = y;
    y = 411009029UL;
    y1 = rt_nrand_Upu32_Yd_f_pw(&y);
    Sensors_DW.NextOutput_k[1] = y1;
    Sensors_DW.RandSeed_ls[1] = y;
    y = 1845526542UL;
    y1 = rt_nrand_Upu32_Yd_f_pw(&y);
    Sensors_DW.NextOutput_k[2] = y1;
    Sensors_DW.RandSeed_ls[2] = y;

    /* InitializeConditions for RandomNumber: '<S98>/White Noise' */
    y = 1689616386UL;
    y1 = rt_nrand_Upu32_Yd_f_pw(&y);
    Sensors_DW.NextOutput_p[0] = y1;
    Sensors_DW.RandSeed_j[0] = y;
    y = 1998225409UL;
    y1 = rt_nrand_Upu32_Yd_f_pw(&y);
    Sensors_DW.NextOutput_p[1] = y1;
    Sensors_DW.RandSeed_j[1] = y;
    y = 1181220867UL;
    y1 = rt_nrand_Upu32_Yd_f_pw(&y);
    Sensors_DW.NextOutput_p[2] = y1;
    Sensors_DW.RandSeed_j[2] = y;
  }
}

/* Model terminate function */
void Sensors_terminate(void)
{
  /* (no terminate code required) */
}
