/*
 * Dynamics_data.c
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

/* Invariant block signals (auto storage) */
const ConstB_Dynamics_T Dynamics_ConstB = {
  1.0
  ,                                    /* '<S15>/q0' */
  0.0
  ,                                    /* '<S15>/q1' */
  0.0
  ,                                    /* '<S15>/q2' */
  0.0
  ,                                    /* '<S15>/q3' */

  { 0.0209, 0.0, 0.0, 0.0, 0.0209, 0.0, 0.0, 0.0, 0.0509 }
  ,                                    /* '<S7>/Selector' */

  { 0.0209, 0.0, 0.0, 0.0, 0.0209, 0.0, 0.0, 0.0, 0.0509 }
  /* '<S7>/Selector2' */
};
