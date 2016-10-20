/*
 * DynModel_data.c
 *
 * Code generation for model "DynModel".
 *
 * Model version              : 1.676
 * Simulink Coder version : 8.8 (R2015a) 09-Feb-2015
 * C source code generated on : Thu Oct 20 20:45:09 2016
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

/* Invariant block signals (auto storage) */
const ConstB_DynModel_T DynModel_ConstB = {
  10.3192456
  ,                                    /* '<S125>/Switch' */
  1.0
  ,                                    /* '<S19>/q0' */
  0.0
  ,                                    /* '<S19>/q1' */
  0.0
  ,                                    /* '<S19>/q2' */
  0.0
  ,                                    /* '<S19>/q3' */
  10.3192456
  ,                                    /* '<S106>/Switch' */

  { 0.0409, 0.0, 0.0, 0.0, 0.0409, 0.0, 0.0, 0.0, 0.0509 }
  ,                                    /* '<S11>/Selector' */

  { 0.0409, 0.0, 0.0, 0.0, 0.0409, 0.0, 0.0, 0.0, 0.0509 }
  /* '<S11>/Selector2' */
};

/* Constant parameters (auto storage) */
const ConstP_DynModel_T DynModel_ConstP = {
  /* Pooled Parameter (Mixed Expressions)
   * Referenced by:
   *   '<S135>/Scale factors & Cross-coupling  errors'
   *   '<S136>/Scale factors & Cross-coupling  errors '
   */
  { 0.98, 0.0, 0.0, 0.0, 0.98, 0.0, 0.0, 0.0, 0.98 }
};
