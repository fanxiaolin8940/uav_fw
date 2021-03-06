/*
 * DynModel_data.c
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

/* Invariant block signals (auto storage) */
const ConstB_DynModel_T DynModel_ConstB = {
  10.423324
  ,                                    /* '<S117>/Switch' */
  1.0
  ,                                    /* '<S17>/q0' */
  0.0
  ,                                    /* '<S17>/q1' */
  0.0
  ,                                    /* '<S17>/q2' */
  0.0
  ,                                    /* '<S17>/q3' */
  10.423324
  ,                                    /* '<S98>/Switch' */

  { 0.0209, 0.0, 0.0, 0.0, 0.0209, 0.0, 0.0, 0.0, 0.0509 }
  ,                                    /* '<S9>/Selector' */

  { 0.0209, 0.0, 0.0, 0.0, 0.0209, 0.0, 0.0, 0.0, 0.0509 }
  /* '<S9>/Selector2' */
};

/* Constant parameters (auto storage) */
const ConstP_DynModel_T DynModel_ConstP = {
  /* Pooled Parameter (Mixed Expressions)
   * Referenced by:
   *   '<S127>/Scale factors & Cross-coupling  errors'
   *   '<S128>/Scale factors & Cross-coupling  errors '
   */
  { 0.98, 0.0, 0.0, 0.0, 0.98, 0.0, 0.0, 0.0, 0.98 }
};
