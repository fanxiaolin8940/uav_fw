/*
 * SYSTEM0_data.c
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

/* Invariant block signals (auto storage) */
const ConstB_SYSTEM0_T SYSTEM0_ConstB = {
  1.0
  ,                                    /* '<S19>/q0' */
  0.0
  ,                                    /* '<S19>/q1' */
  0.0
  ,                                    /* '<S19>/q2' */
  0.0
  ,                                    /* '<S19>/q3' */

  { 0.0209, 0.0, 0.0, 0.0, 0.0209, 0.0, 0.0, 0.0, 0.0509 }
  ,                                    /* '<S11>/Selector' */

  { 0.0209, 0.0, 0.0, 0.0, 0.0209, 0.0, 0.0, 0.0, 0.0509 }
  ,                                    /* '<S11>/Selector2' */
  0.0
  /* '<S119>/Switch' */
};

/* Constant parameters (auto storage) */
const ConstP_SYSTEM0_T SYSTEM0_ConstP = {
  /* Pooled Parameter (Mixed Expressions)
   * Referenced by:
   *   '<S129>/Scale factors & Cross-coupling  errors'
   *   '<S130>/Scale factors & Cross-coupling  errors '
   */
  { 0.98, 0.0, 0.0, 0.0, 0.98, 0.0, 0.0, 0.0, 0.98 }
};
