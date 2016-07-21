/*
 * Sensors_data.c
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

/* Invariant block signals (auto storage) */
const ConstB_Sensors_T Sensors_ConstB = {
  0.0
  ,                                    /* '<S68>/Switch' */
  0.0
  /* '<S49>/Switch' */
};

/* Constant parameters (auto storage) */
const ConstP_Sensors_T Sensors_ConstP = {
  /* Pooled Parameter (Mixed Expressions)
   * Referenced by:
   *   '<S78>/Scale factors & Cross-coupling  errors'
   *   '<S79>/Scale factors & Cross-coupling  errors '
   */
  { 0.98, 0.0, 0.0, 0.0, 0.98, 0.0, 0.0, 0.0, 0.98 }
};
