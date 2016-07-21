/*
 * Sensors_private.h
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

#ifndef RTW_HEADER_Sensors_private_h_
#define RTW_HEADER_Sensors_private_h_
#include "rtwtypes.h"
#include "multiword_types.h"

extern real_T rt_urand_Upu32_Yd_f_pw(uint32_T *u);
extern real_T rt_nrand_Upu32_Yd_f_pw(uint32_T *u);
extern real_T rt_roundd(real_T u);
extern real_T rt_modd(real_T u0, real_T u1);
extern void Sensors_step0(void);
extern void Sensors_step1(void);
extern void Sensors_step2(void);

#endif                                 /* RTW_HEADER_Sensors_private_h_ */
