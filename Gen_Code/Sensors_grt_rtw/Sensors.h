/*
 * Sensors.h
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

#ifndef RTW_HEADER_Sensors_h_
#define RTW_HEADER_Sensors_h_
#include <float.h>
#include <math.h>
#include <string.h>
#ifndef Sensors_COMMON_INCLUDES_
# define Sensors_COMMON_INCLUDES_
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#endif                                 /* Sensors_COMMON_INCLUDES_ */

#include "Sensors_types.h"

/* Shared type includes */
#include "multiword_types.h"

/* Macros for accessing real-time model data structure */
#ifndef rtmCounterLimit
# define rtmCounterLimit(rtm, idx)     ((rtm)->Timing.TaskCounters.cLimit[(idx)])
#endif

#ifndef rtmGetErrorStatus
# define rtmGetErrorStatus(rtm)        ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
# define rtmSetErrorStatus(rtm, val)   ((rtm)->errorStatus = (val))
#endif

#ifndef rtmStepTask
# define rtmStepTask(rtm, idx)         ((rtm)->Timing.TaskCounters.TID[(idx)] == 0)
#endif

#ifndef rtmTaskCounter
# define rtmTaskCounter(rtm, idx)      ((rtm)->Timing.TaskCounters.TID[(idx)])
#endif

/* Block signals (auto storage) */
typedef struct {
  real_T RateTransition;               /* '<S12>/Rate Transition' */
  real_T RateTransition_b;             /* '<S10>/Rate Transition' */
  real_T RateTransition1[3];           /* '<S6>/Rate Transition1' */
  real_T RateTransition_e;             /* '<S6>/Rate Transition' */
  real_T Merge[4];                     /* '<S4>/Merge' */
} B_Sensors_T;

/* Block states (auto storage) for system '<Root>' */
typedef struct {
  real_T NextOutput;                   /* '<S14>/Random Number' */
  real_T NextOutput_l;                 /* '<S12>/Random Number' */
  real_T NextOutput_lh;                /* '<S10>/Random Number' */
  real_T NextOutput_n;                 /* '<S2>/Random Number' */
  real_T NextOutput_o[3];              /* '<S6>/Random Number' */
  real_T NextOutput_h[3];              /* '<S6>/Random Number1' */
  real_T NextOutput_a;                 /* '<S9>/Random Number' */
  real_T NextOutput_lh3;               /* '<S3>/White Noise' */
  real_T NextOutput_k[3];              /* '<S81>/White Noise' */
  real_T NextOutput_p[3];              /* '<S98>/White Noise' */
  uint32_T RandSeed;                   /* '<S14>/Random Number' */
  uint32_T RandSeed_h;                 /* '<S12>/Random Number' */
  uint32_T RandSeed_f;                 /* '<S10>/Random Number' */
  uint32_T RandSeed_fm;                /* '<S2>/Random Number' */
  uint32_T RandSeed_e[3];              /* '<S6>/Random Number' */
  uint32_T RandSeed_i[3];              /* '<S6>/Random Number1' */
  uint32_T RandSeed_p;                 /* '<S9>/Random Number' */
  uint32_T RandSeed_l;                 /* '<S3>/White Noise' */
  uint32_T RandSeed_ls[3];             /* '<S81>/White Noise' */
  uint32_T RandSeed_j[3];              /* '<S98>/White Noise' */
} DW_Sensors_T;

/* Invariant block signals (auto storage) */
typedef struct {
  const real_T Switch_d;               /* '<S68>/Switch' */
  const real_T Switch_b;               /* '<S49>/Switch' */
} ConstB_Sensors_T;

/* Constant parameters (auto storage) */
typedef struct {
  /* Pooled Parameter (Mixed Expressions)
   * Referenced by:
   *   '<S78>/Scale factors & Cross-coupling  errors'
   *   '<S79>/Scale factors & Cross-coupling  errors '
   */
  real_T pooled22[9];
} ConstP_Sensors_T;

/* External inputs (root inport signals with auto storage) */
typedef struct {
  real_T V_e[3];                       /* '<Root>/V_e' */
  real_T X_e[3];                       /* '<Root>/X_e' */
  real_T rpy[3];                       /* '<Root>/rpy' */
  real_T DCM[9];                       /* '<Root>/DCM' */
  real_T V_b[3];                       /* '<Root>/V_b' */
  real_T Omega[3];                     /* '<Root>/Omega ' */
  real_T Omega_dot[3];                 /* '<Root>/Omega_dot' */
  real_T A_b[3];                       /* '<Root>/A_b' */
} ExtU_Sensors_T;

/* External outputs (root outports fed by signals with auto storage) */
typedef struct {
  real_T Temp;                         /* '<Root>/Temp' */
  real_T Press;                        /* '<Root>/Press' */
  real_T rho;                          /* '<Root>/rho' */
  real_T diff_Press;                   /* '<Root>/diff_Press' */
  real_T Baro_alt;                     /* '<Root>/Baro_alt' */
  real_T Lat_meas;                     /* '<Root>/Lat_meas' */
  real_T Lon_meas;                     /* '<Root>/Lon_meas' */
  real_T Alt_meas;                     /* '<Root>/Alt_meas' */
  real_T V_meas[3];                    /* '<Root>/V_meas' */
  real_T COG;                          /* '<Root>/COG' */
  real_T Real_Lat_Lon[2];              /* '<Root>/Real_Lat_Lon' */
  real_T Real_Alt;                     /* '<Root>/Real_Alt' */
  real_T Magn_meas[3];                 /* '<Root>/Magn_meas' */
  real_T RPY[3];                       /* '<Root>/RPY' */
  real_T A_meas[3];                    /* '<Root>/A_meas' */
  real_T Omega_meas[3];                /* '<Root>/Omega_meas' */
  real_T q[4];                         /* '<Root>/q' */
  real_T Sonar_alt;                    /* '<Root>/Sonar_alt' */
} ExtY_Sensors_T;

/* Real-time Model Data Structure */
struct tag_RTM_Sensors_T {
  const char_T *errorStatus;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    struct {
      uint16_T TID[3];
      uint16_T cLimit[3];
    } TaskCounters;

    struct {
      uint16_T TID0_1;
      uint16_T TID0_2;
    } RateInteraction;
  } Timing;
};

/* Block signals (auto storage) */
extern B_Sensors_T Sensors_B;

/* Block states (auto storage) */
extern DW_Sensors_T Sensors_DW;

/* External inputs (root inport signals with auto storage) */
extern ExtU_Sensors_T Sensors_U;

/* External outputs (root outports fed by signals with auto storage) */
extern ExtY_Sensors_T Sensors_Y;
extern const ConstB_Sensors_T Sensors_ConstB;/* constant block i/o */

/* Constant parameters (auto storage) */
extern const ConstP_Sensors_T Sensors_ConstP;

/* Model entry point functions */
extern void Sensors_initialize(void);
extern void Sensors_step(int_T tid);
extern void Sensors_terminate(void);

/* Real-time Model object */
extern RT_MODEL_Sensors_T *const Sensors_M;

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Note that this particular code originates from a subsystem build,
 * and has its own system numbers different from the parent model.
 * Refer to the system hierarchy for this subsystem below, and use the
 * MATLAB hilite_system command to trace the generated code back
 * to the parent model.  For example,
 *
 * hilite_system('simulator_mod/dynamic_model/Sensors')    - opens subsystem simulator_mod/dynamic_model/Sensors
 * hilite_system('simulator_mod/dynamic_model/Sensors/Kp') - opens and selects block Kp
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'simulator_mod/dynamic_model'
 * '<S1>'   : 'simulator_mod/dynamic_model/Sensors'
 * '<S2>'   : 'simulator_mod/dynamic_model/Sensors/Altitude_From_Pressure'
 * '<S3>'   : 'simulator_mod/dynamic_model/Sensors/Band-Limited White Noise1'
 * '<S4>'   : 'simulator_mod/dynamic_model/Sensors/Direction Cosine Matrix  to Quaternions'
 * '<S5>'   : 'simulator_mod/dynamic_model/Sensors/Flat Earth to LLA'
 * '<S6>'   : 'simulator_mod/dynamic_model/Sensors/GPS'
 * '<S7>'   : 'simulator_mod/dynamic_model/Sensors/ISA Atmosphere Model_'
 * '<S8>'   : 'simulator_mod/dynamic_model/Sensors/Invensense MPU-6000'
 * '<S9>'   : 'simulator_mod/dynamic_model/Sensors/Magnetometer'
 * '<S10>'  : 'simulator_mod/dynamic_model/Sensors/Pitot Sensor'
 * '<S11>'  : 'simulator_mod/dynamic_model/Sensors/PixFlow'
 * '<S12>'  : 'simulator_mod/dynamic_model/Sensors/Pressure Sensor'
 * '<S13>'  : 'simulator_mod/dynamic_model/Sensors/Sonar'
 * '<S14>'  : 'simulator_mod/dynamic_model/Sensors/Temperature Sensor'
 * '<S15>'  : 'simulator_mod/dynamic_model/Sensors/Direction Cosine Matrix  to Quaternions/Negative Trace'
 * '<S16>'  : 'simulator_mod/dynamic_model/Sensors/Direction Cosine Matrix  to Quaternions/Positive Trace'
 * '<S17>'  : 'simulator_mod/dynamic_model/Sensors/Direction Cosine Matrix  to Quaternions/trace(DCM)'
 * '<S18>'  : 'simulator_mod/dynamic_model/Sensors/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(1,1)'
 * '<S19>'  : 'simulator_mod/dynamic_model/Sensors/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(2,2)'
 * '<S20>'  : 'simulator_mod/dynamic_model/Sensors/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(3,3)'
 * '<S21>'  : 'simulator_mod/dynamic_model/Sensors/Direction Cosine Matrix  to Quaternions/Negative Trace/diag(DCM)'
 * '<S22>'  : 'simulator_mod/dynamic_model/Sensors/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(1,1)/cos(phi)sin(theta)cos(psi) + sin(phi)sin(psi) -sin(theta)'
 * '<S23>'  : 'simulator_mod/dynamic_model/Sensors/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(1,1)/cos(theta)sin(phi) - (cos(phi)sin(theta)sin(psi) - sin(phi)cos(psi))'
 * '<S24>'  : 'simulator_mod/dynamic_model/Sensors/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(1,1)/cos(theta)sin(psi) + (sin(phi)sin(theta)cos(psi) - cos(phi)sin(psi))'
 * '<S25>'  : 'simulator_mod/dynamic_model/Sensors/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(1,1)/if s~=0; s=0.5//s'
 * '<S26>'  : 'simulator_mod/dynamic_model/Sensors/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(1,1)/u(1) -(u(5)+u(9)) +1'
 * '<S27>'  : 'simulator_mod/dynamic_model/Sensors/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(2,2)/cos(phi)sin(theta)cos(psi) + sin(phi)sin(psi) +sin(theta)'
 * '<S28>'  : 'simulator_mod/dynamic_model/Sensors/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(2,2)/cos(theta)sin(phi) + (cos(phi)sin(theta)sin(psi) - sin(phi)cos(psi))'
 * '<S29>'  : 'simulator_mod/dynamic_model/Sensors/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(2,2)/cos(theta)sin(psi) + (sin(phi)sin(theta)cos(psi) - cos(phi)sin(psi))'
 * '<S30>'  : 'simulator_mod/dynamic_model/Sensors/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(2,2)/if s~=0; s=0.5//s'
 * '<S31>'  : 'simulator_mod/dynamic_model/Sensors/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(2,2)/u(5) -(u(1)+u(9)) +1'
 * '<S32>'  : 'simulator_mod/dynamic_model/Sensors/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(3,3)/cos(phi)sin(theta)cos(psi) + sin(phi)sin(psi) -sin(theta)'
 * '<S33>'  : 'simulator_mod/dynamic_model/Sensors/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(3,3)/cos(theta)sin(phi) + (cos(phi)sin(theta)sin(psi) - sin(phi)cos(psi))'
 * '<S34>'  : 'simulator_mod/dynamic_model/Sensors/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(3,3)/cos(theta)sin(psi) - (sin(phi)sin(theta)cos(psi) - cos(phi)sin(psi))'
 * '<S35>'  : 'simulator_mod/dynamic_model/Sensors/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(3,3)/if s~=0; s=0.5//s'
 * '<S36>'  : 'simulator_mod/dynamic_model/Sensors/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(3,3)/u(9) -(u(1)+u(5)) +1'
 * '<S37>'  : 'simulator_mod/dynamic_model/Sensors/Direction Cosine Matrix  to Quaternions/Positive Trace/cos(phi)sin(theta)cos(psi) + sin(phi)sin(psi) +sin(theta)'
 * '<S38>'  : 'simulator_mod/dynamic_model/Sensors/Direction Cosine Matrix  to Quaternions/Positive Trace/cos(theta)sin(phi) - (cos(phi)sin(theta)sin(psi) - sin(phi)cos(psi))'
 * '<S39>'  : 'simulator_mod/dynamic_model/Sensors/Direction Cosine Matrix  to Quaternions/Positive Trace/cos(theta)sin(psi) - (sin(phi)sin(theta)cos(psi) - cos(phi)sin(psi))'
 * '<S40>'  : 'simulator_mod/dynamic_model/Sensors/Flat Earth to LLA/LatLong wrap'
 * '<S41>'  : 'simulator_mod/dynamic_model/Sensors/Flat Earth to LLA/LatLong wrap1'
 * '<S42>'  : 'simulator_mod/dynamic_model/Sensors/Flat Earth to LLA/LongLat_offset'
 * '<S43>'  : 'simulator_mod/dynamic_model/Sensors/Flat Earth to LLA/pos_deg'
 * '<S44>'  : 'simulator_mod/dynamic_model/Sensors/Flat Earth to LLA/LatLong wrap/Latitude Wrap 90'
 * '<S45>'  : 'simulator_mod/dynamic_model/Sensors/Flat Earth to LLA/LatLong wrap/Wrap Longitude'
 * '<S46>'  : 'simulator_mod/dynamic_model/Sensors/Flat Earth to LLA/LatLong wrap/Latitude Wrap 90/Compare To Constant'
 * '<S47>'  : 'simulator_mod/dynamic_model/Sensors/Flat Earth to LLA/LatLong wrap/Latitude Wrap 90/Wrap Angle 180'
 * '<S48>'  : 'simulator_mod/dynamic_model/Sensors/Flat Earth to LLA/LatLong wrap1/Latitude Wrap 90'
 * '<S49>'  : 'simulator_mod/dynamic_model/Sensors/Flat Earth to LLA/LatLong wrap1/Wrap Longitude'
 * '<S50>'  : 'simulator_mod/dynamic_model/Sensors/Flat Earth to LLA/LatLong wrap1/Latitude Wrap 90/Compare To Constant'
 * '<S51>'  : 'simulator_mod/dynamic_model/Sensors/Flat Earth to LLA/LatLong wrap1/Latitude Wrap 90/Wrap Angle 180'
 * '<S52>'  : 'simulator_mod/dynamic_model/Sensors/Flat Earth to LLA/LongLat_offset/Find Radian//Distance'
 * '<S53>'  : 'simulator_mod/dynamic_model/Sensors/Flat Earth to LLA/LongLat_offset/pos_deg'
 * '<S54>'  : 'simulator_mod/dynamic_model/Sensors/Flat Earth to LLA/LongLat_offset/Find Radian//Distance/Angle Conversion2'
 * '<S55>'  : 'simulator_mod/dynamic_model/Sensors/Flat Earth to LLA/LongLat_offset/Find Radian//Distance/denom'
 * '<S56>'  : 'simulator_mod/dynamic_model/Sensors/Flat Earth to LLA/LongLat_offset/Find Radian//Distance/e'
 * '<S57>'  : 'simulator_mod/dynamic_model/Sensors/Flat Earth to LLA/LongLat_offset/Find Radian//Distance/e^4'
 * '<S58>'  : 'simulator_mod/dynamic_model/Sensors/GPS/Flat Earth to LLA1'
 * '<S59>'  : 'simulator_mod/dynamic_model/Sensors/GPS/Flat Earth to LLA1/LatLong wrap'
 * '<S60>'  : 'simulator_mod/dynamic_model/Sensors/GPS/Flat Earth to LLA1/LatLong wrap1'
 * '<S61>'  : 'simulator_mod/dynamic_model/Sensors/GPS/Flat Earth to LLA1/LongLat_offset'
 * '<S62>'  : 'simulator_mod/dynamic_model/Sensors/GPS/Flat Earth to LLA1/pos_deg'
 * '<S63>'  : 'simulator_mod/dynamic_model/Sensors/GPS/Flat Earth to LLA1/LatLong wrap/Latitude Wrap 90'
 * '<S64>'  : 'simulator_mod/dynamic_model/Sensors/GPS/Flat Earth to LLA1/LatLong wrap/Wrap Longitude'
 * '<S65>'  : 'simulator_mod/dynamic_model/Sensors/GPS/Flat Earth to LLA1/LatLong wrap/Latitude Wrap 90/Compare To Constant'
 * '<S66>'  : 'simulator_mod/dynamic_model/Sensors/GPS/Flat Earth to LLA1/LatLong wrap/Latitude Wrap 90/Wrap Angle 180'
 * '<S67>'  : 'simulator_mod/dynamic_model/Sensors/GPS/Flat Earth to LLA1/LatLong wrap1/Latitude Wrap 90'
 * '<S68>'  : 'simulator_mod/dynamic_model/Sensors/GPS/Flat Earth to LLA1/LatLong wrap1/Wrap Longitude'
 * '<S69>'  : 'simulator_mod/dynamic_model/Sensors/GPS/Flat Earth to LLA1/LatLong wrap1/Latitude Wrap 90/Compare To Constant'
 * '<S70>'  : 'simulator_mod/dynamic_model/Sensors/GPS/Flat Earth to LLA1/LatLong wrap1/Latitude Wrap 90/Wrap Angle 180'
 * '<S71>'  : 'simulator_mod/dynamic_model/Sensors/GPS/Flat Earth to LLA1/LongLat_offset/Find Radian//Distance'
 * '<S72>'  : 'simulator_mod/dynamic_model/Sensors/GPS/Flat Earth to LLA1/LongLat_offset/pos_deg'
 * '<S73>'  : 'simulator_mod/dynamic_model/Sensors/GPS/Flat Earth to LLA1/LongLat_offset/Find Radian//Distance/Angle Conversion2'
 * '<S74>'  : 'simulator_mod/dynamic_model/Sensors/GPS/Flat Earth to LLA1/LongLat_offset/Find Radian//Distance/denom'
 * '<S75>'  : 'simulator_mod/dynamic_model/Sensors/GPS/Flat Earth to LLA1/LongLat_offset/Find Radian//Distance/e'
 * '<S76>'  : 'simulator_mod/dynamic_model/Sensors/GPS/Flat Earth to LLA1/LongLat_offset/Find Radian//Distance/e^4'
 * '<S77>'  : 'simulator_mod/dynamic_model/Sensors/Invensense MPU-6000/Acceleration Conversion'
 * '<S78>'  : 'simulator_mod/dynamic_model/Sensors/Invensense MPU-6000/Three-axis Accelerometer'
 * '<S79>'  : 'simulator_mod/dynamic_model/Sensors/Invensense MPU-6000/Three-axis Gyroscope'
 * '<S80>'  : 'simulator_mod/dynamic_model/Sensors/Invensense MPU-6000/Three-axis Accelerometer/Dynamics'
 * '<S81>'  : 'simulator_mod/dynamic_model/Sensors/Invensense MPU-6000/Three-axis Accelerometer/Random bias'
 * '<S82>'  : 'simulator_mod/dynamic_model/Sensors/Invensense MPU-6000/Three-axis Accelerometer/w x (w x d)'
 * '<S83>'  : 'simulator_mod/dynamic_model/Sensors/Invensense MPU-6000/Three-axis Accelerometer/wdot x d'
 * '<S84>'  : 'simulator_mod/dynamic_model/Sensors/Invensense MPU-6000/Three-axis Accelerometer/Dynamics/No Dynamics'
 * '<S85>'  : 'simulator_mod/dynamic_model/Sensors/Invensense MPU-6000/Three-axis Accelerometer/Dynamics/Second-order Dynamics'
 * '<S86>'  : 'simulator_mod/dynamic_model/Sensors/Invensense MPU-6000/Three-axis Accelerometer/Dynamics/Second-order Dynamics/Transfer Fcn X'
 * '<S87>'  : 'simulator_mod/dynamic_model/Sensors/Invensense MPU-6000/Three-axis Accelerometer/Dynamics/Second-order Dynamics/Transfer Fcn Y'
 * '<S88>'  : 'simulator_mod/dynamic_model/Sensors/Invensense MPU-6000/Three-axis Accelerometer/Dynamics/Second-order Dynamics/Transfer Fcn Z'
 * '<S89>'  : 'simulator_mod/dynamic_model/Sensors/Invensense MPU-6000/Three-axis Accelerometer/w x (w x d)/w x (w x d)'
 * '<S90>'  : 'simulator_mod/dynamic_model/Sensors/Invensense MPU-6000/Three-axis Accelerometer/w x (w x d)/w x d'
 * '<S91>'  : 'simulator_mod/dynamic_model/Sensors/Invensense MPU-6000/Three-axis Accelerometer/w x (w x d)/w x (w x d)/Subsystem'
 * '<S92>'  : 'simulator_mod/dynamic_model/Sensors/Invensense MPU-6000/Three-axis Accelerometer/w x (w x d)/w x (w x d)/Subsystem1'
 * '<S93>'  : 'simulator_mod/dynamic_model/Sensors/Invensense MPU-6000/Three-axis Accelerometer/w x (w x d)/w x d/Subsystem'
 * '<S94>'  : 'simulator_mod/dynamic_model/Sensors/Invensense MPU-6000/Three-axis Accelerometer/w x (w x d)/w x d/Subsystem1'
 * '<S95>'  : 'simulator_mod/dynamic_model/Sensors/Invensense MPU-6000/Three-axis Accelerometer/wdot x d/Subsystem'
 * '<S96>'  : 'simulator_mod/dynamic_model/Sensors/Invensense MPU-6000/Three-axis Accelerometer/wdot x d/Subsystem1'
 * '<S97>'  : 'simulator_mod/dynamic_model/Sensors/Invensense MPU-6000/Three-axis Gyroscope/Dynamics'
 * '<S98>'  : 'simulator_mod/dynamic_model/Sensors/Invensense MPU-6000/Three-axis Gyroscope/Random bias'
 * '<S99>'  : 'simulator_mod/dynamic_model/Sensors/Invensense MPU-6000/Three-axis Gyroscope/Dynamics/No Dynamics'
 * '<S100>' : 'simulator_mod/dynamic_model/Sensors/Invensense MPU-6000/Three-axis Gyroscope/Dynamics/Second-order Dynamics'
 * '<S101>' : 'simulator_mod/dynamic_model/Sensors/Invensense MPU-6000/Three-axis Gyroscope/Dynamics/Second-order Dynamics/Transfer Fcn X'
 * '<S102>' : 'simulator_mod/dynamic_model/Sensors/Invensense MPU-6000/Three-axis Gyroscope/Dynamics/Second-order Dynamics/Transfer Fcn Y'
 * '<S103>' : 'simulator_mod/dynamic_model/Sensors/Invensense MPU-6000/Three-axis Gyroscope/Dynamics/Second-order Dynamics/Transfer Fcn Z'
 * '<S104>' : 'simulator_mod/dynamic_model/Sensors/Magnetometer/Measurement Noise'
 * '<S105>' : 'simulator_mod/dynamic_model/Sensors/Pitot Sensor/Measurement Noise'
 * '<S106>' : 'simulator_mod/dynamic_model/Sensors/PixFlow/Band-Limited White Noise4'
 * '<S107>' : 'simulator_mod/dynamic_model/Sensors/Pressure Sensor/Measurement Noise'
 * '<S108>' : 'simulator_mod/dynamic_model/Sensors/Sonar/Band-Limited White Noise3'
 * '<S109>' : 'simulator_mod/dynamic_model/Sensors/Sonar/MATLAB Function'
 * '<S110>' : 'simulator_mod/dynamic_model/Sensors/Temperature Sensor/Measurement Noise'
 */
#endif                                 /* RTW_HEADER_Sensors_h_ */
