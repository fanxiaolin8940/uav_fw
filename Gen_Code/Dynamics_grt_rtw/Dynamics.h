/*
 * Dynamics.h
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

#ifndef RTW_HEADER_Dynamics_h_
#define RTW_HEADER_Dynamics_h_
#include <math.h>
#include <string.h>
#ifndef Dynamics_COMMON_INCLUDES_
# define Dynamics_COMMON_INCLUDES_
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#endif                                 /* Dynamics_COMMON_INCLUDES_ */

#include "Dynamics_types.h"

/* Shared type includes */
#include "multiword_types.h"

/* Macros for accessing real-time model data structure */
#ifndef rtmGetBlkStateChangeFlag
# define rtmGetBlkStateChangeFlag(rtm) ((rtm)->ModelData.blkStateChange)
#endif

#ifndef rtmSetBlkStateChangeFlag
# define rtmSetBlkStateChangeFlag(rtm, val) ((rtm)->ModelData.blkStateChange = (val))
#endif

#ifndef rtmGetContStateDisabled
# define rtmGetContStateDisabled(rtm)  ((rtm)->ModelData.contStateDisabled)
#endif

#ifndef rtmSetContStateDisabled
# define rtmSetContStateDisabled(rtm, val) ((rtm)->ModelData.contStateDisabled = (val))
#endif

#ifndef rtmGetContStates
# define rtmGetContStates(rtm)         ((rtm)->ModelData.contStates)
#endif

#ifndef rtmSetContStates
# define rtmSetContStates(rtm, val)    ((rtm)->ModelData.contStates = (val))
#endif

#ifndef rtmGetDerivCacheNeedsReset
# define rtmGetDerivCacheNeedsReset(rtm) ((rtm)->ModelData.derivCacheNeedsReset)
#endif

#ifndef rtmSetDerivCacheNeedsReset
# define rtmSetDerivCacheNeedsReset(rtm, val) ((rtm)->ModelData.derivCacheNeedsReset = (val))
#endif

#ifndef rtmGetIntgData
# define rtmGetIntgData(rtm)           ((rtm)->ModelData.intgData)
#endif

#ifndef rtmSetIntgData
# define rtmSetIntgData(rtm, val)      ((rtm)->ModelData.intgData = (val))
#endif

#ifndef rtmGetOdeF
# define rtmGetOdeF(rtm)               ((rtm)->ModelData.odeF)
#endif

#ifndef rtmSetOdeF
# define rtmSetOdeF(rtm, val)          ((rtm)->ModelData.odeF = (val))
#endif

#ifndef rtmGetOdeY
# define rtmGetOdeY(rtm)               ((rtm)->ModelData.odeY)
#endif

#ifndef rtmSetOdeY
# define rtmSetOdeY(rtm, val)          ((rtm)->ModelData.odeY = (val))
#endif

#ifndef rtmGetPeriodicContStateIndices
# define rtmGetPeriodicContStateIndices(rtm) ((rtm)->ModelData.periodicContStateIndices)
#endif

#ifndef rtmSetPeriodicContStateIndices
# define rtmSetPeriodicContStateIndices(rtm, val) ((rtm)->ModelData.periodicContStateIndices = (val))
#endif

#ifndef rtmGetPeriodicContStateRanges
# define rtmGetPeriodicContStateRanges(rtm) ((rtm)->ModelData.periodicContStateRanges)
#endif

#ifndef rtmSetPeriodicContStateRanges
# define rtmSetPeriodicContStateRanges(rtm, val) ((rtm)->ModelData.periodicContStateRanges = (val))
#endif

#ifndef rtmGetZCCacheNeedsReset
# define rtmGetZCCacheNeedsReset(rtm)  ((rtm)->ModelData.zCCacheNeedsReset)
#endif

#ifndef rtmSetZCCacheNeedsReset
# define rtmSetZCCacheNeedsReset(rtm, val) ((rtm)->ModelData.zCCacheNeedsReset = (val))
#endif

#ifndef rtmGetdX
# define rtmGetdX(rtm)                 ((rtm)->ModelData.derivs)
#endif

#ifndef rtmSetdX
# define rtmSetdX(rtm, val)            ((rtm)->ModelData.derivs = (val))
#endif

#ifndef rtmGetErrorStatus
# define rtmGetErrorStatus(rtm)        ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
# define rtmSetErrorStatus(rtm, val)   ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetStopRequested
# define rtmGetStopRequested(rtm)      ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequested
# define rtmSetStopRequested(rtm, val) ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetStopRequestedPtr
# define rtmGetStopRequestedPtr(rtm)   (&((rtm)->Timing.stopRequestedFlag))
#endif

#ifndef rtmGetT
# define rtmGetT(rtm)                  (rtmGetTPtr((rtm))[0])
#endif

/* Block signals (auto storage) */
typedef struct {
  real_T Product[3];                   /* '<S12>/Product' */
  real_T Memory2;                      /* '<S1>/Memory2' */
  real_T Product2[3];                  /* '<S7>/Product2' */
  real_T Add8[3];                      /* '<S1>/Add8' */
  real_T Gain4;                        /* '<S1>/Gain4' */
  real_T Sum[3];                       /* '<S2>/Sum' */
  real_T q0dot;                        /* '<S16>/q0dot' */
  real_T q1dot;                        /* '<S16>/q1dot' */
  real_T q2dot;                        /* '<S16>/q2dot' */
  real_T q3dot;                        /* '<S16>/q3dot' */
  real_T Product_o[4];                 /* '<S5>/Product' */
  real_T Sum2;                         /* '<S45>/Sum2' */
  real_T Sum2_j;                       /* '<S46>/Sum2' */
  real_T Sum2_c;                       /* '<S47>/Sum2' */
  real_T Sum2_p;                       /* '<S48>/Sum2' */
} B_Dynamics_T;

/* Block states (auto storage) for system '<Root>' */
typedef struct {
  real_T Memory2_PreviousInput;        /* '<S1>/Memory2' */
  real_T Product2_DWORK4[9];           /* '<S7>/Product2' */
  real_T discharge;                    /* '<S1>/LiPo Battery' */
  struct {
    int_T IcNeedsLoading;
  } q0q1q2q3_IWORK;                    /* '<S6>/q0 q1 q2 q3' */

  int_T IntegratorSecondOrder_MODE;    /* '<S45>/Integrator, Second-Order' */
  int_T IntegratorSecondOrder_MODE_p;  /* '<S46>/Integrator, Second-Order' */
  int_T IntegratorSecondOrder_MODE_a;  /* '<S47>/Integrator, Second-Order' */
  int_T IntegratorSecondOrder_MODE_pu; /* '<S48>/Integrator, Second-Order' */
  boolean_T IntegratorSecondOrder_DWORK1;/* '<S45>/Integrator, Second-Order' */
  boolean_T IntegratorSecondOrder_DWORK1_k;/* '<S46>/Integrator, Second-Order' */
  boolean_T IntegratorSecondOrder_DWORK1_o;/* '<S47>/Integrator, Second-Order' */
  boolean_T IntegratorSecondOrder_DWORK1_h;/* '<S48>/Integrator, Second-Order' */
} DW_Dynamics_T;

/* Continuous states (auto storage) */
typedef struct {
  real_T q0q1q2q3_CSTATE[4];           /* '<S6>/q0 q1 q2 q3' */
  real_T ubvbwb_CSTATE[3];             /* '<S2>/ub,vb,wb' */
  real_T xeyeze_CSTATE[3];             /* '<S2>/xe,ye,ze' */
  real_T pqr_CSTATE[3];                /* '<S2>/p,q,r ' */
  real_T IntegratorSecondOrder_CSTATE[2];/* '<S45>/Integrator, Second-Order' */
  real_T IntegratorSecondOrder_CSTATE_h[2];/* '<S46>/Integrator, Second-Order' */
  real_T IntegratorSecondOrder_CSTATE_n[2];/* '<S47>/Integrator, Second-Order' */
  real_T IntegratorSecondOrder_CSTATE_d[2];/* '<S48>/Integrator, Second-Order' */
} X_Dynamics_T;

/* State derivatives (auto storage) */
typedef struct {
  real_T q0q1q2q3_CSTATE[4];           /* '<S6>/q0 q1 q2 q3' */
  real_T ubvbwb_CSTATE[3];             /* '<S2>/ub,vb,wb' */
  real_T xeyeze_CSTATE[3];             /* '<S2>/xe,ye,ze' */
  real_T pqr_CSTATE[3];                /* '<S2>/p,q,r ' */
  real_T IntegratorSecondOrder_CSTATE[2];/* '<S45>/Integrator, Second-Order' */
  real_T IntegratorSecondOrder_CSTATE_h[2];/* '<S46>/Integrator, Second-Order' */
  real_T IntegratorSecondOrder_CSTATE_n[2];/* '<S47>/Integrator, Second-Order' */
  real_T IntegratorSecondOrder_CSTATE_d[2];/* '<S48>/Integrator, Second-Order' */
} XDot_Dynamics_T;

/* State disabled  */
typedef struct {
  boolean_T q0q1q2q3_CSTATE[4];        /* '<S6>/q0 q1 q2 q3' */
  boolean_T ubvbwb_CSTATE[3];          /* '<S2>/ub,vb,wb' */
  boolean_T xeyeze_CSTATE[3];          /* '<S2>/xe,ye,ze' */
  boolean_T pqr_CSTATE[3];             /* '<S2>/p,q,r ' */
  boolean_T IntegratorSecondOrder_CSTATE[2];/* '<S45>/Integrator, Second-Order' */
  boolean_T IntegratorSecondOrder_CSTATE_h[2];/* '<S46>/Integrator, Second-Order' */
  boolean_T IntegratorSecondOrder_CSTATE_n[2];/* '<S47>/Integrator, Second-Order' */
  boolean_T IntegratorSecondOrder_CSTATE_d[2];/* '<S48>/Integrator, Second-Order' */
} XDis_Dynamics_T;

/* Invariant block signals (auto storage) */
typedef struct {
  const real_T q0;                     /* '<S15>/q0' */
  const real_T q1;                     /* '<S15>/q1' */
  const real_T q2;                     /* '<S15>/q2' */
  const real_T q3;                     /* '<S15>/q3' */
  const real_T Selector[9];            /* '<S7>/Selector' */
  const real_T Selector2[9];           /* '<S7>/Selector2' */
} ConstB_Dynamics_T;

#ifndef ODE4_INTG
#define ODE4_INTG

/* ODE4 Integration Data */
typedef struct {
  real_T *y;                           /* output */
  real_T *f[4];                        /* derivatives */
} ODE4_IntgData;

#endif

/* External inputs (root inport signals with auto storage) */
typedef struct {
  real_T h;                            /* '<Root>/h' */
  real_T rho;                          /* '<Root>/rho' */
  real_T pwm1;                         /* '<Root>/pwm1' */
  real_T pwm2;                         /* '<Root>/pwm2' */
  real_T pwm3;                         /* '<Root>/pwm3' */
  real_T pwm4;                         /* '<Root>/pwm4' */
} ExtU_Dynamics_T;

/* External outputs (root outports fed by signals with auto storage) */
typedef struct {
  real_T V_e[3];                       /* '<Root>/V_e' */
  real_T X_e[3];                       /* '<Root>/X_e' */
  real_T rpy[3];                       /* '<Root>/rpy' */
  real_T DCM_be[9];                    /* '<Root>/DCM_be' */
  real_T V_b[3];                       /* '<Root>/V_b' */
  real_T Omega[3];                     /* '<Root>/Omega' */
  real_T Omega_dot[3];                 /* '<Root>/Omega_dot' */
  real_T A_b[3];                       /* '<Root>/A_b' */
  real_T Forces[3];                    /* '<Root>/Forces' */
  real_T Torques[3];                   /* '<Root>/Torques' */
  real_T Thursts[4];                   /* '<Root>/Thursts' */
  real_T Rotor_Speeds[4];              /* '<Root>/Rotor_Speeds' */
} ExtY_Dynamics_T;

/* Real-time Model Data Structure */
struct tag_RTM_Dynamics_T {
  const char_T *errorStatus;
  RTWSolverInfo solverInfo;

  /*
   * ModelData:
   * The following substructure contains information regarding
   * the data used in the model.
   */
  struct {
    X_Dynamics_T *contStates;
    int_T *periodicContStateIndices;
    real_T *periodicContStateRanges;
    real_T *derivs;
    boolean_T *contStateDisabled;
    boolean_T zCCacheNeedsReset;
    boolean_T derivCacheNeedsReset;
    boolean_T blkStateChange;
    real_T odeY[21];
    real_T odeF[4][21];
    ODE4_IntgData intgData;
  } ModelData;

  /*
   * Sizes:
   * The following substructure contains sizes information
   * for many of the model attributes such as inputs, outputs,
   * dwork, sample times, etc.
   */
  struct {
    int_T numContStates;
    int_T numPeriodicContStates;
    int_T numSampTimes;
  } Sizes;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    uint32_T clockTick0;
    uint32_T clockTickH0;
    time_T stepSize0;
    uint32_T clockTick1;
    uint32_T clockTickH1;
    boolean_T firstInitCondFlag;
    SimTimeStep simTimeStep;
    boolean_T stopRequestedFlag;
    time_T *t;
    time_T tArray[2];
  } Timing;
};

/* Block signals (auto storage) */
extern B_Dynamics_T Dynamics_B;

/* Continuous states (auto storage) */
extern X_Dynamics_T Dynamics_X;

/* Block states (auto storage) */
extern DW_Dynamics_T Dynamics_DW;

/* External inputs (root inport signals with auto storage) */
extern ExtU_Dynamics_T Dynamics_U;

/* External outputs (root outports fed by signals with auto storage) */
extern ExtY_Dynamics_T Dynamics_Y;
extern const ConstB_Dynamics_T Dynamics_ConstB;/* constant block i/o */

/* Model entry point functions */
extern void Dynamics_initialize(void);
extern void Dynamics_step(void);
extern void Dynamics_terminate(void);

/* Real-time Model object */
extern RT_MODEL_Dynamics_T *const Dynamics_M;

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
 * hilite_system('simulator_mod/dynamic_model/Dynamics')    - opens subsystem simulator_mod/dynamic_model/Dynamics
 * hilite_system('simulator_mod/dynamic_model/Dynamics/Kp') - opens and selects block Kp
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'simulator_mod/dynamic_model'
 * '<S1>'   : 'simulator_mod/dynamic_model/Dynamics'
 * '<S2>'   : 'simulator_mod/dynamic_model/Dynamics/6DoF (Quaternion)'
 * '<S3>'   : 'simulator_mod/dynamic_model/Dynamics/LiPo Battery'
 * '<S4>'   : 'simulator_mod/dynamic_model/Dynamics/Subsystem'
 * '<S5>'   : 'simulator_mod/dynamic_model/Dynamics/transfer_function'
 * '<S6>'   : 'simulator_mod/dynamic_model/Dynamics/6DoF (Quaternion)/Calculate DCM & Euler Angles'
 * '<S7>'   : 'simulator_mod/dynamic_model/Dynamics/6DoF (Quaternion)/Calculate omega_dot'
 * '<S8>'   : 'simulator_mod/dynamic_model/Dynamics/6DoF (Quaternion)/Determine Force,  Mass & Inertia'
 * '<S9>'   : 'simulator_mod/dynamic_model/Dynamics/6DoF (Quaternion)/Vbxw'
 * '<S10>'  : 'simulator_mod/dynamic_model/Dynamics/6DoF (Quaternion)/Velocity Conversion'
 * '<S11>'  : 'simulator_mod/dynamic_model/Dynamics/6DoF (Quaternion)/Velocity Conversion1'
 * '<S12>'  : 'simulator_mod/dynamic_model/Dynamics/6DoF (Quaternion)/transform to Inertial axes '
 * '<S13>'  : 'simulator_mod/dynamic_model/Dynamics/6DoF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to  Direction Cosine Matrix'
 * '<S14>'  : 'simulator_mod/dynamic_model/Dynamics/6DoF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to Rotation Angles'
 * '<S15>'  : 'simulator_mod/dynamic_model/Dynamics/6DoF (Quaternion)/Calculate DCM & Euler Angles/Rotation Angles to Quaternions'
 * '<S16>'  : 'simulator_mod/dynamic_model/Dynamics/6DoF (Quaternion)/Calculate DCM & Euler Angles/qdot'
 * '<S17>'  : 'simulator_mod/dynamic_model/Dynamics/6DoF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to  Direction Cosine Matrix/A11'
 * '<S18>'  : 'simulator_mod/dynamic_model/Dynamics/6DoF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to  Direction Cosine Matrix/A12'
 * '<S19>'  : 'simulator_mod/dynamic_model/Dynamics/6DoF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to  Direction Cosine Matrix/A13'
 * '<S20>'  : 'simulator_mod/dynamic_model/Dynamics/6DoF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to  Direction Cosine Matrix/A21'
 * '<S21>'  : 'simulator_mod/dynamic_model/Dynamics/6DoF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to  Direction Cosine Matrix/A22'
 * '<S22>'  : 'simulator_mod/dynamic_model/Dynamics/6DoF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to  Direction Cosine Matrix/A23'
 * '<S23>'  : 'simulator_mod/dynamic_model/Dynamics/6DoF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to  Direction Cosine Matrix/A31'
 * '<S24>'  : 'simulator_mod/dynamic_model/Dynamics/6DoF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to  Direction Cosine Matrix/A32'
 * '<S25>'  : 'simulator_mod/dynamic_model/Dynamics/6DoF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to  Direction Cosine Matrix/A33'
 * '<S26>'  : 'simulator_mod/dynamic_model/Dynamics/6DoF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to  Direction Cosine Matrix/Create Transformation Matrix'
 * '<S27>'  : 'simulator_mod/dynamic_model/Dynamics/6DoF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to  Direction Cosine Matrix/Quaternion Normalize'
 * '<S28>'  : 'simulator_mod/dynamic_model/Dynamics/6DoF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to  Direction Cosine Matrix/Quaternion Normalize/Quaternion Modulus'
 * '<S29>'  : 'simulator_mod/dynamic_model/Dynamics/6DoF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to  Direction Cosine Matrix/Quaternion Normalize/Quaternion Modulus/Quaternion Norm'
 * '<S30>'  : 'simulator_mod/dynamic_model/Dynamics/6DoF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to Rotation Angles/Quaternion Normalize'
 * '<S31>'  : 'simulator_mod/dynamic_model/Dynamics/6DoF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to Rotation Angles/Quaternion Normalize/Quaternion Modulus'
 * '<S32>'  : 'simulator_mod/dynamic_model/Dynamics/6DoF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to Rotation Angles/Quaternion Normalize/Quaternion Modulus/Quaternion Norm'
 * '<S33>'  : 'simulator_mod/dynamic_model/Dynamics/6DoF (Quaternion)/Calculate omega_dot/I x w'
 * '<S34>'  : 'simulator_mod/dynamic_model/Dynamics/6DoF (Quaternion)/Calculate omega_dot/I x w1'
 * '<S35>'  : 'simulator_mod/dynamic_model/Dynamics/6DoF (Quaternion)/Calculate omega_dot/wx(Iw)'
 * '<S36>'  : 'simulator_mod/dynamic_model/Dynamics/6DoF (Quaternion)/Calculate omega_dot/wx(Iw)/Subsystem'
 * '<S37>'  : 'simulator_mod/dynamic_model/Dynamics/6DoF (Quaternion)/Calculate omega_dot/wx(Iw)/Subsystem1'
 * '<S38>'  : 'simulator_mod/dynamic_model/Dynamics/6DoF (Quaternion)/Vbxw/Subsystem'
 * '<S39>'  : 'simulator_mod/dynamic_model/Dynamics/6DoF (Quaternion)/Vbxw/Subsystem1'
 * '<S40>'  : 'simulator_mod/dynamic_model/Dynamics/Subsystem/Multicopter_SFunction'
 * '<S41>'  : 'simulator_mod/dynamic_model/Dynamics/Subsystem/SYSTEM_MADE_WITH_BLOCKS'
 * '<S42>'  : 'simulator_mod/dynamic_model/Dynamics/Subsystem/multicopter'
 * '<S43>'  : 'simulator_mod/dynamic_model/Dynamics/Subsystem/Multicopter_SFunction/__InputSSForSFun__'
 * '<S44>'  : 'simulator_mod/dynamic_model/Dynamics/Subsystem/Multicopter_SFunction/__OutputSSForSFun__'
 * '<S45>'  : 'simulator_mod/dynamic_model/Dynamics/transfer_function/Linear Second-Order Actuator'
 * '<S46>'  : 'simulator_mod/dynamic_model/Dynamics/transfer_function/Linear Second-Order Actuator1'
 * '<S47>'  : 'simulator_mod/dynamic_model/Dynamics/transfer_function/Linear Second-Order Actuator2'
 * '<S48>'  : 'simulator_mod/dynamic_model/Dynamics/transfer_function/Linear Second-Order Actuator3'
 */
#endif                                 /* RTW_HEADER_Dynamics_h_ */
