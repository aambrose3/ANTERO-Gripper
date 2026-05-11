#pragma once
/* ||| ANTERO Finger Controller |||
	Author: A. B. Ambrose
	Date Created: 11/13/2024
	Date Update: 11/06/2025
	Decription: Custom library for controlling the ANTERO Finger controllers
	
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
===============================================
*/

#ifndef ANTERO_Unified_h
#define ANTERO__Unified_h

#include <stdint.h>

// Macros
#define sgn(x) ((x) < 0 ? -1 : ((x) > 0 ? 1 : 0))
#define E1A 0
#define E1B 1
#define E2A 2
#define E2B 3
#define DXL_SERIAL Serial5
#define ZERO -180
#define E4T_PPR 2880
#define XC330_PPR 61439

#ifndef FABS
#define FABS(x) ((x) < 0 ? -(x) : (x))
#endif

#ifndef M_PI_F
#define M_PI_F 3.14159265358979323846f
#endif

// Safety floors for contact solver to keep denominators reasonable
#define EPS_POS      (1.0e-6f)
#define MIN_A_B      (1.0e-5f)
#define MIN_STEP     (1.0e-12f)
#define BACKTRACK_BETA 0.5f
#define BACKTRACK_C    1.0e-4f
#define MAX_BACKTRACK  10

#define FILTER_LENGTH 7
#define POSE 0
#define VEL 1
#define CURR 2

// ID for transmitted CAN messages from this controller 
// Each finger should be different so comment out unused IDs
extern uint32_t MY_CAN_ID; // Finger ID

// Dynamixel/Actuator global variables
extern const uint8_t DXL_ID;
extern const float DXL_PROTOCOL_VERSION;
extern const int DXL_DIR_PIN;
extern float DXL_OFFSET;
extern float DXL_ZERO; // -67.96 degrees
extern float T1_ZERO; // 57.9 degrees
extern float T2_ZERO; // 20 degrees
extern float HOME_POSE; // home position of the DXL

// Timers and state machine variabels
extern const float f1; // TIMER_1 frequency in Hz (must be >1)
extern const float f2; // TIMER_2 frequency in Hz (must be >1)
extern volatile uint32_t idx;
extern volatile bool newData; // Flag for computing the control updates
extern volatile bool sendFlag; // Flag to send commands to dynamixel 
extern volatile bool errorFlag; // Flag for dynamixel error detection
extern volatile bool reportCard; // flag for routinely reporting data over serial

// Initialize Force Esitmation Variables
extern float K[2]; // Spring quadratic stiffness coefficients

// volatile variables tied to hardware/interrupt service routines
extern volatile float t1; // angle 1 form Encoder 1 (t1)
extern volatile float t2; // angle 2 from encoder 2 (t2)
extern volatile float t3; // angle 3 from actuator (Ta)
extern volatile float l1; // Encoder 1 tick counter
extern volatile float l2 ; // Encoder 2 tick counter
extern volatile float l3; // Actuator encoder 1 tick counter
extern volatile bool GPIO_IFG; // GPIO_ISR flag
extern volatile bool idle; // init idle flag. true idles until RPi commands otherwise
extern volatile bool rxFlag; // CAN receive interrupt flag
extern volatile uint8_t rxData; // buffer element for received data;

// Initialize Control Parameters
extern float gain; // gain for force transmission during maintain grasp
extern float MAX_FORCE; // max grasp force limit (N)
extern float k; // approximate actual spring stiffness (N/mm) — fallback only
extern float k_spring; // instantaneous spring stiffness (N/mm) — updated each cycle by getSpring()
extern float Ref; // grasp force preload ----------------------
// Unified admittance controller virtual mechanical parameters
// U = E*(1/c_v + 1/(J_v*s) + (k/k_v - 1)/k * s)
extern float c_v; // virtual damping     (N / raw-vel-unit) → Kp = 1/c_v
extern float J_v; // virtual inertia     (N·s / raw-vel-unit) → Ki = 1/J_v
extern float k_v; // virtual stiffness gain (dimensionless): effective stiffness = k_v*k → Kd = (1/k_v - 1)/k
extern float E[2]; // Error term (N)
extern float E_dot; // Derivative of Error w.r.t time (N/s)
extern float E_int; // Integral of Error w.r.t time (Ns)
extern float beta; // Tunable Alhpa Filter parameter (alpha already in use). 
extern float deadBand; // deadband (N)
extern float Up[2]; // position command prior
extern float U; // actuator command
extern uint32_t N;

//=============================
// Declare Global Structures
//=============================

//=====================================
// Output structure from forceEstimator
//=====================================
typedef struct {
    float   N1;   // contact force magnitude on proximal pad [N]
    float   N2;   // contact force magnitude on distal pad   [N]
	float   N;    // total enveloping grasp force
    int32_t flag; // status code (see FINGER_* codes below)
} ContactForces;

//===========================
// Parameter block (getParam)
//===========================
typedef struct {
    // Fixed link lengths (meters)
    float a;     // 0.081715
    float b;     // 0.03175
    float d;     // 0.046018
    float e;     // 0.0381
    float t;     // 0.0127 (pad offset thickness)

    // Fixed angles (radians)
    float alpha; // deg2rad(165)
    float gamma; // deg2rad(75)

    // Palm center (unused in core solve)
    float O_x;   // -47.3271/1000
    float O_y;   // 7.9769/1000

    // Constraints / springs
    float TA_max;  // 6 Nm (cap for actuator torque check)
    float l1_r;    // 0.118342 m (compression spring rest length)
    float l1_o;    // minimal compression length (rest - 0.015 m)
    float K2;      // 0.03 Nm/rad (torsion spring about joint 2)
    float l2_r;    // -deg2rad(20) (torsion spring rest angle)
    float l2_o;    // ~110 deg (unused here but kept for parity)
} Params;

//=============================
// Geometry (updateParam parts)
//=============================
typedef struct {
    // derived geometry
    float p1, p2;       // meters (after scale)
    float tC;           // global angle of spring/link c
    float c;            // length of spring/link c
    float n1_x, n1_y;   // contact normal at pad 1 (toward link)
    float n2_x, n2_y;   // contact normal at pad 2 (toward link)
    float nC_x, nC_y;   // unit vector along link C (from O4→O3)
} Geometry;

typedef struct {
    float L;   // link length (m)
    float t;   // offset t (m)
    float Ox;  // O(1) (m)
    float Oy;  // O(2) (m)
} ContactParams;

typedef struct {
    float p1, p2, a, b;   // unknowns to solve
    int32_t iterations;   // iterations used
    float residual_norm;  // ||F||_2 at exit
    int32_t converged;    // 1 on success, 0 otherwise
} ContactSolution;

typedef struct {
    float p1, p2, a, b; // optional initial guess (all > 0)
    int valid_guess;       // set to 1 if fields are valid, else 0
} InitContact;


extern volatile Params P;
extern volatile Geometry G;
extern volatile ContactForces F;
extern ContactSolution sol;
extern InitContact init;
extern ContactParams CP;



//// Function Declarations

// Setup the Teensy 4.0 (pins, communication protocols, etc.)
extern void Initialize(uint8_t card_presence);
// While 1 loop where the code should default back to
extern void ReturnLoop(void);
// Begin the homing process of the fingers on boot up
extern void GripperInit(void);
// Function to start the PID velocity-force control of the fingers
extern void Run(void);
// Stops updating the control of the gripper and send the gripper home
extern void Idle(void);
// Sends gripper fingers to home position ASAP
extern void EStop(void);
// Function for updating the control parameters
extern void ControlUpdate(void);
// High-frequency timer interrupt service routine - Sensing
extern void TIMER_1_ISR(void);
// LOW-frequency timer interrupt service routine - Commanding
extern void TIMER_2_ISR(void);
// Read the encoders and the dynamixel position (flag from TIMER_1_ISR)
extern void ReadStates(void);
// Send the updated dynamixel commands (flag from TIMER_2_ISR)
extern void sendDXLCommand(void);
// GPIO interrupt service routine (button pressed)
extern void GPIO_ISR(void);
// CAN receive interrupt service routine
extern void canSniff(const CAN_message_t &msg);
// CAN transmit function for feedback, if requested
extern void canTransmit(void);
// reads and checks the position (ID=0), velocity (ID=1), or current (ID=2) readings from the dyanmixel
extern float safeReadDXL(uint8_t ID);
// used to find the median of an array -> sorts an array from lowest to highest
extern int cmpfunc(const void *a, const void *b);

//======================================================
// Public API: compute contact force magnitudes (main.m)
//======================================================
extern void forceEstimator(float t1, float t2, float t3, volatile ContactForces *F);

extern ContactSolution ellipseContactSolve(float t1, float t2,
                              const ContactParams* P,
                              const InitContact* init,
                              int32_t max_iters,
                              float f_tol,
                              float step_tol);

#endif // ANTERO_Unified_h
