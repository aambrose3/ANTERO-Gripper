/*
	|| ANTERO Finger Controller ||
	Author: A. B. Ambrose
	Date Created: 11/13/2024
	Date Update: 11/06/2025
	Decription: Custom library for the embedded control 
	for each ANTERO Finger. See ANTERO header file for 
	descriptions of functions.
*/
// Teensy-Arduino Libraries
#include <Arduino.h>
#include "imxrt.h" // Teensy 4.0 SoC register definitions
#include <Metro.h> // Include the Metro timing library
//Standard Open Source Libraries
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <wiring.h> // grants access to min/max/rad2deg/etc. functions
// Hardware Sprecific Libraries (Open source and customized)
#include <FlexCAN_T4.h> // FlexCAN Library for Teensy
#include <QuadEncoder.h> // Teensy 4.1 Quadrature Encoder library
#include <Dynamixel2Arduino.h>
// must include this last
#include "ANTERO_Unified.h" // This ANTERO Gripper source code

//============================
// Initialize Global Variables
//============================

// This controller will resond to CAN commands. IDs for 
// the transmitted CAN messages from the four fingers:
//uint32_t MY_CAN_ID = 0x01; // Finger 1
//uint32_t MY_CAN_ID = 0x02; // Finger 2
//uint32_t MY_CAN_ID = 0x03; // Finger 3
uint32_t MY_CAN_ID = 0x04; // Finger 4

// declare sleep mode functions -- UNTESTED --
static inline void dsb() { __asm__ volatile("dsb"); }
static inline void isb() { __asm__ volatile("isb"); }
static inline void wfi() { __asm__ volatile("wfi"); }
// set_arm_clock() is in Teensy core startup.c (C linkage) — declare it here
extern "C" uint32_t set_arm_clock(uint32_t frequency);

// Dynamixel/Actuator global variables
const uint8_t DXL_ID = 1;
const float DXL_PROTOCOL_VERSION = 2.0;
const int DXL_DIR_PIN = 22;
float DXL_OFFSET = 0;
// NOTE: Zero positions depend on how hard you run into hardstop.
// These need to be tuned again if you change the gains of harstop design.
//float DXL_ZERO = -11859; // -69.5 // need to tune until grasp force is correct (-11731 = 68.75 deg))
//float T1_ZERO = 454; // 58.4 degrees // tuned manually (468) 
//float T2_ZERO = 180; // 22.0 degrees // tuned manually (180)
float DXL_ZERO = -12459; // -69.5 degrees // -11784; // -69.1 degrees //
float T1_ZERO = 441; // 56.5 degrees //  452
float T2_ZERO = 162; // 20.2500 degrees // 180
float HOME_POSE = 0; // home position of the DXL

// Timers and state machine variabels
const float f1 = 500; // TIMER_1 frequency in Hz (must be >1)
const float f2 = 100; // TIMER_2 frequency in Hz (must be >1)
volatile uint32_t idx = 0;
volatile bool newData = false; // Flag for computing the control updates
volatile bool sendFlag = false; // Flag to send updated data over CAN
volatile bool errorFlag = false; // Flag for catching DXL communication error
volatile bool reportCard = false; // Flag to report over Serial interface

// Initialize Force Esitmation Variables
//float K[2] = {-0.014841223, 1.958607036}; // Spring quadratic stiffness coefficients (actual 1.8 mm spring)
float K[2] = {-0.044880367, 7.730617089}; // Spring quadratic stiffness coefficients (actual 2.36 mm spring)

// volatile variables tied to hardware/interrupt service routines (ISR)
volatile float t1 = 0; // angle 1 from Encoder 1 (t1)
volatile float t2 = 0; // angle 2 from encoder 2 (t2)
volatile float t3 = 0; // angle 3 from actuator (Ta)
volatile float l1 = 0; // Encoder 1 tick counter
volatile float l2 = 0; // Encoder 2 tick counter
volatile float l3 = 0; // Actuator encoder 1 tick counter
volatile bool GPIO_IFG = true; // GPIO_ISR flag
volatile bool idle = false; // init idle flag. true idles until RPi commands otherwise
volatile bool rxFlag = false; // CAN receive interrupt flag
volatile uint8_t rxData = 0; // buffer element for received CAN data;

// Initialize Control Parameters
float gain = 150; // gain for force transmission during admittance control
float MAX_FORCE = 50; // max grasp force limit (N)
float MIN_FORCE = 5; // minimum grasp force limit (N)
// For spring thickness 0.00181 → k ≈ 1.8
// For spring thickness 0.00236 → k ≈ 7.2
float k = 7.2; // approximate actual spring stiffness (N/mm) — fallback only
float k_spring = 7.2f; // instantaneous spring stiffness (N/mm) — updated every cycle by getSpring()
float Ref = 10; // grasp force preload ---------------------------------------------------
// Unified admittance controller virtual mechanical parameters.
// Transfer function: U = E * (1/c_v + 1/(J_v*s) + (k/k_v - 1)/k * s)
// Maps to PID gains:  Kp = 1/c_v,  Ki = 1/J_v,  Kd = (k/k_v - 1)/k
float c_v = 5.0f; // virtual damping     (N / raw-vel-unit)
float J_v = 10.0f;  // virtual inertia     (N·s / raw-vel-unit)
float k_v = 0.8f;  // virtual stiffness gain (dimensionless): effective stiffness = k_v*k; k_v=1 → Kd=0
float E[2] = {1, 1}; // Error term (N)
float E_dot = 0; // Derivative of Error w.r.t time (N/s)
float E_int = 0; // Integral of Error w.r.t time (Ns)
float deadBand = 0.4; // deadband (N)
float Up[2] = {0}; // position command prior
float U = 0; // actuator command
uint32_t N  = 0;

// Alpha-Beta Filter for force error derivative estimation
static const float AB_ALPHA       = 0.1f;
static const float AB_DT          = 1.0f / 500.0f; // = 1/f1
static const float AB_BETA        = (AB_ALPHA * AB_ALPHA) / (2.0f - AB_ALPHA);
static const float AB_BETA_OVR_DT = AB_BETA / AB_DT;

typedef struct {
    float x_k_hat; // filtered force error estimate
    float v_k_hat; // filtered force error derivative estimate
} AlphaBetaFilter_t;

AlphaBetaFilter_t forceFilter = {0.0f, 0.0f};

// Forward declarations (implementations at bottom of file)
static void AB_Init(AlphaBetaFilter_t *filter, float init_x, float init_v);
static void AB_Update(AlphaBetaFilter_t *filter, float measurement);

//=============================
// Initialize Global Structures
//=============================

// Declared in Header File ANTERO.h
volatile Params P{};
volatile Geometry G{};
volatile ContactForces F{};
ContactSolution sol{};
InitContact init{};
ContactParams CP{};

//==========================
// Initialize Gloabl Objects
//===========================

// extern Metro sysTimer;
extern IntervalTimer TIMER_1; // IntervalTimer TIMER_1;
extern IntervalTimer TIMER_2; // IntervalTimer TIMER_2;

//extern object for the CAN interface
extern FlexCAN_T4<CAN3, RX_SIZE_256, TX_SIZE_16> Can3;
extern CAN_message_t msg;

//extern quadrature encoder objects
extern QuadEncoder E1; // QuadEncoder E1(1, E1A, E1B, 0);
extern QuadEncoder E2; // QuadEncoder E2(2, E2A, E2B, 0);

extern Dynamixel2Arduino dxl; // Dynamixel2Arduino dxl(DXL_SERIAL, DXL_DIR_PIN)
using namespace ControlTableItem; // Here or in Main????????

//===============================================
// Finger model parameter generation (genParam.m)
//===============================================
static void getParams(void) 
{
    P.a = 81.715/1000.0; // phalange link lengths (m)
    P.b = 31.75/1000.0; // 5-bar link b length (m)
    P.d = 46.018/1000.0; // 5-bar link d length (m)
    P.e = 38.1/1000.0; // 5-bar link e length (m)
    P.t = 12.7/1000.0; // approximate thickness of the phalanges (m)

    P.alpha = 165.0 * (M_PI_F/180.0f); // angle of the link b (rad)
    P.gamma =  75.0 * (M_PI_F/180.0f); // angle of link e (rad)

    P.O_x = -47.3271/1000.0; // O: origin to palm center (m)
    P.O_y =   7.9769/1000.0;

    P.TA_max = 4.0; // max allowable actuator torque (Nm)
    P.l1_r   = 0.118342; // resting length of K1 (m)
    P.l1_o   = P.l1_r - 0.0127; // minimum allowable spring length (m)
    P.K2     = 0.005; // 0.0353; // Approximate torsion spring stiffness (Nm/rad)
    P.l2_r   = -20.0 * (M_PI_F/180.0f); // resting angle of torsion spring (rad)
    P.l2_o   = 11.0/18.0 * M_PI_F; // maximum allowable torsion spring deflection ~110 deg (rad)
}

//==========================
// Initialize the controller
//==========================
void Initialize(uint8_t state)
{
	// Initialize Generic Stuff
	delay(10); // extra boot up time
	Serial.begin(460800); // UART over USB for debugging
	//while(!Serial); // wait forever until Serial port is opened
	if (CrashReport) // Print Crash Report if needed
	{
		Serial.println("\n" __FILE__ " " __DATE__ " " __TIME__); // File Information
		Serial.println(CrashReport);
	}
	// Set up the push button pin with pull up resistor
	pinMode(9, INPUT_PULLUP);
	attachInterrupt(9, GPIO_ISR, FALLING); // LOW means pressed
	// Setup pins 14 and 15 as debug pins
	pinMode(14, OUTPUT);
	pinMode(15, OUTPUT);
	digitalWrite(14, LOW);
	digitalWrite(15, LOW);
	
	if (state == true) // Acuator is Connected
	{
		// Start Communications with the Acuator (XC330-T181-T)
		pinMode(13, OUTPUT); delayMicroseconds(10); // Actuator power mosfet gate pin. Also LED
		digitalWrite(13, HIGH); delay(500); // wait for Actuator to boot
		dxl.begin(1000000); delayMicroseconds(10); // bps of dynamixel communication
		dxl.setPortProtocolVersion(DXL_PROTOCOL_VERSION); delayMicroseconds(10);
		// Drop to 24 MHz (crystal, no PLL) while waiting for actuator boot.
		// UART clock comes from PLL3 — unaffected by set_arm_clock().
		set_arm_clock(24000000);
		delayMicroseconds(10);
		bool pingFlag = false;
		while (pingFlag == false) // Wait forever for a ping from actuator
		{
			pingFlag = dxl.ping(DXL_ID);
			if (pingFlag == true)
			{
				Serial.println("Actuator Detected");
				break;
			}
			// 5 Hz polling at 24 MHz — low power without WFI (WFI breaks USB serial on Teensy 4.x)
			delay(200);
		}
		// Restore 600 MHz before the rest of Initialize runs.
		set_arm_clock(600000000);
		delayMicroseconds(10);
		// Configuring the XC330-T181-T 
		// Never put actuator in standard position control mode -> messes up encoder readings
		dxl.torqueOff(DXL_ID); delayMicroseconds(10); // turn off actuator to unlock control registers
		dxl.ledOff(DXL_ID); delayMicroseconds(10); // turn off LED to indicate Torque is disabled
		dxl.setOperatingMode(DXL_ID, OP_EXTENDED_POSITION); delayMicroseconds(10); // put dynamixel into position mode
		dxl.writeControlTableItem(CURRENT_LIMIT, DXL_ID, 100); delayMicroseconds(10); // current limit in mA (~167 mA per Nm)
		dxl.writeControlTableItem(VELOCITY_LIMIT, DXL_ID, 66); delayMicroseconds(10); // set the maximum velocity (~4.36/66 per RPM)
		dxl.writeControlTableItem(POSITION_P_GAIN, DXL_ID, 100); delayMicroseconds(10); // set gain (should be low to start)
		dxl.writeControlTableItem(POSITION_I_GAIN, DXL_ID, 0); delayMicroseconds(10); // no integral gain
		dxl.writeControlTableItem(POSITION_D_GAIN, DXL_ID, 0); delayMicroseconds(10); // no derivative gain
		dxl.writeControlTableItem(MOVING_THRESHOLD, DXL_ID, 1); delayMicroseconds(10); // set moving threshold to 0.229/15 RPM
		dxl.ledOn(DXL_ID); delayMicroseconds(10); // turn on LED to indicate Torque is enabled
		dxl.torqueOn(DXL_ID); delayMicroseconds(10); // enable torque to lock EEPROM registers
	}
	// Set timer priorities
	// Setting these priorities caused headaches. First tried 0 and 1
	// I am guessing other libraries used interrupts use similar priorities
	// So I needed to make these timers low priority. Could test to increase in the future...
	TIMER_1.priority(250); // High Frequency Timer 
	TIMER_2.priority(251); // Low Frequency Timer
	
	//Setup the Quadrature Encoders:
	//No index pin, 
	//5 samples before trigger,
	//127 clock cycle filter bandwidth
	//pinMode(E1A, INPUT);
	//pinMode(E1B, INPUT);
	E1.setInitConfig();
	E1.EncConfig.IndexTrigger = DISABLE; // Might be disabled by default
	E1.EncConfig.INDEXTriggerMode = DISABLE; // Might be disabled by default
	//E1.EncConfig.filterCount = 0x2; // range from 0 to 7 (DEFAULT: 0 (3 samples)_
	//E1.EncConfig.filterSamplePeriod = 0x80; // range from 0 - 255 (DEFAULT: 0 (1 clock cycle))
	E1.init(); delayMicroseconds(10);
	
	E2.setInitConfig();
	E2.EncConfig.IndexTrigger = DISABLE; // Might be disabled by default
	E2.EncConfig.INDEXTriggerMode = DISABLE; // Might be disabled by default
	//E2.EncConfig.filterCount = 0x2; // range from 0 to 7 (DEFAULT: 0 (3 samples)_
	//E2.EncConfig.filterSamplePeriod = 0x80; // range from 0 - 255 (DEFAULT: 0 (1 clock cycle))
	E2.init(); delayMicroseconds(10);
	// Init the finger parameter stucture
	getParams();
	// Update the broadcast variable for the contact solver
	CP.L = P.a;
	CP.t = P.t;
	CP.Ox = P.O_x;
	CP.Oy = P.O_y;
	// Seed the initial guess for the contact solver too
	init.p1 = 0.040f;
	init.p2 = init.p1;
	init.a = init.p1;
	init.b = init.p1;
	init.valid_guess = 1;
	delay(100); // pause a little
}


//==========================
// while 1 loop with checks
//==========================
void ReturnLoop(void)
{
	while(1)
	{
		// Perform control computation if there is new data to analyze
		if (newData)
		{
			ControlUpdate(); // compute control updates
		}
		if (sendFlag)
		{
			sendDXLCommand(); // send command to actuator
		}
		// Check for GPIO interrupt flag to reset the ISR
		if (GPIO_IFG == false)
		{
			// reattach the button interrupt after debouncing
			attachInterrupt(9, GPIO_ISR, FALLING); 
			GPIO_IFG = true; // reset flag
		}
		// Transmit of Serial for debugging
		if (reportCard == true)
		{
			float time = ((float)idx)/f1;
			Serial.print(time, 4); Serial.print("\t");
			Serial.print(t1, 4); Serial.print("\t");
			Serial.print(t2, 4); Serial.print("\t");
			Serial.print(G.c, 4); Serial.print("\t");
			Serial.print(F.N1, 3); Serial.print("\t");
			Serial.print(F.N2, 3); Serial.print("\t");
			Serial.print(F.N, 3); Serial.print("\t");
			Serial. print(F.flag, 1); Serial.print("\t");
			Serial.println(Ref, 3);
			reportCard = false;
		}
		// Check for actuator hardware error -> restart state machine
		if (errorFlag == true)
		{
			TIMER_1.end();
			TIMER_2.end();
			Serial.println("Error Detected");
			// Check for ping again
			dxl.reboot(DXL_ID, 100);
			delay(20);
			bool pingFlag = false;
			while(pingFlag == false) // Wait for ping from actuator
			{
				pingFlag = dxl.ping(DXL_ID);
				delay(2);
				if (pingFlag == true)
				{
					Serial.println("Actuator Detected");
					delay(2);
				}
			}
			errorFlag = false;
			//check actuator angle and compare it to reading before reset
			// read actuator angle (4095 ppr)
			float l3temp = (float) dxl.getPresentPosition(DXL_ID, UNIT_RAW) - DXL_OFFSET; delayMicroseconds(10);
			DXL_OFFSET = DXL_OFFSET - l3 + l3temp; // update DXL_OFFSET
			// Restart state machine
			Idle();
		}
	}
}

//========================
// Initialize the actuator
//========================
void GripperInit(void)
{
	TIMER_1.end(); // stop timer interrupts
	TIMER_2.end();
	Serial.println("Homing...");
	// Put motor(s) in very weak velocity control mode...
	dxl.torqueOff(DXL_ID); delayMicroseconds(10);
	dxl.ledOff(DXL_ID); delayMicroseconds(10);
	dxl.setOperatingMode(DXL_ID, OP_VELOCITY); delayMicroseconds(10); // put dynamixel into position mode
	dxl.setGoalVelocity(DXL_ID, 0, UNIT_RAW); delayMicroseconds(10); // reset the goal velocity before turning back on
	dxl.writeControlTableItem(CURRENT_LIMIT, DXL_ID, 100); delayMicroseconds(10); // current limit in mA (~167 mA per Nm)
	dxl.writeControlTableItem(VELOCITY_LIMIT, DXL_ID, 493); delayMicroseconds(10); // set the maximum velocity (~66 per RPM)
	dxl.writeControlTableItem(PROFILE_ACCELERATION, DXL_ID, 0); delayMicroseconds(10); // set the acceleration rate (per 214 rev/min/min)
	dxl.writeControlTableItem(VELOCITY_P_GAIN, DXL_ID, 200); delayMicroseconds(10); // set gain (200 for normal)
	dxl.writeControlTableItem(VELOCITY_I_GAIN, DXL_ID, 0); delayMicroseconds(10); // set gain
	dxl.writeControlTableItem(MOVING_THRESHOLD, DXL_ID, 1); delayMicroseconds(10); // set moving threshold to 0.229/15 RPM
	dxl.ledOn(DXL_ID); delayMicroseconds(10); // turn on LED to indicate Torque is enabled
	dxl.torqueOn(DXL_ID); delayMicroseconds(10); // turn on actuator to lock control registers

	// Set velocity goal to drive system into the zeroing jig
	dxl.setGoalVelocity(DXL_ID, 120, UNIT_RAW); 
	delay(1000); // delay enought so the motor starts moving
	// move until the hardstop is hit
	volatile float vel = 10; // init velocity measurement
	DXL_OFFSET = 0;
	while (abs(vel) > 1) // hardstop has been hit: vel -> 0 RPM
	{
		vel = dxl.getPresentVelocity(DXL_ID, UNIT_RAW); delayMicroseconds(10);
		//Serial.println(vel, 0);
		delay(5); // about 200 Hz
	}
	dxl.torqueOff(DXL_ID); delayMicroseconds(10);
	dxl.ledOff(DXL_ID); delayMicroseconds(10);
	// put dynamixel into position mode
	dxl.setOperatingMode(DXL_ID, OP_EXTENDED_POSITION); delayMicroseconds(10);
	dxl.writeControlTableItem(CURRENT_LIMIT, DXL_ID, 300); delayMicroseconds(10); // current limit in mA (~1 mA)
	Serial.println("Found Home..."); delay(100);
	// set zero position (-68 deg)
	DXL_OFFSET = dxl.getPresentPosition(DXL_ID, UNIT_RAW) - DXL_ZERO; delayMicroseconds(10); 
	float test = dxl.getPresentPosition(DXL_ID, UNIT_RAW) - DXL_ZERO; delayMicroseconds(10); 
	// the homing offset value will need to be offset from the present position 
	// check the home position
	if (abs(DXL_OFFSET -test) > 10)
	{
		Serial.print("DXL_OFFSET Error\t");
		Serial.print(DXL_OFFSET, 1); Serial.print("\t");
		Serial.println(test, 1);
		DXL_OFFSET = dxl.getPresentPosition(DXL_ID, UNIT_RAW) - DXL_ZERO; delayMicroseconds(10);
		Serial.print("-->"); Serial.println(DXL_OFFSET);
	}
	// from here on out when commanding and reading position.
	// like this: dxl.getPresentPosition(DXL_ID, UNIT_RAW) - DXL_OFFSET
	E1.write(T1_ZERO); delayMicroseconds(10); // update the encoder positions 
	E2.write(T2_ZERO); delayMicroseconds(10);
	// send to the next state
	Idle();
}

//=================================================
// Send finger home and wait -> start state machine
//=================================================
void Idle(void)
{
	TIMER_1.end();
	TIMER_2.end();
	
	// send the actuator to the home/open gripper position.
	float error = -14503 - (dxl.getPresentPosition(DXL_ID, UNIT_RAW) - DXL_OFFSET);
	delayMicroseconds(10);
	if(abs(error)>10) // Send home
	{		
		Serial.println("Going Home...");
		dxl.torqueOff(DXL_ID); delayMicroseconds(10);
		dxl.ledOff(DXL_ID); delayMicroseconds(10);
		dxl.setOperatingMode(DXL_ID, OP_VELOCITY); delayMicroseconds(10); // put dynamixel into position mode
		// Clear stale goal velocity from Run mode so torqueOn doesn't chase
		// the previous positive (closing) command against loaded springs.
		dxl.setGoalVelocity(DXL_ID, 0, UNIT_RAW); delayMicroseconds(10);
		dxl.writeControlTableItem(CURRENT_LIMIT, DXL_ID, 400); delayMicroseconds(10); // current limit in mA about 1:1 (~167 mA per Nm)
		dxl.writeControlTableItem(VELOCITY_LIMIT, DXL_ID, 300); delayMicroseconds(10); // set the maximum velocity (~66 per RPM)
		dxl.writeControlTableItem(PROFILE_ACCELERATION, DXL_ID, 0); delayMicroseconds(10); // set the acceleration rate (per 214 rev/min/min)
		dxl.writeControlTableItem(VELOCITY_P_GAIN, DXL_ID, 800); delayMicroseconds(10); // set gain
		dxl.writeControlTableItem(VELOCITY_I_GAIN, DXL_ID, 80); delayMicroseconds(10); // set gain
		dxl.ledOn(DXL_ID); delayMicroseconds(10); // turn on LED to indicate Torque is enabled
		dxl.torqueOn(DXL_ID); delayMicroseconds(10); // turn on actuator to lock control registers
		dxl.setGoalVelocity(DXL_ID, -200, UNIT_RAW); delayMicroseconds(10); // start to move
		while(error<0)
		{
			// Check CAN bus for messages
			//Can3.events(); // unnecessary to update this every iteration
			if (rxFlag)
			{
				rxFlag = false;
				// Check message from master
				if (msg.buf[0] >= 127)
				{ // if dead man switch is pressed
					delayMicroseconds(10);
					idle = false;
					Run();
					ReturnLoop();
				} 
			}
			//int e1 = E1.read(); // home e1 position is 390 -> crank link ~ -85 deg
			//Serial.println(e1); // optional print out the encoder reading
			//error = -14503 - safeReadDXL(POSE, false); delayMicroseconds(10);
			error = -14503 - (dxl.getPresentPosition(DXL_ID, UNIT_RAW) - DXL_OFFSET); delay(2);
			//error = abs(error);
			//Serial.println(error);
			delayMicroseconds(10);
		}
	}
	// Update the home position in case an actuator error occurs
	HOME_POSE = dxl.getPresentPosition(DXL_ID, UNIT_RAW) - DXL_OFFSET; delayMicroseconds(10);
	// pause actuator
	dxl.torqueOff(DXL_ID); delayMicroseconds(10);
	dxl.ledOff(DXL_ID); delayMicroseconds(10);
	// Print to Console/UART
	Serial.println("System is Idle...");
	idle = true;
	// Wait for CAN command from the higher level controller
	//uint16_t idx = 0; // dummy variable for debugging
	while(idle == true)
	{
		// UNTESTED -> enter low power mode until CAN interrupt
		// Ensure all memory ops complete before sleeping
		//dsb();   // enter light sleep
		//wfi();   // CPU idles here until **any** IRQ (CAN RX) occurs
		//isb();   // Sync after wake
		
		// Can3.events(); // push CAN queue
		// Check for high level controller ID and activation code
		if (rxFlag)
		{
			CAN_message_t local_msg = msg;
			//Serial.print(local_msg.id, HEX); Serial.print(":\t");
			//Serial.println(local_msg.buf[0]);
			rxFlag = false;
			// Check message from master
			if (local_msg.buf[0] >= 127)
			{
				idle = false;
				Run();
				ReturnLoop();
			}
		} // if not the right ID or activation code -> wait again
		
		// -------------- Debugging Only --------------
		/*
		if (idx > 5000) // roughly 1000 per second
		{
			idle = false;
			Run();
			ReturnLoop();
		}
		// Counter for bypassing idle loop
		idx = idx+1;
		if (idx > 1000000) // do not roll over
		{
			idx = 0;
		}
		// Print outs to check encoders
		l1 = (float) E1.read(); // read the current encoder tick value (2880 ppr)
		l2 = (float) E2.read();
		Serial.print(l1, 0);
		Serial.print("\t");
		Serial.println(l2, 0);
		*/
		delay(5); // 200 Hz
	}
}

//======================================
// Close fingers to achieve stable grasp
//======================================
void Run(void)
{
	digitalWrite(13, HIGH);
	k_v = 1.0; c_v = 5; J_v = 10; // from simulink
	// Put XC330 into main velocity mode
	dxl.torqueOff(DXL_ID); delayMicroseconds(10);
	dxl.ledOff(DXL_ID); delayMicroseconds(10);
	dxl.setOperatingMode(DXL_ID, OP_VELOCITY); delayMicroseconds(10); // put dynamixel into position mode
	dxl.writeControlTableItem(CURRENT_LIMIT, DXL_ID, 700); delayMicroseconds(10); // current limit in mA (~167 mA per Nm)
	dxl.writeControlTableItem(VELOCITY_LIMIT, DXL_ID, 450); delayMicroseconds(10); // set the maximum velocity (~66 per RPM)
	// Can't add position limits in velocity mode -> add limitations in control update
	//dxl.writeControlTableItem(MIN_POSITION_LIMIT, DXL_ID, (int32_t)DXL_OFFSET-15400); delayMicroseconds(10); // set soft limit on position (-91 deg)
	//dxl.writeControlTableItem(MAX_POSITION_LIMIT, DXL_ID, (int32_t)DXL_OFFSET+6825); delayMicroseconds(10); // set soft limit on position (40 deg)
	dxl.writeControlTableItem(PROFILE_ACCELERATION, DXL_ID, 0); delayMicroseconds(10); // set the acceleration rate (per 214 rev/min/min)
	dxl.writeControlTableItem(VELOCITY_P_GAIN, DXL_ID, 800); delayMicroseconds(10); // set gain
	dxl.writeControlTableItem(VELOCITY_I_GAIN, DXL_ID, 5); delayMicroseconds(10); // set gain
	dxl.ledOn(DXL_ID); delayMicroseconds(10); // turn on LED to indicate Torque is enabled
	dxl.torqueOn(DXL_ID); delayMicroseconds(10); // turn on actuator to lock control registers
	AB_Init(&forceFilter, 0.0f, 0.0f); // reset alpha-beta filter before control starts
	Serial.println("Transitioning to Grasp..."); delayMicroseconds(10);
	TIMER_1.begin(TIMER_1_ISR, round((1./f1)*pow(10, 6)));
	TIMER_2.begin(TIMER_2_ISR, round((1./f2)*pow(10, 6)));
	// sit in while(1) in main
	//Serial.print(k_v, 4); Serial.print ("\t");
	//Serial.print(c_v, 4); Serial.print ("\t");
	//Serial.println(J_v, 4);
	ReturnLoop();
}


//==========================
// End state machine -> Idle
//==========================
void EStop(void)
{
	Serial.println("EMERGENCY STOP DETECTED!");
	Idle();
	return;
}

//==========================================================
// Update the force estimator and compute the control update
//==========================================================
void ControlUpdate(void)
{	
	ReadStates(); // read encoders and DXL -> updates l1, l2, l3
	//// calculate the latest force estimate
	t1 = l1/((float) E4T_PPR)*2*PI; // convert encoder count to radians
	t2 = l2/((float) E4T_PPR)*2*PI; // convert encoder count to radians
	t3 = l3/((float) XC330_PPR)*2*PI; // convert actuator count to radians
	// Solve for the contact locations and contact forces (elliptical objects)
	forceEstimator(t1, t2, t3, &F);
	// if the sum of t1 and t2 is less than 90 degrees -> object is too large -> hardcode velocity
	
	// Force Observer Flag Handler
	if (F.flag != 0)
	{
		E[1] = Ref - F.N; // destroy error history
		E[0] = E[1]; // destroy derivative
		E_int = 0; // destroy integral
		AB_Init(&forceFilter, E[1], 0.0f); // reset filter velocity
	}
	else
	{
		// Normal Control Computation (called from main after TIMER_1 executes)
		E[0] = E[1]; // update the history of the error term
		E[1] = Ref - F.N; // update the most current force error
	}
	//// Unified Admittance Controller
	// U = E*(1/c_v + 1/(J_v*s) + (k/k_v - 1)/k * s)
	// Pre-contact guard: force estimator is unreliable until both pads engage
	static bool was_in_guard = true; // start in guard at boot
	bool in_guard = (t1+t2 < 1.57);  // distal pad not yet perpendicular to palm
	if (in_guard)
	{
		E[1] = Ref + 1; // destroy error history
		E[0] = E[1];
		E_int = 0;
		E_dot = 0;
		// Filter is intentionally not updated while in guard
	}
	else
	{
		if (was_in_guard)
		{
			// First cycle out of guard: seed filter with current measurement so
			// there is no synthetic step for it to chase (which would spike E_dot
			// negative and, with Kd>0, drive U_si the wrong way).
			AB_Init(&forceFilter, E[1], 0.0f);
			E_dot = 0;
		}
		else
		{
			AB_Update(&forceFilter, E[1]);
			E_dot = forceFilter.v_k_hat;
		}
		E_int = E_int + E[1]/f1; // update approximate integral, rectangular method
		if (sgn(E[1]) != sgn(E[0]))
		{ // Integral Zero Error Cross Over Reset
			E_int = 0; // reset integral
		}
	}
	was_in_guard = in_guard;

	// global variables for debugging
	//reportFlag = F.flag;
	//reportForce = F.N;
	// NaN guard: drive finger away if force estimate is invalid
	if (isnan(E[1]))
	{
		E[1] = -1;
		E_int = 0;
		AB_Init(&forceFilter, E[1], 0.0f);
		Serial.println("NAN Force Measured");
	}
	// Compute PID gains from virtual mechanical parameters (SI units: U_si in rad/s)
	float Kp_v = 1.0f / c_v;
	float Ki_v = 1.0f / J_v;
	// Spring compression Jacobian: dcmp_mm/dt3 [mm/rad]
	// From virtual work (TA = P.d * F_cmp * cos(phi), phi = tC - t3 - pi/2):
	//   dc/dt3 = -P.d * cos(phi) = -P.d * sin(tC - t3)  [m/rad]
	//   dcmp_mm/dt3 = P.d * 1000 * |sin(tC - t3)|       [mm/rad]
	// Effective joint stiffness: k_joint = k_spring [N/mm] * J_spring [mm/rad]  [N/rad]
	float J_spring = P.d * 1000.0f * fabsf(sinf(G.tC - t3)); // [mm/rad]
	if (J_spring < 1.0f) J_spring = 1.0f;  // safety floor — avoids divide-by-zero near alignment
	float k_joint  = k_spring * J_spring;   // [N/rad]: effective joint stiffness
	float Kd_v = (1.0f / k_v - 1.0f) / k_joint; // [rad/N] ✓ — was [mm/N] before
	float U_si = Kp_v*E[1] + Ki_v*E_int + Kd_v*E_dot; // rad/s
	// Deadband
	if (abs(E[1]) <= deadBand)
	{
		U_si = 0;
	}
	// Clamp velocity command in rad/s (0.640 rad/s == 400 raw == ~61 RPM)
	if (U_si > 0.6f)
	{
		U_si = 0.6f;
	} else if (U_si < -0.6f)
	{
		U_si = -0.6f;
	}
	// Convert rad/s -> Dynamixel raw units (1 raw = 0.015266667 RPM = 0.015266667*2*pi/60 rad/s)
	const float DXL_UNITS_PER_RADS = 60.0f / (0.015266667f * 2.0f * PI); // ~625.45
	U = U_si * DXL_UNITS_PER_RADS;

	// Check to make sure the fingers joint are not hitting soft upper limits
	//if (t1 > 1.85 || t2 > 1.85 || t3 > 0.69) // 1.85 <-> 100 deg, 0.69 <-> 40 deg
	if (t3 > 0.69) // 0.69 <-> 40 deg do not really care about t1 and t2
	{
		// upper bounds on joint limits have been reached
		if (U > 0) 
		{
			// Cannot drive finger farther. Override command
			U = 0;
		}
	}

	if (t3 < -1.57 && U < 0)
	{
		U = 0;
	}
	
	newData = false; // wait back in while(1) until next interrupt records new data.
	//digitalWrite(14, LOW); // debug pins
	//digitalWrite(15, LOW);
}

//==================================================
// High Frequency Timer for Sensing -> ControlUpdate
//==================================================
void TIMER_1_ISR(void)
{
	//digitalWrite(14, HIGH); // toggle pin for debugging
	// Flag to begin computing the control updates
	idx = idx+1;
	newData = true;
}

//=====================================================================
// Lower Frequency Timer for update control/communicating with actuator
//=====================================================================
void TIMER_2_ISR(void)
{
	//digitalWrite(15, HIGH); // toggle pin for debugging
	// Report Card over Serial
	if (N > 9) // @ f/10 Hz
	{
		reportCard = true;
		N = 0;
	}
	N=N+1;
	// flag to request info be sent to the actuator
	sendFlag = true;
}

//=========================================
// Read the Encoders and Dynamixel Position
//=========================================
void ReadStates(void)
{
	delayMicroseconds(150); // minimal delay (0.15 ms) for actuator TX/RX buffer to clear
	// Measure Encoder Peripherals
	l1 = (float) E1.read(); // read the current encoder tick values (2880 ppr)
	l2 = (float) E2.read();
	
	// Check for Dynamixel Errors (blocking function call)
	int32_t dxlError = dxl.readControlTableItem(HARDWARE_ERROR_STATUS, DXL_ID); delayMicroseconds(10);		
	if (dxlError > 0)
	{
		// Report an Error over Serial
		Serial.print("Error:"); Serial.print("\t");
		Serial.print(dxlError); Serial.print("\t");
		Serial.println(l3);
		errorFlag = true;
		U = 0;
		return;
	} else
	{ // This may need to be placed outside of timer ISR
		l3 = (float) dxl.getPresentPosition(DXL_ID, UNIT_RAW) - DXL_OFFSET; delayMicroseconds(10);	// read actuator angle (4095 ppr)
	}
	
	//// check for new CAN message
	//Can3.events(); // Update CAN Queue regularly
	if (rxFlag)
	{
		uint8_t ii; // loop over packet bytes
		for (ii = 0; ii<msg.len; ii++)
		{
			rxData = msg.buf[ii];
			// Check bytes individually from master
			if (rxData < 127)
			{ // deadman switch released
				EStop();
			} else
			{ // update commands
				//Serial.println(((float) rxData-128)/127*(MAX_FORCE-5)+5, 4);
				Ref = ((float) rxData-128)/127*(MAX_FORCE-MIN_FORCE) + MIN_FORCE; // set reference grasp force
				//k_d = ((float) rxData-128)/127*(k-0.72) + 0.72; // set emulated stiffness
				//d_m = k_d/100; // set emulated damping
				//Updating the desired force and stiffness simultaneously causes issues
			}
		}
		rxFlag = false;
	}
	delayMicroseconds(150); // minimal delay (0.15 ms) for TX/RX buffer to clear
	return;
}

//=====================================================================
// Send Control Commands to the Dynamixel Actuator
//=====================================================================
void sendDXLCommand(void)
{
	delayMicroseconds(150); // minimal delay (0.15 ms) for TX/RX buffer to clear
	dxl.setGoalVelocity(DXL_ID, U, UNIT_RAW); // velocity command
	delayMicroseconds(10);
	delayMicroseconds(150); // minimal delay (0.15 ms) for TX/RX buffer to clear
}
	
//=================================================
// Interrupt Service Routine for button (debugging)
//=================================================
void GPIO_ISR(void)
{
	detachInterrupt(9);
	TIMER_1.end(); // stop timer interrupts
	TIMER_2.end();
	// Debouncing logic
	delay(5); // wait until switch chatter is over
	volatile bool switchState = digitalRead(9); // recheck the button
	while(!switchState)  // wait for button release - high
	{
		switchState = digitalRead(9);
	}
	//delay(5); // wait until switch chatter is over
	GPIO_IFG = false;
	Serial.println("Reporting Data");
	Serial.print("t1:"); Serial.print("\t"); Serial.print(t1, 4);
	Serial.print("t2:"); Serial.print("\t"); Serial.print(t2, 4);
	Serial.print("t3:"); Serial.print("\t"); Serial.print(t3, 4);
	Serial.print("p1:"); Serial.print("\t"); Serial.print(G.p1, 4);
	Serial.print("p2:"); Serial.print("\t"); Serial.print(G.p2, 4);
	Serial.print("F"); Serial.print("\t"); Serial.print(F.N, 4);
}

//=====================================
// Interrupt Service Routine for CAN RX
//=====================================
void canSniff(const CAN_message_t &RX_msg)
{
	msg = RX_msg;
	if (msg.id == 0x10) // ID for CAN control command
	{
		rxFlag = true;
	} else if (msg.id == 0x20) // ID for CAN feedback request
	{
		canTransmit();
		return;
	}
}

//================================================
// Executable for transmitting state info over CAN
//================================================
void canTransmit(void) // Need to update this...
{
	/*
	CAN_message_t msgWrite;
	msgWrite.id = MY_CAN_ID; // ID related to the finger
	uint8_t ii = 0;
	msgWrite.buf[0] = t1; // encoder 1
	msgWrite.buf[1] = t2; // encoder 2
	msgWrite.buf[2] = t3; // encoder 3
	union // decompose the float 32 bit "F" into four bytes
	{
		float f;
		uint8_t bytes[4];
	} float_union;
	float_union.f = F.N;
	for (ii=0;ii<4;ii++)
	{
		// package F into message buffer
		msgWrite.buf[3+ii] = float_union.bytes[ii]; // little-endian ordering
	}
	msgWrite.buf[7] = 0;
	if (Can3.write(msgWrite))
	{
		Serial.println("Sent Data to HLC");
	} else
	{
		Serial.println("SEND FAILED");
	}
	*/
}

//=====================
// Utility math helpers
//=====================
static float wrap2Pi(float x) 
{
    float y = fmodf(x, 2.0f*M_PI_F);
    if (y < 0.0) y += 2.0f*M_PI_F;
    return y;
}

static float clampASine(float x) 
{
    if (x >  1.0) x =  1.0;
    if (x < -1.0) x = -1.0;
    return asinf(x);
}

//==================
// 2x2 linear solver
//==================
static int32_t solve_2x2(const float A[2][2], const float b[2], float x[2]) 
{
    float det = A[0][0]*A[1][1] - A[0][1]*A[1][0];
    if (fabsf(det) < 1e-9f) return 0; // singular
    float inv00 =  A[1][1]/det;
    float inv01 = -A[0][1]/det;
    float inv10 = -A[1][0]/det;
    float inv11 =  A[0][0]/det;
    x[0] = inv00*b[0] + inv01*b[1];
    x[1] = inv10*b[0] + inv11*b[1];
    return 1;
}

//===============================================
// 4th order polynomial approximation of Cosine()
//===============================================
static float cosine(float x) 
{
  x = fmodf(x, 2*PI); // wrap angle to be in range (0, 2*pi)
  if (x<0)
  {
    x+=2*PI;
  }
  // 4th order polynomial fit of cosine in varying range (0, 2pi)
  float A[5] = {0, 0, 0, 0, 0}; // Default initialization

  if (x < 1.571f && x >= 0.0f) {
      float temp[5] = {0.02863769f, 0.02385776f, 
  -0.51558711f, 0.00366414f, 0.99980876f};
      memcpy(A, temp, sizeof(A));
  } else if (x >= 1.571f && x < 3.142f) {
      float temp[5] = {-0.02863769f, 0.38372957f, 
  -1.4051229f, 1.02232707f, 0.54800935f};
      memcpy(A, temp, sizeof(A));
  } else if (x >= 3.142f && x < 4.712f) {
      float temp[5] = {-0.02863769f, 0.33601405f, 
  -0.95541475f, -0.39780106f, 2.05051237f};
      memcpy(A, temp, sizeof(A));
  } else if (x >= 4.712f && x <= 6.284f) {
      float temp[5] = {0.02863769f, -0.74360137f, 
  6.71754479f, -24.76454405f, 31.21932936f};
      memcpy(A, temp, sizeof(A));
  }

  return A[4] + A[3] * x + A[2] * powf(x, 2) + A[1] * 
	powf(x, 3) + A[0] * powf(x, 4);
}

//=============================================
// 4th order polynomial approximation of Sine()
//=============================================
static float sine(float x) 
{
  x = fmodf(x, 2*PI); // wrap angle to be in range (0, 2*pi)
  if (x<0)
  {
    x+=2*PI;
  }
  // 4th order polynomial fit of sine in varying range (0, 2pi)
  float B[5] = {0, 0, 0, 0, 0}; // Default initialization

  if (x < 1.571f && x >= 0.0f) {
      float temp[5] = {0.02863769f, -0.20379366f, 
  0.02080391f, 0.99552651f, 0.00021991f};
      memcpy(B, temp, sizeof(B));
  } else if (x >= 1.571f && x < 3.142f) {
      float temp[5] = {0.02863769f, -0.15607814f, 
  -0.20405017f, 1.35605475f, -0.1962264f};
      memcpy(B, temp, sizeof(B));
  } else if (x >= 3.142f && x < 4.712f) {
      float temp[5] = {-0.02863769f, 0.56366547f, 
  -3.63736986f, 8.72106917f, -6.1864614f};
      memcpy(B, temp, sizeof(B));
  } else if (x >= 4.712f && x <= 6.284f) {
      float temp[5] = {-0.02863769f, 0.51594995f, 
  -2.96280763f, 5.53494123f, -1.15868086f};
      memcpy(B, temp, sizeof(B));
  }

  return B[4] + B[3] * x + B[2] * powf(x, 2) + B[1] * 
	powf(x, 3) + B[0] * powf(x, 4);
}

//=========================================
// Fitted contact locations (updateParam.m)
//=========================================
static void getContactPoints(float t1, float t2, volatile float *p1, volatile float *p2) 
{
	/*
    // These are the outdated approximate polynomials from updateParam.m (MATLAB) — see comments there.
    // p1 -> f(t1,t2)
    float p1_mm = 14.4501062925132
                 + 28.6147335594148*t1
                 +  9.54112998810529*t2
                 + (-14.3693994611985)*t1*t2
                 +   2.7779706576232*t1*t1
                 +  (-2.91159799728414)*t2*t2;

    // p2 -> f(t1,t2)
    float p2_mm =  70.5221095238136
                 + (-17.4187426901454)*t1
                 + (-28.6747238947982)*t2
                 +   16.0445700081864*t1*t2
                 +  (-17.171428241453)*t1*t1*t2
                 +   19.7847306722473*t1*t2*t2;    // Convert to meters (mm → m)
    if (p1) *p1 = p1_mm * 1e-3; // convert mm → m
    if (p2) *p2 = p2_mm * 1e-3; // convert mm → m
	*/
	// Warm-start from previous solution if available, otherwise use fixed initial guess
	if (sol.converged)
	{
		init = {.p1=sol.p1, .p2=sol.p2, .a=sol.a, .b=sol.b, .valid_guess=1};
	} else
	{
		init = {.p1=0.04f, .p2=0.04f, .a=0.04f, .b=0.04f, .valid_guess=1};
	}
	sol = ellipseContactSolve(t1, t2, &CP, &init, 30, 1e-6f, 1e-7f);
	// Bound contact locations 
	if (sol.p1 < 0.01)
	{
		sol.p1 = 0.01;
	} else if (sol.p1 > P.a*0.9)
	{
		sol.p1 = P.a*0.9;
	}
	if (sol.p2 < 0.01)
	{
		sol.p2 = 0.01;
	} else if (sol.p2 > P.a*0.9)
	{
		sol.p2 = P.a*0.9;
	}
	*p1 = sol.p1;
	*p2 = sol.p2;
}

//======================================
// compute virtual links (updateParam.m)
//======================================
static void computeKinematics(volatile Params *P, float t1, float t2, float t3, volatile Geometry *G) 
{
    // p1, p2 from fits
    getContactPoints(t1, t2, &G->p1, &G->p2);

    // Virtual links g, h (all in meters / radians)
    // g depends on actuator angle t3
    float g = sqrtf(P->e*P->e + P->d*P->d - 2.0*P->e*P->d*cosine(M_PI_F - P->gamma - t3));
    float tG = wrap2Pi( clampASine( P->d/g * sine(M_PI_F - P->gamma - t3) ) - P->gamma );

    // h depends on t2
    float h = sqrtf(P->a*P->a + P->b*P->b - 2.0*P->a*P->b*cosine(M_PI_F - P->alpha + t2));
    float tH = wrap2Pi( t1 - clampASine( P->b/h * sine(M_PI_F - P->alpha + t2) ) );

    // Spring/link c length and angle
    G->c  = sqrtf(g*g + h*h - 2.0*g*h*cosine(tH - tG));
    G->tC = M_PI_F + tG - clampASine( h/G->c * sine(tH - tG) );

    // Contact normals (pointing toward the links)
    G->n1_x = cosine(t1 - M_PI_F/2.0);  G->n1_y = sine(t1 - M_PI_F/2.0);
    G->n2_x = cosine(t1 + t2 - M_PI_F/2.0);  G->n2_y = sine(t1 + t2 - M_PI_F/2.0);

    // Direction of spring/link C
    G->nC_x = cosine(G->tC);  G->nC_y = sine(G->tC);
}

//==============================================
// Jacobians for P1, P2, and PC (getJacobians.m)
//==============================================
static void getJacobians(volatile Params *P, volatile Geometry *G, float t1, float t2,
                      float J1[2][2], float J2[2][2], float JC[2][2]) 
{
    // J1 (2x2): velocity of P1 wrt [t1;t2]
    // [-p1*sine(t1) - t*cosine(t1),   0
    //   p1*cosine(t1) - t*sine(t1),   0]
    J1[0][0] = -G->p1*sine(t1) - P->t*cosine(t1);  J1[0][1] = 0.0;
    J1[1][0] =  G->p1*cosine(t1) - P->t*sine(t1);  J1[1][1] = 0.0;

    // J2 (2x2): MATLAB-form as provided (includes +a terms in the first column)
    // J2 = [ -a*sine(t1) - p2*sine(t1+t2) - t*cosine(t1+t2),   -p2*sine(t1+t2) - t*cosine(t1+t2);
    //         a*cosine(t1) + p2*cosine(t1+t2) - t*sine(t1+t2),    p2*cosine(t1+t2) - t*sine(t1+t2) ];
    float t12 = t1 + t2;
    J2[0][0] = -P->a*sine(t1) - G->p2*sine(t12) - P->t*cosine(t12);
    J2[0][1] =                 - G->p2*sine(t12) - P->t*cosine(t12);
    J2[1][0] =  P->a*cosine(t1) +  G->p2*cosine(t12) - P->t*sine(t12);
    J2[1][1] =                   G->p2*cosine(t12) - P->t*sine(t12);

    // JC (2x2): point C at angle beta = t1 + t2 - alpha
    float zeta = t1 + t2 - P->alpha;
    JC[0][0] = -P->a*sine(t1) - P->b*sine(zeta);   JC[0][1] = -P->b*sine(zeta);
    JC[1][0] =  P->a*cosine(t1) + P->b*cosine(zeta);   JC[1][1] =  P->b*cosine(zeta);
}

//==================================
// Spring force/limits (getSpring.m)
//==================================
static int32_t getSpring(volatile Params *P, volatile Geometry *G, float t3,
                                 float *F_cmp_out, float *k_spring_out)
{
    // Check spring extension limits first
    if (G->c < P->l1_o) { // compressed beyond min length → invalid
        if (F_cmp_out)   *F_cmp_out   = 0.0f;
        if (k_spring_out) *k_spring_out = k; // fallback to approximate
        return 1;
    } else if (G->c > P->l1_r) { // longer than rest length → no spring force
        if (F_cmp_out)   *F_cmp_out   = 0.0f;
        if (k_spring_out) *k_spring_out = K[1]; // linear stiffness at zero compression
        return -1;
    }

    // Compression in millimeters (to match MATLAB fit)
    float cmp_mm = 1000.0f * (P->l1_r - G->c);

    // Nonlinear Hooke's law: F = K[0]*cmp^2 + K[1]*cmp  (Newtons)
    float F_cmp = K[0]*cmp_mm*cmp_mm + K[1]*cmp_mm;

    // Instantaneous stiffness: dF/d(cmp) = 2*K[0]*cmp + K[1]  (N/mm)
    if (k_spring_out) *k_spring_out = 2.0f*K[0]*cmp_mm + K[1];

    // Optional actuator torque cap check
    float phi = G->tC - (t3 + M_PI_F/2.0f);
    float TA  = P->d * F_cmp * cosine(phi);
    if (TA > P->TA_max) {
        if (F_cmp_out)   *F_cmp_out   = 0.0f;
        if (k_spring_out) *k_spring_out = k; // fallback to approximate
        return 1;
    }

    if (F_cmp_out) *F_cmp_out = F_cmp;
    return 0; // OK
}

//======================================================
// Public API: compute contact force magnitudes (main.m)
//======================================================
void forceEstimator(float t1, float t2, float t3, volatile ContactForces *F) 
{
    computeKinematics(&P, t1, t2, t3, &G);

    // Spring state, force, and instantaneous stiffness
    float F_cmp = 0.0f;
    int32_t sflag = getSpring(&P, &G, t3, &F_cmp, &k_spring);
    if (sflag == 1) { // invalid/out-of-range
		F->N1 = NAN;
		F->N2 = NAN;
		F->N = NAN;
		F->flag = 1;
		return;
    } else if (sflag == -1) {
        // Spring extended past rest length — no compression force, no contact.
        // Torsion spring (Tau_h) is an internal restoring moment, not a contact load;
        // continuing to solve produces a phantom N1 driven by K2*(t2-l2_r). Short-circuit.
        F->N1 = 0.0f;
        F->N2 = 0.0f;
        F->N  = 0.0f;
        F->flag = -1;
        return;
    }

    // Jacobians
    float J1[2][2], J2[2][2], JC[2][2];
    getJacobians(&P, &G, t1, t2, J1, J2, JC);

    // Build A = [J1^T*n1, J2^T*n2], b = -(Tau_h + Tau_C)
    // Column 1: J1^T * n1
    float c1_0 = J1[0][0]*G.n1_x + J1[1][0]*G.n1_y; // row=0
    float c1_1 = J1[0][1]*G.n1_x + J1[1][1]*G.n1_y; // row=1 (will be 0)

    // Column 2: J2^T * n2
    float c2_0 = J2[0][0]*G.n2_x + J2[1][0]*G.n2_y;
    float c2_1 = J2[0][1]*G.n2_x + J2[1][1]*G.n2_y;

    float A[2][2] = { { c1_0, c2_0 }, { c1_1, c2_1 } };

    // Tau_h
    float Tau_h[2] = { 0.0, -P.K2 * (t2 - P.l2_r) }; // updated the sign -> F = -kx

    // Tau_C = JC^T * nC * F_cmp
    float JCt_nC_0 = JC[0][0]*G.nC_x + JC[1][0]*G.nC_y; // row 0 of JC^T times nC
    float JCt_nC_1 = JC[0][1]*G.nC_x + JC[1][1]*G.nC_y; // row 1 of JC^T times nC
    float Tau_C[2] = { JCt_nC_0 * F_cmp, JCt_nC_1 * F_cmp };

    // b = -(Tau_h + Tau_C)
    float bvec[2] = { -(Tau_h[0] + Tau_C[0]), -(Tau_h[1] + Tau_C[1]) };

    // Solve A*[N1;N2] = b
    float x[2] = {0,0};
    if (!solve_2x2(A, bvec, x)) {
		F->N1 = NAN;
		F->N2 = NAN;
		F->N = 200;
		F->flag = 2;
		Serial.println("Singular Position Found!");
		return;
        //return (ContactForces){ .N1 = NAN, .N2 = NAN, .flag = 2 };
    }
	// Assemble Contact Force Struct
	F->N1 = x[0];
	F->N2 = x[1];
	float y = sqrt(x[0]*x[0] + x[1]*x[1]);
	F->N = y; // Compute the total enveloping force (y)
    // Discard negative contact solutions
	// Soft positive contact constraint (allows N > -1 N)
    if (x[0] < -1 || x[1] < -1) { 
		F->flag = 3;
        // flag = 3 -> Small errors in
    } else
	{
		F->flag = 0;
		// flag = 0 -> All is good
	}
	// Flag handling
	int8_t flag = F->flag;
	switch (flag){
		case -1: // Spring Extension -> check to see if the encoders are reading correctly
			F->N = -1;
			break;
		case 1: // Over Compression -> can damage the springs
			F->N = 200; // report very high force -> drives finger away from limb
			break;
		case 2: // Singular -> TO DO -> what is the protocol for this? Rely on previous info?
			F->N = 200;
			break;
		case 3: // negative contacts -> check to see if the encoders are reading correctly
			F->N = -1;
			break;
		default:
			break; // everything is good
	}
	//return (ContactForces){ .N1 = x[0], .N2 = x[1], .N = y, .flag = 0 };
}

//====================
// Helper x^2 function
//====================
static inline float sqr(float x){ return x*x; }

//============================================
// System of Equation Solver (ellipseSolver.m)
//============================================
// Solve 4x4 linear system A x = b in-place using partial pivoting.
// A is 4x4 (row-major), b is length-4. Overwrites A and b. Returns 1 on success.
static int solve4(float A[16], float b[4]) 
{
    //int piv[4] = {0,1,2,3};
    for (int k = 0; k < 4; ++k) {
        // Pivot
        int p = k;
        float amax = fabsf(A[k*4 + k]);
        for (int r = k+1; r < 4; ++r) {
            float val = fabsf(A[r*4 + k]);
            if (val > amax) { amax = val; p = r; }
        }
        if (amax < 1e-20f) return 0; // singular

        if (p != k) {
            // swap rows p and k
            for (int c = 0; c < 4; ++c) {
                float tmp = A[k*4 + c];
                A[k*4 + c] = A[p*4 + c];
                A[p*4 + c] = tmp;
            }
            float tb = b[k]; b[k] = b[p]; b[p] = tb;
        }

        float diag = A[k*4 + k];
        float invd = 1.0f / diag;
        // normalize row
        for (int c = k; c < 4; ++c) A[k*4 + c] *= invd;
        b[k] *= invd;

        // eliminate
        for (int r = 0; r < 4; ++r) {
            if (r == k) continue;
            float f = A[r*4 + k];
            if (f == 0.0f) continue;
            A[r*4 + k] = 0.0f;
            for (int c = k+1; c < 4; ++c) A[r*4 + c] -= f * A[k*4 + c];
            b[r] -= f * b[k];
        }
    }
    return 1;
}

//====================
// Residual Calculator
//====================
// Compute residual F(u) and Jacobian J(u) at u=[p1,p2,a,b]
static void residualJacobian(const ContactParams* CP,
                             float t1, float t2,
                             const float u[4],
                             float F[4], float J[16])
{
    const float p1 = fmaxf(u[0], EPS_POS);
    const float p2 = fmaxf(u[1], EPS_POS);
    const float a  = fmaxf(u[2], MIN_A_B);
    const float b  = fmaxf(u[3], MIN_A_B);

    // Precompute trig
    const float c1   = cosf(t1);
    const float s1   = sinf(t1);
    const float t12  = t1 + t2;
    const float c12  = cosf(t12);
    const float s12  = sinf(t12);

    // Points
    const float x1 = p1*c1 - CP->t * s1;            // cos(t1 - pi/2) = sin(t1)
    const float y1 = p1*s1 + CP->t * c1;            // sin(t1 - pi/2) = -cos(t1)
    const float x2 = CP->L*c1 + p2*c12 - CP->t * s12;
    const float y2 = CP->L*s1 + p2*s12 + CP->t * c12;

    const float xc = CP->Ox;
    const float yc = CP->Oy + b;

    const float dx1 = x1 - xc;
    const float dy1 = y1 - yc;
    const float dx2 = x2 - xc;
    const float dy2 = y2 - yc;

    const float inva2 = 1.0f / (a*a);
    const float invb2 = 1.0f / (b*b);

    // Residuals
    // F1: ellipse at contact 1
    F[0] = inva2*sqr(dx1) + invb2*sqr(dy1) - 1.0f;

    // F2: ellipse at contact 2
    F[1] = inva2*sqr(dx2) + invb2*sqr(dy2) - 1.0f;

    // F3: tangency at contact 1
    F[2] = (dx1*inva2)*c1 + (dy1*invb2)*s1;

    // F4: tangency at contact 2
    F[3] = (dx2*inva2)*c12 + (dy2*invb2)*s12;

    // Jacobian entries (row-major, J[i*4 + j] = dFi/duj)
    // Helpers
    const float dxdp1_1 = c1,  dydp1_1 = s1;
    const float dxdp2_2 = c12, dydp2_2 = s12;

    // dF1/dp1
    J[0*4 + 0] = 2.0f*inva2*dx1*dxdp1_1 + 2.0f*invb2*dy1*dydp1_1;
    // dF1/dp2
    J[0*4 + 1] = 0.0f;
    // dF1/da
    J[0*4 + 2] = -2.0f * sqr(dx1) / (a*a*a);
    // dF1/db: via invb2 and yc
    J[0*4 + 3] = (-2.0f * sqr(dy1) / (b*b*b)) - (2.0f * invb2 * dy1);

    // dF2/dp1
    J[1*4 + 0] = 0.0f;
    // dF2/dp2
    J[1*4 + 1] = 2.0f*inva2*dx2*dxdp2_2 + 2.0f*invb2*dy2*dydp2_2;
    // dF2/da
    J[1*4 + 2] = -2.0f * sqr(dx2) / (a*a*a);
    // dF2/db
    J[1*4 + 3] = (-2.0f * sqr(dy2) / (b*b*b)) - (2.0f * invb2 * dy2);

    // F3
    // dF3/dp1
    J[2*4 + 0] = inva2 * c1 * c1 + invb2 * s1 * s1;
    // dF3/dp2
    J[2*4 + 1] = 0.0f;
    // dF3/da
    J[2*4 + 2] = -2.0f * c1 * dx1 / (a*a*a);
    // dF3/db  (through invb2 and dy1)
    J[2*4 + 3] = (-2.0f * s1 * dy1 / (b*b*b)) - (s1 * invb2);

    // F4
    // dF4/dp1
    J[3*4 + 0] = 0.0f;
    // dF4/dp2
    J[3*4 + 1] = inva2 * c12 * c12 + invb2 * s12 * s12;
    // dF4/da
    J[3*4 + 2] = -2.0f * c12 * dx2 / (a*a*a);
    // dF4/db
    J[3*4 + 3] = (-2.0f * s12 * dy2 / (b*b*b)) - (s12 * invb2);
}

//=====================
// Helper norm function
//=====================
static float norm2_4(const float v[4])
{
    return sqrtf(v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + v[3]*v[3]);
}

//===============================
// Helper compare matrix function
//===============================
int cmpfunc(const void *a, const void *b)
{
	float fa = *(const float *)a;
    float fb = *(const float *)b;
    return (fa > fb) - (fa < fb);
}

//============================================================================
// Public API: Compute the contact locations assuming the object is an ellipse
//============================================================================
ContactSolution ellipseContactSolve(float t1, float t2,
                              const ContactParams* CP,
                              const InitContact* init,
                              int32_t max_iters,
                              float f_tol,
                              float step_tol)
{
    ContactSolution out;
    out.p1 = out.p2 = out.a = out.b = 0.0f;
    out.iterations = 0;
    out.residual_norm = INFINITY;
    out.converged = 0;

    // Initial guess
    float u[4];
    if (init && init->valid_guess) {
        u[0] = fmaxf(init->p1, EPS_POS);
        u[1] = fmaxf(init->p2, EPS_POS);
        u[2] = fmaxf(init->a,  MIN_A_B);
        u[3] = fmaxf(init->b,  MIN_A_B);
    } else {
        // Heuristic: moderate positive values scaled to geometry
        const float scale = fmaxf(0.5f*CP->L + 2.0f*CP->t, 0.02f);
        u[0] = fmaxf(0.5f*CP->L, EPS_POS);
        u[1] = fmaxf(0.5f*CP->L, EPS_POS);
        u[2] = fmaxf(0.8f*scale, MIN_A_B);
        u[3] = fmaxf(0.6f*scale, MIN_A_B);
    }

    for (int32_t k = 0; k < max_iters; ++k) {
        float F[4], J[16];
        residualJacobian(CP, t1, t2, u, F, J);
        float Fn = norm2_4(F);
        out.iterations = k + 1;

        if (Fn < f_tol) {
            out.converged = 1;
            out.residual_norm = Fn;
            break;
        }

        // Solve J * delta = -F
        float A[16];
        float rhs[4] = { -F[0], -F[1], -F[2], -F[3] };
        for (int i = 0; i < 16; ++i) A[i] = J[i];
        if (!solve4(A, rhs)) {
            // Jacobian singular: bail
            out.converged = 0;
            out.residual_norm = Fn;
            break;
        }

        // Backtracking line search
        float alpha = 1.0f;
        float stepn = norm2_4(rhs);
        if (stepn < step_tol) {
            out.converged = 1;
            out.residual_norm = Fn;
            break;
        }

        float u_try[4], F_try[4];
        for (int bt = 0; bt < MAX_BACKTRACK; ++bt) {
            u_try[0] = fmaxf(u[0] + alpha*rhs[0], EPS_POS);
            u_try[1] = fmaxf(u[1] + alpha*rhs[1], EPS_POS);
            u_try[2] = fmaxf(u[2] + alpha*rhs[2], MIN_A_B);
            u_try[3] = fmaxf(u[3] + alpha*rhs[3], MIN_A_B);

            float Jtmp[16];
            residualJacobian(CP, t1, t2, u_try, F_try, Jtmp);
            float Fn_try = norm2_4(F_try);

            if (Fn_try <= (1.0f - BACKTRACK_C*alpha) * Fn) {
                // sufficient decrease
                u[0]=u_try[0]; u[1]=u_try[1]; u[2]=u_try[2]; u[3]=u_try[3];
                break;
            }
            alpha *= BACKTRACK_BETA;
            if (alpha < 1.0e-6f) {
                // take tiny step
                u[0]=u_try[0]; u[1]=u_try[1]; u[2]=u_try[2]; u[3]=u_try[3];
                break;
            }
        }
        // continue iterating
        if (k+1 == max_iters) {
            float Fend[4], Jend[16];
            residualJacobian(CP, t1, t2, u, Fend, Jend);
            out.residual_norm = norm2_4(Fend);
        }
    }

    out.p1 = u[0];
    out.p2 = u[1];
    out.a  = u[2];
    out.b  = u[3];
    return out;
}

//===================================================
// Alpha-Beta Filter: Initialize
//===================================================
static void AB_Init(AlphaBetaFilter_t *filter, float init_x, float init_v)
{
    filter->x_k_hat = init_x;
    filter->v_k_hat = init_v;
}

//===================================================
// Alpha-Beta Filter: Update (call every AB_DT = 1/f1)
//===================================================
static void AB_Update(AlphaBetaFilter_t *filter, float measurement)
{
    float x_pred   = filter->x_k_hat + filter->v_k_hat * AB_DT;
    float v_pred   = filter->v_k_hat;
    float residual = measurement - x_pred;
    filter->x_k_hat = x_pred + AB_ALPHA * residual;
    filter->v_k_hat = v_pred + AB_BETA_OVR_DT * residual;
}
