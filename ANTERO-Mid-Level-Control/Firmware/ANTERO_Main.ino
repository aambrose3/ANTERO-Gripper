#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <wiring.h> // grants access to min/max/rad2deg/etc. functions
#include <Metro.h> // Include the Metro library
#include <FlexCAN_T4.h>
#include <QuadEncoder.h>
#include <Dynamixel2Arduino.h>
// must include this last
#include "ANTERO_Unified.h" // Custom Library for the ANTERO Gripper controller source code
// #include "ANTERO.h" // Old library with 2 controllers

//Timer objects
IntervalTimer TIMER_1;
IntervalTimer TIMER_2;

// Encoder 1 object
QuadEncoder E1(1, E1B, E1A, 0);
// Encoder 2 object (See notes)
/* 
	|::::::|                    |::::::|
	|::::::|                    |::::::|
	|::::::|                    |::::::|
	|::::::|<--4            1-->|::::::|
	|::::::|                    |::::::|
	|::::::|      Top View      |::::::|
	|::::::: :::::::: :::::::: ::::::::|
	|::::::: :::::::: :::::::: ::::::::|
	|::::::: :::::::: :::::::: ::::::::|
	|::::::: :::::::: :::::::: ::::::::|     
	         |::::::| |::::::|
	         |::::::| |::::::|
	     3-->|::::::| |::::::|<--2
	         |::::::| |::::::|
	         |::::::| |::::::|
	         |::::::| |::::::|	
*/
// Quadrature encoder object //////// ONLY FINGERS 1 & 2 ////////
QuadEncoder E2(2, E2A, E2B, 0); // comment out for programming fingers 3 & 4

// Quadrature encoder object //////// ONLY FINGERS 3 & 4 ////////
//QuadEncoder E2(2, E2B, E2A, 0); // comment out for programming fingers 1 & 2

// Dynamixel actuator object
Dynamixel2Arduino dxl(DXL_SERIAL, DXL_DIR_PIN);

// CAN controller object
FlexCAN_T4<CAN3, RX_SIZE_256, TX_SIZE_16> Can3;
CAN_message_t msg;

void setup() 
{
	Initialize(true); // actuator is connected
  delay(1);
  // The CAN setup needs to happen in this main file, not source code.
  //Start CAN tranceiver
	Can3.begin(); delay(1);
	pinMode(29, OUTPUT); digitalWrite(29, HIGH); delay(1);
	Can3.setBaudRate(250000); delay(1);
	digitalWrite(29, LOW); delay(1);
	// Start CAN ISR
	Can3.setMaxMB(16); delay(1);
	Can3.enableFIFO(); delay(1);
	Can3.enableFIFOInterrupt(); delay(1);
	Can3.onReceive(canSniff); delay(1); // CAN ISR is in the RPTG.cpp source code
	Can3.mailboxStatus(); delay(50);  // no CAN messages can be received before any of this setup 
  delay(5);
}

void loop()
{
	GripperInit();
	ReturnLoop();
	Serial.println("Wrong Place Wrong Time!");
}
