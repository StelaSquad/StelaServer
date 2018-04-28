
#include "STELA_arduino.h"
#include <Time.h>


STELA_arduino stela(2);
motor AzMotor("az",13,A0,A3,A2,A1);
motor AltMotorA("alt",12,11,8,9,10);
motor AltMotorB("alt",5,6,2,3,4);


void setup() {
  // Setup arduino, specify motors 

  Serial.begin(9600);
  delay(20);

  stela.AzMotor = AzMotor;
  stela.AltMotorA = AltMotorA;
  stela.AltMotorB = AltMotorB;

}

void loop() {
  // Keep updating and responding accordingly
  stela.SerialUpdate();
  stela.Move();

}
