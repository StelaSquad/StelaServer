
#include "STELA_arduino2.h"
#include <Time.h>
//#include <TimerThree.h>


STELA_arduino stela(2);
motor AzMotor("az",13,A0,A3,A2,A1);
motor AltMotorA("alt",12,11,8,9,10);
motor AltMotorB("alt",5,6,2,3,4);


void setup() {
  // put your setup code here, to run once:
  Serial.begin(9600);
  delay(50);

  stela.AzMotor = AzMotor;
  stela.AltMotorA = AltMotorA;
  stela.AltMotorB = AltMotorB;
}

double lastruntime;
double dt;
float spd = 1;
double targaz;
int dir = 1;
float max_spd = 10;

void loop() {
  // set_targ:[180,0]
  
  dt = (micros()-lastruntime)/pow(10,6);
  lastruntime = micros();

  spd = 0.00416666;
  targaz += spd*dt;
  if (targaz>360){
    targaz-=360;
    //Serial.println("full rotation!");  
  } else if(targaz<0){
    targaz += 360;
    //Serial.println("full rotation!");
  }
  spd += dir * max_spd*dt/10;
  
  if(abs(spd)> max_spd){
    spd *= -1/abs(spd);
    dir *= -1;
  }
  
  stela.send_update();
  
  stela.set_targ_man(targaz,targaz);
  stela.SerialUpdate();
  stela.Move();

}
