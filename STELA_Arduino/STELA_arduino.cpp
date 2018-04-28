#include "Arduino.h"
#include "STELA_arduino.h"
#include "TimerOne.h"

STELA_arduino::STELA_arduino(int in) : AzMotor("az",1,2,3,4,5), AltMotorA("alt",5,6,7,8,9), AltMotorB("alt",1,2,3,4,5) {
    Serial.begin(9600);
    Timer1.initialize(1000);
}

// Master code for reading the serial port information
void STELA_arduino::SerialUpdate(){

    // Read message
    STELA_arduino::read();

    // Check messages for commands, if one is recognized, respond somehow
    if (inmsg == "info"){
        // Request for position and target update.
        STELA_arduino::send_update();

    } else if (inmsg.substring(0,8) == "set_targ"){
        // Input format: "set_targ:[{azimuth},{altitude}]"
        STELA_arduino::set_targ();
    
    } else if (inmsg.substring(0,7) == "set_pos"){

        STELA_arduino::set_pos();

    } else if (inmsg == "setup"){
        // Returns message with name
        Serial.println("STELA Microcontroller 1.0");

    } else if (inmsg.length() > 0){
        // If it doesn't match any protocols, return error
        Serial.print("Command not recognized: ");
        Serial.println(inmsg);
    }

}

void STELA_arduino::Move(){
    AzMotor.targ = targ[0];
    AzMotor.move();
    pos[0] = AzMotor.pos;

    //Serial.println(targ[0])

    AltMotorA.targ = targ[1];
    AltMotorB.targ = targ[1];
    AltMotorA.move();
    AltMotorB.move();
    pos[1] = AltMotorA.pos;
}

void STELA_arduino::read() // Extracts the incoming serial message
{   
    inmsg = "";
    while (Serial.available()){
        inmsg += char(Serial.read());
    }
}

void STELA_arduino::send_update() // Sends position and target info through serial port
{
    outmsg = ("Pos: [" + String(pos[0],6) + ", " + String(pos[1],6) + "]"  + 
              "\tTarg: [" + String(targ[0],6)+", " + String(targ[1],6)+"]" +
              "\tGear: [" + String(AzMotor.gear) + ", " + String(AzMotor.gear) + "]");
    Serial.println(outmsg);
}

void STELA_arduino::set_targ() // Command: "set_targ:[{az},{alt}]" to set targets
{
    msginds[0] = inmsg.indexOf("[")+1;
    msginds[1] = inmsg.indexOf(",");
    msginds[2] = inmsg.indexOf("]");

    targ[0] = inmsg.substring(msginds[0],msginds[1]).toFloat();
    targ[1] = inmsg.substring(msginds[1]+1,msginds[2]).toFloat();

}

void STELA_arduino::set_pos()
{
    msginds[0] = inmsg.indexOf("[")+1;
    msginds[1] = inmsg.indexOf(",");
    msginds[2] = inmsg.indexOf("]");

    pos[0] = inmsg.substring(msginds[0],msginds[1]).toFloat();
    pos[1] = inmsg.substring(msginds[1]+1,msginds[2]).toFloat();

    targ[0] = pos[0];
    targ[1] = pos[1];

}

void STELA_arduino::set_targ_man(float az,float alt)
{
    targ[0] = az;
    targ[1] = alt;
}





/////////////////////////////////////////
///////   Motor Control Class  //////////
/////////////////////////////////////////

motor::motor(String type, int set_dir_pin, int set_step_pin, int set_MS1_pin, int set_MS2_pin, int set_MS3_pin)
{
    // Initialize and set class variables for the needed pins
    dir_pin = set_dir_pin;
    step_pin = set_step_pin;
    MS1_pin = set_MS1_pin;
    MS2_pin = set_MS2_pin;
    MS3_pin = set_MS3_pin;

    pinMode(dir_pin,OUTPUT);
    pinMode(step_pin,OUTPUT);
    pinMode(MS1_pin,OUTPUT);
    pinMode(MS2_pin,OUTPUT);
    pinMode(MS3_pin,OUTPUT);
}

void motor::move() // Calculates when to step the stepper
{
    looptime = (micros()-lastloop)/pow(10,6);
    targspeed = (targ-lasttarg)/looptime;
    lasttarg = targ;

    //Serial.println(targ);
    // This part allows the telescope to find the shortest path (around 360->0)
    mid = proportion*(targ-pos);
    hi = mid + proportion*360;
    lo = mid - proportion*360;

    // Choose the fastest path
    if (abs(mid)<abs(hi) and abs(mid)<abs(lo)){
        newspd = mid;
    } else if (abs(hi)<abs(lo)){
        newspd = hi;
    } else{
        newspd = lo;
    }

    spd = newspd + 0.1*(targspeed - spd);

    // Calculate the looptime and the time since the last step
    //looptime = micros()-lastrun;
    //lastrun = micros();

    //for 200 steps per revolution (1.8deg/stp), 1step/2ms = 500 stps/sec maximum
    // 1 -> 1.8 deg/stp         max -> 
    // 2 -> 0.9 deg/stp
    // 4 -> 0.45 deg/stp
    // 8 -> 0.225 deg/stp
    // 16 -> 0.1125 deg/stp

    maxstepsize = 1.8/gearRatio; //1.63;
    gear = 5;
    gearval = 16;
    stepsize = maxstepsize/gearval;
    maxspd = maxstepspeed*stepsize;

    
    while( (abs(spd) > maxspd) && (gear > 1)){
        gear--;
        gearval = pow(2,int(gear-1));
        stepsize = maxstepsize/gearval;
        maxspd = maxstepspeed*stepsize;
    }

    // Lower speed if it's over the limit
    if (spd > maxspd){
        spd = maxspd;
    }
    if (spd < -maxspd){
        spd = -maxspd;
    }

    if (abs(targ-pos) < stepsize){
        spd = 0;
    }


    //Serial.println(azspd);
    
    if ((micros()-laststeptime) / pow(10,6) * abs(spd) > stepsize){
        motor::set_speed(gear);

        //Serial.println("STEP!");
        motor::step(int(spd/abs(spd)));

        pos += stepsize * spd/abs(spd);

        if (pos > 360){
            pos -= 360;
        } else if (pos < 0){
            pos += 360;
        }

        laststeptime = micros();
    }

}

void motor::step(int dir)
{
    // Steps the motor once
    digitalWrite(dir_pin,(dir*azplus+1)/2);
    
    digitalWrite(step_pin,HIGH);
    digitalWrite(step_pin,LOW);

    digitalWrite(dir_pin,0);
}

void motor::set_speed(int gear)
{
    // Changes the speed
    // Gear variable:
    //  1 -> full step
    //  2 -> half step
    //  3 -> quarter step
    //  4 -> eigth step
    //  5 -> ninth step

    gear--;

    digitalWrite(MS1_pin, gearspeeds[gear][0]);
    digitalWrite(MS2_pin, gearspeeds[gear][1]);
    digitalWrite(MS3_pin, gearspeeds[gear][2]);

}
