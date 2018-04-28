/*
    STELA_arduino.h - Library for STELA arduino operations
    Created by Duncan Rocha for the STELA team, March 22, 2018
*/
#ifndef STELA_arduino_h
#define STELA_arduino_h

#include "Arduino.h"
#include "TimerOne.h"

class motor
{
    public:
        // Functions
        motor(String type, int dir_pin, int step_pin,int MS1_pin, int MS2_pin, int MS3_pin);
        void move();
        void step(int dir);
        void set_speed(int gear);
        //static void stepup(int dir);

        // Variables
        int dir_pin,step_pin, MS1_pin,MS2_pin,MS3_pin;
        String type;
        int azplus=1;

        double targ;
        double pos, hi, mid, lo, spd, newspd;
        double proportion = 20;
        double stepsize,maxstepsize,maxspd;
        float maxstepspeed = 1/0.004;

        double laststeptime = 0;

        float gearRatio = 4;
        int gear;
        float gearval;

    private:

        double looptime,targspeed,lasttarg,lastloop;
        // Speed settings (works for 1,2,4,8,16 gearspeed driver)
        int gearspeeds[5][3] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{1,1,1}};

};


class STELA_arduino
{

    public:
        STELA_arduino(int in);

        motor AzMotor;
        motor AltMotorA;
        motor AltMotorB;

        void SerialUpdate();
        void set_targ_man(float az, float alt);
        void Move();
        
        String inmsg;
        String outmsg;
        double pos[2] = {0,0};
        double targ[2] = {0,0};

        void send_update();
        
    private:
        void set_targ();
        void set_pos();
        void read();
        char inchar;
        int msginds[3] = {0,0,0};


};

#endif
