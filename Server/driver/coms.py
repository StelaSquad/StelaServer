import serial
import time
import sys

def testFunc():
    print("called testFunc", file=sys.stderr)

def serial():
    ser = serial.Serial('/dev/ttyS8')
    c=0

    while c<=5:
	    ser.write('2')
	    time.sleep(1)
	    ser.write('1')
	    time.sleep(1)
	    c += 1


