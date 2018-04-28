# StelaServer
Small Server meant to run on a raspberry pi and handle requests that will correspond to commands for the Stela Project Telescope.

## Installation
Currently the code is set to run on either a laptop or a raspberry pi. Generally on a laptop, the usb ports must be unlocked by changing the permissions of /dev/ttyS* to allow access from the script. 

To comunicate with the server, one can use the python script "telecon.py" in the highest folder. The right ip code must be specified, the ip of the machine running the server. 
Generally, if the request just times out or hangs, the firewall or the network is preventing the signal from going through. It's usually better to just try and run it on a pi or something with no security measures. 

The Arduino scripts are located in the Arduino folder. There are two files that house the classes used, and one that runs the thing.


 ###Python dependencies 
 	The following commnads should install the required python
 	dependicies in your system. 
 	
 	*pip install Flask
 	*pip intall flask-autodoc
	*pip install Serial
	*pip install numpy
	*pip install astropy
	*pip install astroquery
	*pip install pyzmq

## Testing
###Running Unit Tests
	To test the motors, upload STELA_test.ino to the arduino, but first make sure that the pins are properly specified. The motors should start at a slow speed, gradually speed up, and then do the same process in reverse after 10 seconds.

	There is currently no test for STELA, nor the server. The program must be debugged manually. 

##
Tested on python 2.7
