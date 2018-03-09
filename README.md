# StelaServer
Small Server meant to run on a raspberry pi and handle requests that will correspond to commands for the Stela Project Telescope.

## Installation
 ###Ngrok
 We currently use Ngrok. This service sets up a secure tunnel between the server and a public url endpoint. 
 Download here: https://ngrok.com/download 
 This will give a command line application. Move to its directory and 
 run 
 	./ngrok http 5000
 This will setup a tunnel through port 5000 (which is the port our server uses by default).

 ###Python dependencies 
 	The following commnads should install the required python
 	dependicies in your system. 
 	
 	*pip install Flask
 	*pip intall flask-autodoc
	*pip install Serial

## Testing
###Running Unit Tests
	The file contains (or will contain soon) unit-tests for all endpoints. 
	* python index_test.py

##
Tested on python 3.6.2 and python 3.6.4
