from flask import Flask, jsonify, request
from flask import render_template
from flask_autodoc import Autodoc

# TODO: move these imports along with the 
# serail coms functionality to another file
import serial
import time
import STELA as st
import zmq
import json

port = "5122"
context = zmq.Context()
socket = context.socket(zmq.PAIR)
socket.connect("tcp://localhost:%s" % port)
socket.SNDTIMEO = 1000

app = Flask(__name__)
auto = Autodoc(app)
app.config['DEBUG'] = False

st.count = 0
print "Setup Complete"

calib_coors=[[],[],[]]


# Setup code. Currently only gets reference stars and returns them
@app.route('/setup',methods=['POST'])
def setup():
    data = request.get_json(force=True)
    data = json.loads(request.data)
    data["command"] = "setup"
    socket.send(json.dumps(data))
    return socket.recv()

@app.route("/set_calib")
def set_calib():
    dic = {"command": "set_calib"}
    socket.send(json.dumps(dic))
    return "ok"
    
@app.route("/calibrate")
def calibrate():
    dic = {"command": "calibrate"}
    socket.send(json.dumps(dic))
    return "ok"

@app.route('/search',methods=['POST'])
def search():
    data = json.loads(request.data)
    data["command"] = "search"
    socket.send(json.dumps(data))
    response = json.loads(socket.recv())
    return json.dumps(response)
    

# Move command. Can be used to request movement in [azimuth,altitude] directions.
@app.route('/move',methods=['POST'])
def move():
    data = json.loads(request.data)
    data["command"] = "move"
    socket.send(json.dumps(data))
    return 'ok'

# Set the target coordinate. Haywire unless calibrated already.
@app.route('/set_targ',methods=['POST'])
def set_targ():
    data = json.loads(request.data)
    data["command"] = "set_targ"
    socket.send(json.dumps(data))
    return 'ok'

# Returns the current azimuth altitude coordinates of the telescope. Currently useless.
@app.route('/get_pos')
def get_pos():
    data = {"command": "get_pos"}
    socket.send(json.dumps(data))
    #time.sleep(1)
    data = json.loads(socket.recv())
    return json.dumps(data)



# Set the users location as [longitude,latitude] and save it into a file
@app.route('/set_location',methods=['POST'])
def set_location():
    data = json.loads(request.data)
    data["command"] = 'set_location'
    socket.send(json.dumps(data))
    return "ok"

@app.route('/documentation')
def documentation():
    return auto.html(author='StelaSquad',
                     date="Last Update: 2018-02-18")



