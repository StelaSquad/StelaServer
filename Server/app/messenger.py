# Python middleman for controlling STELA telescope. Meant to be run in unison with the app.
# We need a middleman because the target on the arduino continuously needs to be updated when 
# star tracking, which means a while True loop needs to be implemented
# 
# Duncan Rocha, April 2018

import zmq
import time
import json
import STELA
import astropy.units as u

# Set up communication between app and this script
port = "5122"
context = zmq.Context()
socket = context.socket(zmq.PAIR)
socket.bind("tcp://*:%s" % port)

# Start STELA library
st = STELA.STELA()
st.setup_cats()
st.setup_serial()

calib_obs = [[],[],[]]
c=0
pi = 3.1415

# Loop it up
while True:
    
    # Check for messages from the app
    msg = json.loads(socket.recv())
    
    if len(msg) > 0:
        # Check command
        command = msg["command"]
        print command
        
        # List of possible commands and the appropriate responses.
        if command == "setup":
            
            st.load()
            st.set_time(msg["time"])
            coors = st.get_ref_stars(representation='list')
            dic = {"calib_coors": coors}
            socket.send(json.dumps(dic))
            
            print st.cel_calib.to_string()
        
        if command == "set_calib":
            try:
                n = msg["calib_num"]
            except:
                n=c 
            calib_obs[n] = st.get_pos()
            c+=1
            
        if command == "calibrate":
            
            print "Triangulating..."
            
            v2_v1 = [(calib_obs[1][i]-calib_obs[0][i])*pi/180 for i in [0,1]]
            v3_v2 = [(calib_obs[2][i]-calib_obs[1][i])*pi/180 for i in [0,1]]
            st.triangulate(v2_v1,v3_v2)
            
            pos = st.cel_calib[2].transform_to(st.home)
            st.set_pos([pos.az.deg,pos.alt.deg])
            
            print "Telescope position triangulated to: " + str(st.home_coors)      
            
        if command == "set_targ":
            targ = msg["targ"]
            st.set_targ(targ[0],targ[1])
            print "Target Set: " + str(targ)
        
        if command == "get_pos":
            print "getting position"
            pos = st.get_pos()
            dic = {"pos": pos.tolist()}
            socket.send(json.dumps(dic))
            print "Position Returned: " + str(pos)
            
        if command == "move":
            mv = msg["increment"]
            ard_targ = st.get_pos(return_targ = True)[1]
            newtarg = [ard_targ[i] + mv[i] for i in [0,1]]
            st.set_targ(newtarg[0],newtarg[1])
            print "Position target incremented by: " + str(mv)
            
        if command == "set_location":
            st.location = [msg["location"][i]*u.Unit(msg["units"][i]) for i in [0,1]]
            st.save()
            
            
            
            
            
            
            
            
            
            
            
