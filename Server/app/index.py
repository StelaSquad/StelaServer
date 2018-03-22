from flask import Flask, jsonify, request
from flask import render_template
from flask_autodoc import Autodoc

# TODO: move these imports along with the 
# serail coms functionality to another file
import serial
import time
import STELA as st
s = st.STELA()
s.setup_cats()
s.setup_serial()

app = Flask(__name__)
auto = Autodoc(app)
app.config['DEBUG'] = False

stage = 0
# Setup code. Currently only gets reference stars and returns them
@app.route('/setup')
def setup():

	if stage == 0:
		try:
			s.load()
		except:
			return 'No location stored. Please send location to "/set_location"'
	
		coors = s.get_ref_stars('list')
		data = {'altaz_calib': coors}
		return jsonify(data)

# Move command. Can be used to request movement in [azimuth,altitude] directions.
@app.route('/move',methods=['POST'])
def move():
	data = request.get_json()
	[mvaz,mvalt] = data['move']
	s.get_pos()
	print mvaz, s.ard_az
	s.set_targ(mvaz+s.ard_az, mvalt + s.ard_alt)
	return 'ok'

# Set the target coordinate. Haywire unless calibrated already.
@app.route('/set_targ',methods=['POST'])
def set_targ():
	data = request.get_json()
	[az,alt] = data['targ']
	s.get_pos()
	s.set_targ(az,alt)
	return 'ok'

# Returns the current azimuth altitude coordinates of the telescope. Currently useless.
@app.route('/get_pos')
def get_pos():
	az,alt,targaz,targalt = s.get_pos(return_targ=True)
	data = {'pos':[az,alt],'targ':[targaz,targalt]}
	return jsonify(data)

# Set the users location as [longitude,latitude] and save it into a file
@app.route('/set_location',methods=['POST'])
def set_location():
	data = request.get_json()
	s.location = data['location']
	s.save()

@app.route('/documentation')
def documentation():
	return auto.html(author='StelaSquad',
					date="Last Update: 2018-02-18")



