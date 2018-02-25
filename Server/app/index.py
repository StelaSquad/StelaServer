from flask import Flask, jsonify, request
from flask import render_template
from flask_autodoc import Autodoc

app = Flask(__name__)
auto = Autodoc(app)
# TODO: use blueprints to seperate api across files


# Endpoint for obtaining info about the API
# Check flask_autodoc for usage
# TODO: Add a costumized template
@app.route('/documentation')
def documentation():
	return auto.html(author='StelaSquad',
					date="Last Update: 2018-02-18")


# Current coordinates of the telescope
# TODO: Create Coordinate Object
# TODO: probs want to store previous and current coordinates.
coordinates = {'azimuth' : [0, 0], 'equitorial': [0, 0]}

# Get current coordinates
@app.route('/coordinates')
@auto.doc()
def get_coordinates():
	"""Takes in coordinates and will return them for now"""
	return jsonify(coordinates)

# TODO: figure out how we will implement coordinates
# POST: coordinates needed for calibration
@app.route('/calibration', methods=['POST'])
@auto.doc()
def add_calibration_coordinate():
	"""Starts calibration procedure"""
    # request = request.get_json() - will get us the post info
	return '', 204
