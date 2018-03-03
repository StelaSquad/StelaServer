import sys
sys.path.append('../')

from flask import Flask, jsonify, request
from flask import render_template
from flask_autodoc import Autodoc
from coms import coms
from scope import telescope, coordinates

tele = telescope.Telescope()
app = Flask(__name__)
auto = Autodoc(app)
# TODO: use blueprints to seperate api across files


# Endpoint for obtaining info about the API
# Check flask_autodoc for usage
# TODO: Add a costumized template
@app.route('/documentation')
def documentation():
    return auto.html(author='StelaSquad', date="Last Update: 2018-03-02")


# Current coordinates of the telescope
# Get current coordinates
@app.route('/coordinates')
@auto.doc()
def get_coordinates():
    """Takes in coordinates and will return them for now"""
    response = jsonify({
        'azimuth': tele.coordinates.azimuth,
        'altitude': tele.coordinates.altitude})

    return response

# TODO: figure out how we will implement coordinates
# POST: coordinates needed for calibration
@app.route('/calibration', methods=['POST'])
@auto.doc()
def add_calibration_coordinate():
    """Starts calibration procedure"""
    coms.testFunc()
    # request = request.get_json() - will get us the post info
    return '', 204

@app.route('/testPost', methods=['POST'])
@auto.doc()
def testPost():
    """Test Post endpoint"""
    return request.get_json()

@app.route('/testGet')
@auto.doc()
def testGet():
    """Test Post endpoint"""
    response = jsonify({
        'message': "Success"})
    return response, 200
