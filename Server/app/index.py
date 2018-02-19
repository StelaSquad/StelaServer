from flask import Flask, jsonify, request

app = Flask(__name__)
# TODO: use blueprints to seperate api across files

# Endpoint for obtaining info about the API
# TODO: For better readability make this into a website
info = {'Last Update': '2018-02-18', 'authors': 'Stela Squad'}

@app.route('/info')
def get_info():
    return jsonify(info)

# Current coordinates of the telescope
# TODO: Create Coordinate Object
# TODO: probs want to store previous and current coordinates.
coordinates = {'azimuth' : [0, 0], 'equitorial': [0, 0]}

# Get current coordinates
@app.route('/coordinates')
def get_coordinates():
  return jsonify(coordinates)

# TODO: figure out how we will implement coordinates
# POST: coordinates needed for calibration
@app.route('/calibration', methods=['POST'])
def add_calibration_coordinate():
    # request = request.get_json() - will get us the post info
  return '', 204
