import os
import index
import unittest
import tempfile
from flask import jsonify, json

class IndexTestCase(unittest.TestCase):

    def setUp(self):
        """Sets up a small test client"""
        index.app.testing = True
        self.app = index.app.test_client()
        index.app.app_context()
        
        

    def tearDown(self):
        """ Nothing to tear down so far"""


    def test_calibration(self):
        data = {"azimuth": 1, "altitude": 2 }
        a=self.app.post('/calibration', data=data, 
            follow_redirects=True)
        
        assert b'altitude\": \"2\"' in a.data




    def test_get_coordinates(self):
        """Checks whether get_coordinates returns
            what is expected"""
        rv = self.app.get('/coordinates')
        expected_json = b'azimuth'
        assert expected_json in rv.data

    # def test_manual_control(self):
    #     test_client = index.app.test_client()

    #     response=self.app.post('/manual', 
    #                    data=json.dumps(dict(foo='Hey')),
    #                    content_type='application/json')

    #     assert b"Hey" in response.data

if __name__ == '__main__':
    unittest.main()