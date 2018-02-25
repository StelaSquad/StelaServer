import os
import index
import unittest
import tempfile
from flask import jsonify

class IndexTestCase(unittest.TestCase):

    def setUp(self):
        """Sets up a small test client"""
        index.app.testing = True
        self.app = index.app.test_client()
        index.app.app_context()
        

    def tearDown(self):
        """ Nothing to tear down so far"""



    def test_get_coordinates(self):
        """Checks whether get_coordinates returns
            what is expected"""
        rv = self.app.get('/coordinates')
        expected_json = b'azimuth'
        assert expected_json in rv.data

if __name__ == '__main__':
    unittest.main()