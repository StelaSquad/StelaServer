import requests
import time
import json


class telecon:
    def __init__(self):
        # Specify flask ip adress here (this one works if flask is run on your own computer)
        self.url = 'http://192.168.0.100:5123'

    def set_location(self,lon_lat,units):
        data = {'location': lon_lat,'units':units}
        resp = requests.post(self.url+'/set_location',data=json.dumps(data))
        response = resp.content
        return response
        
    def set_targ(self,az_alt):
        # Set the target position in the arduino (notice dictionary protocol)
        data = {'targ': az_alt}
        resp = requests.post(self.url+'/set_targ',data=json.dumps(data))
        return resp.content
        
    def get_pos(self):
        # Get the current position of the arduino (again, protocols)
        resp = requests.get(self.url+'/get_pos')
        dic = json.loads(resp.content)
        az_alt = dic['pos']
        return az_alt
    
    def setup(self):
        
        t = time.localtime()
        string = ''
        for i in t[:6]:
            bit = str(i)
            if len(bit)==1:
                bit = '0'+bit
            string += bit + '-'
        string += '000'
        # Currently nonfuntional
        data = {'time': string}
        
        # Currently nonfuntionalipython
        a = requests.post(self.url+'/setup',data=json.dumps(data))
        return json.loads(a.content)["calib_coors"]
    
    def move(self,az_alt):
        # move a small incriment
        data = {'increment': az_alt}
        req = requests.post(self.url+'/move',data=json.dumps(data))
        response = req.content
        return response
        
    def set_calib(self,n=None):
        if n==None:
            requests.post(self.url+'/set_calib')
        else:
            data = {"calib_num": n }
            req = requests.post(self.url+'/set_calib',data=json.dumps(data))
            response = req.content
            
    def calibrate(self):
        requests.get(self.url+'/calibrate')
        
