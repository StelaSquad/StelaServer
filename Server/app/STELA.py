import numpy as np
import astroquery.simbad as aq
import astropy.coordinates as cp
import astropy.time as time
import astropy.units as u
import astropy.table as table
import os
import serial
import json
import warnings
import sys
<<<<<<< HEAD
import urllib2
from astropy.utils import iers
iers.conf.auto_download = False
Table = table.Table
=======
import datetime

def blockPrint():
    sys.stdout = open(os.devnull, 'w')

def enablePrint():
    sys.stdout = sys.__stdout__
>>>>>>> 45b2e9b5dfdb2af0c496fa36a0c1f761c525ca29

class STELA():
    """
    Class that governs the alignment of telescope, star positions, catalogs, and star identification
    
    Quantities and objects:
        STELA.naked: catalog of nearest stars, brightes in sky
        STELA.ard_pos: last stored position of arduino
        STELA.ard_targ: last stored arduino target
    
    Functions:
        STELA.setup_cats:
            Setup catalog files.
        STELA.setup_serial:
            Initiate serial connection with Arduino.
        STELA.get_ref_stars:
            Given an estimated latitude and longitude, returns the three brightest stars in the sky to align.
        STELA.gen_mock_obs:
            For testing triangulation
        STELA.triangulate:
            After calling get_ref_stars, and measuring the differences in alt-az coordinates of the three
            points, STELA.triangulate can locate the new latitude and longitude that accounts for the error
            in telescope positioning.
        STELA.set_targ:
            Sends target alt azimuth coordinate to arduino.
        STELA.get_pos:
            Reads current arduino poisition.
        
    """
    
    def __init__(self):
        self.savefile = 'steladata/savefile.sv'
        self.simbad = aq.Simbad()
        self.simbad.TIMEOUT = 1000000
        self.simbad.remove_votable_fields('coordinates')
        self.simbad.add_votable_fields('id(NAME)', 'ids','id(NGC)','id(M)','id(HD)','ra','dec','otype(V)',
                                       'sp','plx','z_value', 'flux(V)')
        self.reset_cats = False
        self.triangulation_class = Triangulate()
        self.calibrated = False
        cp.builtin_frames.utils.iers.conf.remote_timeout = 0.5

    
    def setup_cats(self):
        """ 
        Sets up the necessary catalogs, prints them to a file. (No parameters)
        
        """
        
        #os.chdir("./app")
        self.online = self.connect()
        
        if self.online == False:
            warnings.warn('No internet connection. Running offline mode.',Warning)
        
        # If user requests fresh import, do it
        if self.reset_cats == True and self.online == True:
            print "Deleting old catalog..."
            os.system('rm -r steladata')
        elif self.reset_cats == True: 
            print "Requested refresh of catalogs but no internet connection. Find one!"
            
        if os.path.exists('steladata') == False:
                os.mkdir('steladata')
                
        # Download necessary catalogs
        catlabels = ["GJ","New Galactic Catalog","Messier","Henry Draper"]
        cats = ["gj","ngc","m","hr"]
        catobjs = []
        
        for i in range(len(cats)):
            
            if os.path.exists('./steladata/' + cats[i].lower() + '.dat') == False:
                if self.online == False:
                    raise RuntimeError("No saved catalogs, run once with internet")
                    
                print "Downloading " + catlabels[i] +" data..."
                cat = self.simbad.query_catalog(cats[i].upper())
                select = np.array(cat["RA"] != '') *  np.array(cat["DEC"] != '')
                catobjs += [cat[select]]
                catobjs[-1].write('./steladata/'+ cats[i].lower() + '.dat',format='ascii')
                print "Done!"
                
            else:
                print catlabels[i] + " catalog file found."
                catobjs += [Table.read('./steladata/'+cats[i].lower()+'.dat',format='ascii')]
                
        [self.gj, self.ngc, self.m, self.hd] = catobjs
        catobjs = None
        
        
        # Create catalog of naked eye stars (used in calibration)
        if os.path.exists('./steladata/naked.dat') == False:
            
            print "Setting up naked eye catalogs"
            # remove objects with no recorded magnitude
            select = np.ones(len(self.gj),dtype='bool')
            select[np.where(np.isnan(np.array(self.gj['FLUX_V'])))[0]] = False
            query_naked_rm_nans = self.gj[select][:200]
            
            self.naked = Table(np.unique(query_naked_rm_nans))
            self.naked.sort("FLUX_V")
            
            # Write it to the catalogs folder
            self.naked.write('./steladata/naked.dat',format='ascii')
        else:
            self.naked = Table.read('./steladata/naked.dat',format='ascii')        
        
            
    def setup_serial(self):
        """ Searches for a serial connection on the coms. ttyS* ports corresponding to usb COM numbers in /dev/ 
        folder must be fully read/write/ex permissable. If failed, returns a Runtime Error. Some ports just
        don't work so switching USB ports might solve any problems."""
        
        ports_closed = []
        portsopen=0
        
        for i in range(20):
            # New path to try
            path = "/dev/ttyS" + str(i)
            
            try:
                # Send message, if recieve correct response, keep serial information
                ser = serial.Serial(path)
                ser.write('setup')
                
                if ser.readline()[:5] == 'STELA':
                    print 'Found STELA arduino running on on COM' + str(i)
                    self.COM = 5
                    self.ser = ser
                    portsopen+=1
                    break
                    
            except serial.SerialException as err:
                # Otherwise, check if error has to do with permissions.
                if err.args[0] == 13:
                    ports_closed += [i]
                
            except IOError:
                pass
          
        if os.path.exists("/dev/serial/by-id"):
            for i in os.listdir("/dev/serial/by-id"):
                print "trying other paths..."

                try:
                    ser = serial.Serial("/dev/serial/by-id/" + i)
                    ser.write("setup")
                    if ser.readline()[:5] == "STELA":
                        print "Found connection with: " + i
                        self.COM = "Unknown"
                        self.ser = ser
                        portsopen = 1
                        break
                except:
                    pass
                
        # If no serial port, raise error and if permission issues were found.
        if portsopen==0:
            if len(ports_closed) > 0:
                msg = "Connection Failed. Unable to access ports: " + str(ports_closed) + ". Try changing permissions."
            else:
                msg = "Connection Failed. Try different usb port."
                
            raise RuntimeError(msg)
                
    def search(self,string):
        
<<<<<<< HEAD
        print "Searching for: " + string
        self.online = self.connect(verbose=False)
        
        result = None
        matches = []
        
        try:
            # search all catalogs
            for cat in [self.m,self.gj,self.ngc,self.hr]:
                names = cat["IDS"]
                arr = names.view(type=np.recarray)

                low = np.array([s.lower() for s in arr])
                true = np.array([string.lower() in l for l in low])

                match = np.where(true==True)[0]

                if len(match) > 0:
                    matches += [cat[match]]

            matches = table.vstack(matches)

            #This algorithm is jenky but works somewhat
            scores = []
            for ma in matches:
                name_l = ma["IDS"].tolist().split("|")
                a = [string.lower() in name.lower() for name in name_l]
                possible = np.array(name_l)[a]
                lens = np.array([abs(len(string) - len(s)) for s in possible])
                scores += [ lens[np.where(lens==min(lens))] ]

            scores=np.array(scores)
            result = matches[np.where(scores == min(scores))[0]]
            
        except:
            pass
        
        try:
            coors = cp.get_body(string,time.Time.now())
            try:
                altaz = coors.transform_to(self.home)
                self.set_targ([coors.az.deg,coors.az.deg])
                return {"Name": string}
            except:
                return {"Name": string}
        except:
            pass
            
        # Check connectivity
        if self.online == True and type(result) == type(None):
            result = self.simbad.query_object(string)

        if type(result) == type(None):
            return {"Error":"No object found"}
        [main, name, ids, idngc, idm, 
         idhd, ra, dec, otype, sptype, 
         plx, redshift, mag]          = result.as_array()[0]
            
        
        usenames=[]
        for i in [name,idm,idngc,idhd]:
            if i != '':
                usenames += [i]
                
        data = {"Name": usenames}
        if mag != None:
            data["Mag"] = str(mag)
        if redshift != None:
            data["Redshift"] = str(redshift)
        if otype != None:
            data["Otype"] = otype
        else:
            data["Otype"] = "Unknown"
        if sptype != '':
            data["Sptype"] = str(sptype)
        if plx != None:
            data["Plx"] = str(plx)
            distpc = 1000./plx
            distpcrnd = round(distpc,3)
            distly = round(distpc*3.2616,3)
            data["Distance"] = str(distpcrnd) + " pc, " + str(distly) + " ly"
                
        if plx != None and mag != None:
            Msun = 4.74
            absmag = mag - 5 *(np.log10(distpc) - 1)
            L = 10 ** (1./2.5 * (Msun - absmag))
            data["Luminosity"] = str(L)
            
        if self.calibrated == True:
            celcoor = cp.SkyCoord(ra,dec,unit=[u.hourangle,u.deg])
            altazcoor = celcoor.transform_to(self.home)
            self.set_targ([altazcoor.az.deg,altazcoor.alt.deg])
                
        return data
            
        
    def set_time(self,datetime):
        l = datetime.split("-")
        date_str = l[1]+l[2]+l[3]+l[4]+l[0]+"."+l[5]
        if os.environ["USER"] == 'pi':
            os.system("sudo date " + date_str)
        else:
            print "Not setting time."

    def set_targ(self, array):
        """ Sends the target coordinates in the arduino. """
        
        msg = 'set_targ:' + str(array)
=======
    def set_time(timestr):
        import ctypes
        import ctypes.util
        import time


        timestr = [int(x) for x in timestr.strip().split('-')]



        # /usr/include/linux/time.h:
        #
        # define CLOCK_REALTIME                     0
        CLOCK_REALTIME = 0

        # /usr/include/time.h
        #
        # struct timespec
        #  {
        #    __time_t tv_sec;            /* Seconds.  */
        #    long int tv_nsec;           /* Nanoseconds.  */
        #  };
        class timespec(ctypes.Structure):
            _fields_ = [("tv_sec", ctypes.c_long),
                        ("tv_nsec", ctypes.c_long)]

        librt = ctypes.CDLL(ctypes.util.find_library("rt"))

        ts = timespec()
        ts.tv_sec = int( time.mktime( datetime.datetime( *time_tuple[:6]).timetuple() ) )
        ts.tv_nsec = time_tuple[6] * 1000000 # Millisecond to nanosecond

        # http://linux.die.net/man/3/clock_settime
        librt.clock_settime(CLOCK_REALTIME, ctypes.byref(ts))


     
    def set_targ(self, az, alt):
        """ Sends the target coordinates in the arduino. """
        msg = 'set_targ:' + str([az,alt])
>>>>>>> 45b2e9b5dfdb2af0c496fa36a0c1f761c525ca29
        self.ser.write(msg)
        
    def get_pos(self,return_targ=False):
        """ Gets the current arduino position and target. """
        now = time.Time.now()
        
        # Send request to arduino, read response
        self.ser.write('info')
        msg = self.ser.readline()
        
        # Find alt-az information
        pos_str = msg[ msg.find('[')+1 : msg.find(']') ]
        targ_str = msg[ msg.find('[',11)+1 : msg.find(']', msg.find(']') + 1) ] 
        
        azalt = np.fromstring(pos_str,sep=', ')
        targ = np.fromstring(targ_str, sep=', ')
        
        self.ard_pos = azalt
        self.ard_targ = targ
        self.update_time = now

        if return_targ==True:
            return self.ard_pos, self.ard_targ
        
        return self.ard_pos
    
    def set_pos(self,newpos):
        msg = 'set_pos:' + str(newpos)
        self.ser.write(msg)
        
    def get_ref_stars(self,representation='SkyCoord'):
        """
        Given an estimation of longitude and latitude, identifies 3 target stars to use as 
        triangulation coordinates
        
        Input:
            ------
            lon_est: float
                The current longitude at which the telescope is set up
            
            lat_est: float
                The current latitude at which the telescope is set up
        
        Optional:
            ------
            representation: string
                The preferred representation of the coordinates ('SkyCoord' or 'String')
        Output:
            ------
            altaz_calib: list [args, shape=(3,2), dtype=float]
                The estimated positions in the altitude azimuth coordinate frame. Can be used to point
                telescope to approximate position of stars.
        """
        
        if self.connect(verbose=False) == False:
            cp.builtin_frames.utils.iers.conf.auto_download=False 
        [lon_est, lat_est] = self.location
        
        # Set observation time for calibration
        self.time = time.Time.now()
        
        # Set up approximate earth frame
        loc_est = cp.EarthLocation(lon_est,lat_est)
        earth_n = cp.AltAz(0*u.deg,90*u.deg,location=loc_est,obstime=self.time)
        
        # Coordinates of all visible stars
        ra= self.naked["RA"].data
        dec = self.naked["DEC"].data
        
        cat = cp.SkyCoord(ra = ra,dec = dec,unit=['hourangle','deg'])
        
        print 1
        earth_n_cel = earth_n.transform_to(cat)
    
        print 2 

        # Choose stars that are at least 30 degrees above horizon
        select = earth_n_cel.separation(cat) < cp.Angle(60*u.deg)
        cat_close = cat[select]
        
        print 3

        # Select the first three
        coors=[cat_close[0]]
        for i in cat_close[1:]:
            if sum([i.separation(c) < 10*u.deg for c in coors]) == 0:
                
                coors+=[i]
                if len(coors) >= 3:
                    break
        print 4
                    
        cel_calib = cp.SkyCoord(coors)
        
        # coordinates in earth frame
        altaz_calib = cel_calib.transform_to(earth_n)
        
        # set class objects
        self.altaz_calib = altaz_calib
        self.cel_calib = cel_calib
        
        #print earth_n.transform_to(cat)
        if representation == 'SkyCoord':
            return altaz_calib
        elif representation == 'String' or representation == 'string':
            return altaz_calib.to_string()
        elif representation == 'list':
            d = np.array([altaz_calib.az.deg,altaz_calib.alt.deg]).T
            return d.tolist()
    
    def gen_mock_obs(self):
        """
        Generates fake observation data based on possible errors of +/- 5deg.
        
        Output:
            ------
            obs: list [args,len=2,dtype=float]
                Difference in altitude and azimuth for three points, used in triangulation.
        
        """
        
        # creates class object based on random errors in telescope placement
        self.mock_home = [self.lon_est + np.random.uniform(-5,5)*u.deg,
                          self.lat_est + np.random.uniform(-5,5)*u.deg]
        
        # creates objects based on errors
        mock_loc = cp.EarthLocation(self.mock_home[0],self.mock_home[1])
        surf = cp.AltAz(location=mock_loc,obstime=self.time)
        pts = self.altaz_calib.transform_to(surf)
        
        # calculate angular differences in altitude and azimuth for points
        self.v2_v1 = [pts[1].az.rad - pts[0].az.rad, pts[1].alt.rad - pts[0].alt.rad]
        self.v3_v2 = [pts[2].az.rad - pts[1].az.rad, pts[2].alt.rad - pts[1].alt.rad]
        
        return [self.v2_v1, self.v3_v2]
    
    def triangulate(self, v2_v1,v3_v2,iterations=5):
        """
        Used to triangulate the true latitude and longitude corresponding to the norm of the telescope position.
        
        Input
            ------
            v2_v1: list [args,len=2,dtype=float]
                The difference in [azimuth,altitude] between object 2 and object 1
               
            v3_v2: list [args,len=2,dtype=float]
                The difference in [azimuth,altitude] between object 3 and object 2
        """
        
        ra = cp.Angle(self.cel_calib.ra, unit=u.hourangle)
        dec = cp.Angle(self.cel_calib.dec, unit=u.deg)

        v = np.array([ra.rad,dec.rad]).T
        
        try:
            out = self.triangulation_class.triangulate(v[0],v[1],v[2],v2_v1,v3_v2,iterations=iterations)
        except Exception as e:
            print e
            

        n = cp.SkyCoord(out[0][0],out[0][1],unit=u.rad,frame='icrs')

        npr = cp.EarthLocation(lon=0*u.deg,lat=90*u.deg,height=0*u.m)

        n.location=npr
        n.obstime=self.time

        [self.lon,self.lat] = [180-n.altaz.az.deg, n.altaz.alt.deg]
        self.home_coors = [self.lon,self.lat]
        #self.save()
        
        h = cp.EarthLocation(self.home_coors[0]*u.deg,self.home_coors[1]*u.deg)
        self.home = cp.AltAz(location=h,obstime=self.time)
        v3 = cp.SkyCoord(v[2][0],v[2][1],unit=u.rad)
        self.tel_pos = v3
        self.set_pos(v3)
        
        self.calibrated = True

        
    def save(self):
        [lon,lat] = self.location 
        savedata = {'location': [lon.value,lat.value], 
                    'location_units': [lon.unit.to_string(), lat.unit.to_string()]}
        
        with open(self.savefile,'w') as file:
            file.write(json.dumps(savedata))
            
    def load(self):
        with open(self.savefile,'r') as file:
            data = json.loads(file.read())
            
        loc = data['location']
        loc_u = data['location_units']
        
        self.location = [loc[0]*u.Unit(loc_u[0]),loc[1]*u.Unit(loc_u[1])]
        
        
    def connect(self,verbose=True):
         
        if verbose == True:
            print "Testing internet connection..."

        try:
            urllib2.urlopen('http://216.58.192.142', timeout=1)
            if verbose == True:
                print "Connection succeeded!"
            return True
        except urllib2.URLError as err:
            if verbose == True:
                print "Connection failed."
            return False

    
    
    
    
    
    
    
  
    
    
    
    
    
    
    
    
    
    
    
    
    

class Triangulate():
    """
    This class is used to triangulate the location of the telescope based on three celestial objects. 
    
    It operates based on a latitude and longitude equatorial coordinate system, with celestial coordinates
    given in that frame. The telescope can be placed a floor that is not perfectly normal to the surface 
    of the earth. The use a (technically) false latitude and longitude, one for which the normal vector is the
    same as the normal vector of the base of the telescope, can account for this error. This is the method
    used, but note that the outputted earth coordinates will likely differ slightly than the true location of
    the telescope.
    
    Objects:
        Triangulate.obserrs: 
                Expected errors in measured altitude and azimuth (default=1')
        Triangulate.comp_power: T
                he calculative potential of the computer being used (default=100)
    
    Functions:
        Triangulate.triangulate:
                Computes the normal vector of the telescope based on three objects and their observed angular
                differences.
        
    Helper Functions:
        Triangulate.xp, Triangulate.yp, Triangulate.zp: 
                The coordinate vectors of a frame on the surface of the earth, at location [a,b].
        Triangulate.vec_prime: 
                Transforms a vector in equatorial frame into an earth alt-az frame at location [a,b].
        Triangulate.find_valid: 
                Given two objects and their angular difference, calculates probability the distribution across
                latitude and longitude space.
        Triangulate.match:
                Given one coordinate grid and three 2d probability distributions, calculates the total 
                probability distribution combining the three points.
        
    Outputs:
        Triangulate.lon: 
                Calculated longitude
        Triangulate.lat: 
                Calculated latitude
        Triangulate.errs: 
                The errors of the latitude and longitude [errlon,errlat]
    
    
    """
    def __init__(self):
        self.comp_power = 100
        self.chain=[]
        self.dist=[]
        self.dth=0
        self.dph=0
        self.obserrs = [np.pi/180/30,np.pi/180/30]

    
    def xp(self,a,b):
        """
        Returns x-prime, the vector representing the primed x-axis in the celestial frame. 
        Points toward the sky, normal to surface of earth.
        
        Input:
            ------
            a: float 
                Alpha, the angle that the primed coordinate system is spun about its z-axis.
               
            b: float
                Beta, angle that the primed coordinate system is pitched about its y-axis.
                
        Output:
            ------
            x_prime: ndarray [arg, dtype=float, shape=(3)
                3-vector that represents the x-axis of the prime frame in celestial xyz coordinates.
        """
        
        x_prime = np.array([np.cos(a)*np.cos(b),np.sin(a)*np.cos(b),np.sin(b)])
        return x_prime

    def yp(self,a,b):
        """
        Returns y-prime, the vector representing the primed x-axis in the celestial frame. 
        Points toward the east, tangent to surface of earth.
        
        Input:
            ------
            a: float 
                Alpha, the angle that the primed coordinate system is spun about its z-axis.
               
            b: float
                Beta, angle that the primed coordinate system is pitched about its y-axis.
                
        Output:
            ------
            y_prime: ndarray [arg, dtype=float, shape=(3)
                3-vector that represents the y-axis of the prime frame in celestial xyz coordinates.
        """
    
        y_prime = np.array([-np.sin(a),np.cos(a),0])
        return y_prime
    
    def zp(self,a,b):
        """Returns z-prime, the vector representation of the primed z-axis in the celestial frame.
        Points north, tangent to surface of earth
        
        Input:
            ------
            a: float 
                Alpha, the angle that the primed coordinate system is spun about its z-axis.
               
            b: float
                Beta, angle that the primed coordinate system is pitched about its y-axis.
                
        Output:
            ------
            z_prime: ndarray [arg, dtype=float, shape=(3)
                3-vector that represents the z-axis of the prime frame in celestial xyz coordinates.
        """

        z_prime = np.array([-np.cos(a)*np.sin(b),-np.sin(a)*np.sin(b),np.cos(b)])
        return z_prime


    def vec_prime(self, a, b, v, form='xyz'):
        """
        Transforms a vector in equatorial celestial frame into an earth alt-az frame at location [a,b].
        
        Inputs:
            ------
            a: float 
                Alpha, the angle that the primed coordinate system is spun about its z-axis.
               
            b: float
                Beta, angle that the primed coordinate system is pitched about its y-axis.
                
            v: ndarray or list [arg, dtype = float, shape=(2) or shape=(3)]
                A celestial vector, in either a xyz or theta-phi representation.
                
        Optional:
            ------
            rep: string, arg = 'xyz' or 'th-ph'
                Specify whether the output data should be in xyz coordinates or theta-phi coordinates. 
                
        Output:
            v_prime: ndarray [arg, dtype=float, shape=(3) or shape=(2)]
                The vector represented in the primed coordinate system.
        """
        
        v = np.array(v)
        
        if len(v) == 2:
            v = self.xyz(v[0],v[1])
            
            
        x = np.dot(self.xp(a,b),v)
        z = np.dot(self.zp(a,b),v)
        y = np.dot(self.yp(a,b),v)
        
        v_xyz = np.round(np.array([x,y,z]),12)

        
        if form == 'xyz':
            v_out = v_xyz
            
        elif form == 'th-ph':
            th = np.arctan2(v_xyz[1],v_xyz[2])
            ph = np.arcsin(v_xyz[0])
            
            
            if th < 0:
                th = 2*np.pi+th
            
            v_out = np.array([th,ph])
            
        else:
            raise ValueError("Requested representation not understood. Use either 'xyz' or 'th-ph")
            
        return v_out


    def xyz(self,th,ph):
        """Transforms spherical coordinates into xyz coordinates in the celestial frame."""
        return np.array([np.cos(th)*np.cos(ph), np.sin(th)*np.cos(ph), np.sin(ph)])
    
    
    def gen_mock(self,errs=[0,0]):
        """
        Generates a mock triangulation data set, to test the program.
                
        Optional:
            ------
            errs: list or ndarray [args, dtype=float shape=(2)] 
                list the error in the the generated data of the observed angle differences
        
        Output:
            v_prime: ndarray [arg, dtype=float, shape=(3) or shape=(2)]
                The vector represented in the primed coordinate system.
        """
        
        # Generate random truth values
        a = np.random.rand()*2*np.pi
        b = (np.random.rand()-0.5)*np.pi
        
        # Generate three random triangulation vectors
        v1 = [np.random.uniform(0,2*np.pi),np.random.uniform(-np.pi/2,np.pi/2)]
        v2 = [np.random.uniform(0,2*np.pi),np.random.uniform(-np.pi/2,np.pi/2)]
        v3 = [np.random.uniform(0,2*np.pi),np.random.uniform(-np.pi/2,np.pi/2)]
        
        # Calculate the true difference in alt-azimuth, and add a little random error
        v2_v1 = self.vec_prime(a,b,v2,form='th-ph') - self.vec_prime(a,b,v1,form='th-ph') + np.random.uniform(-errs[0],errs[0])
        v3_v2 = self.vec_prime(a,b,v3,form='th-ph') - self.vec_prime(a,b,v2,form='th-ph') + np.random.uniform(-errs[1],errs[1])
                
        
        return [a,b],v1,v2,v3,v2_v1,v3_v2

    
    

    def find_valid(self, obj1coor, obj2coor, obs, lims=[[0,2*np.pi],[-np.pi/2,np.pi/2]]):
        """
        Calculates the probability distribution of lattitude and longitude points within the given limits.
        
        Inputs:
            ------
            obj1coor: list or ndarray [args, dtype=float, shape=(2)]
                Celestial angle coordinates for the first object
                
            obj2coor: list or ndarray [args, dtype=float, shape=(2)]
                Celestial angle coordinates for the second object
                
            obs: list or ndarray [args, dtype=float, shape=(2)]
                the observed difference in [azimuth, altitude] between object 1 and 2
                
        Optional:
            ------
            lims: list or ndarray [args, dtype=float, shape=(2,2)]
                The limits in which to look for valid longitude and lattitudes. Default is the entire space.
            
        Output:
            ------
            grid: list [args,dtype=float,shape=(2,100,100)]
                The coordinates of each longitude and lattitude tested
            
            dist_norm: ndarray [args, dtype=float, shape=(100,100)]
                The normalized probability distribution over all points in the grid.
                
        """
        
        #print obs*np.pi/180
        
        n=self.comp_power
    
        obsth = obs[0]
        obsph = obs[1]
        
        # Initalize variable space
        A = np.linspace(lims[0][0],lims[0][1],n)
        B = np.linspace(lims[1][0],lims[1][1],n)
        grid = np.meshgrid(A,B)
                    
        # calculate observation vectors
        obj1_xyz = self.xyz(obj1coor[0],obj1coor[1]) # xyz(obj1[0],obj1[1])
        obj2_xyz = self.xyz(obj2coor[0],obj2coor[1]) # (obj2[0],obj2[1])

        # For each potential A and B vector, calculate the theoretical change in theta and phi
        
        diff = lambda A1,B1: (np.array(self.vec_prime(A1,B1,obj2_xyz,form="th-ph")) -
                              np.array(self.vec_prime(A1,B1,obj1_xyz,form="th-ph")))
        
        [thp,php] = np.array(map(diff,grid[0].flatten(),grid[1].flatten())).T
        
        #back from column to 2d grid
        thp = np.array(thp).reshape(n,n)
        php = np.array(php).reshape(n,n)
        
        # Start real small with the binsize, extremely restrictive
        stdth = np.std(thp.flatten())
        stdph = np.std(php.flatten())
                
        mod = (self.obserrs[0] + self.obserrs[1])/2
        
        while mod*stdth < 3*self.obserrs[0] and mod*stdph < 3*self.obserrs[1]:
            
            mod*=1.05
        
        dA = A[1]-A[0]
        dB = B[1]-B[0]
        
        thdist = np.exp(- ((thp-obsth)/(mod*stdth))**2 )
        phdist = np.exp(- ((php-obsph)/(mod*stdph))**2 )
                
        dist = (thdist)*(phdist)
        dist_norm = dist/np.sum(dist*dA*dB)
        
        #self.chain+=[dist_norm]
        
        return grid,dist_norm
    
    def match(self,grid, c1,c2,c3):
        """
        Combines three probability distributions to isolate points for which the lattitude and longitude 
        values match observation
        
        Inputs
            ------
            c1: ndarray [args,shape=(n,n),dtype=float]
                Normalized distribution as outputted by find_valid.
        
        Outputs:
            stat_array: ndarray [args,shape=(2,2),dtype=float]
                A list consisting in the expected values of the latitude and longitude and the errors.
        """
        
        # We want to compare two functions whose data points are not necessarily the same. Thus, we have to 
        # bin both into a new, global data set.
        
        comb = c1*c2*c3
        
        comb /= np.sum(comb)
                
        self.dist+=[[grid,comb]]
        
        thn = np.sum(comb,axis=0)
        phn = np.sum(comb,axis=1)
        
        th_ax = grid[0][0]
        ph_ax = grid[1][:,0]
        
        th_avg = np.sum(thn*th_ax)
        ph_avg = np.sum(phn*ph_ax)
        
        
        th_std = np.sqrt(np.sum( (th_ax - th_avg)**2 * thn))# * Np/(Np-1) )/ np.sqrt(n)
        ph_std = np.sqrt(np.sum( (ph_ax - ph_avg)**2 * phn))# * Np/(Np-1) )/ np.sqrt(n)
        
        return np.array([[th_avg,ph_avg],[th_std,ph_std]])

    
    def triangulate(self, v1, v2, v3, obs_v1_v2, obs_v2_v3, iterations = 5, obserrs=[None,None]):
        """
        Given v1,v2,v3, equatorial celestial coordinates for three objects, and observed changes in altitude 
        and azimuth of those objects from the ground, returns coordinates of normal vector. i.e. latitude
        and longitude
        
        Input:
            ------
            v1: list [args,shape=(2),dtpe=float]
                The celestial coordinates of the first triangulation object
                
            v2: list [args,shape=(2),dtpe=float]
                The celestial coordinates of the second triangulation object
                
            v3: list [args,shape=(2),dtpe=float]
                The celestial coordinates of the third triangulation object
                
            v2_v1: list [args,shape=(2),dtpe=float]
                The difference in azimuth and altitude between object 2 and object 1
                
            v3_v2: list [args,shape=(2),dtpe=float]
                The difference in azimuth and altitude between object 3 and object 2
        """
        
        if sum(np.array(obserrs)==None) == 0:
            self.obserrs = obserrs
        
        # Calculate difference between v1 and v3
        obs_v1_v3 = [obs_v1_v2[0]+ obs_v2_v3[0], obs_v1_v2[1]+ obs_v2_v3[1]]
        
        lims = [[0,2*np.pi],[-np.pi/2,np.pi/2]]
        
        i=0
        while i < iterations:
            print "Running for lims: " + str(np.round(lims,5).tolist())
            
            # find the probability distributions for each observation
            grid, c1 = self.find_valid(v1, v2, obs_v1_v2, lims=lims)
            _, c2 = self.find_valid(v1, v3, obs_v1_v3, lims=lims)
            _, c3 = self.find_valid(v2, v3, obs_v2_v3, lims=lims)
            
            if np.sum(np.isnan(c1*c2*c3)) > 0 :
                print "Error in triangulation. Increasing obserrs in attempt to mend situation."
                self.obserrs = [self.obserrs[i]*10 for i in [0,1]]
                
            else:

                # Matches all three
                [av,acc] = self.match(grid,c1,c2,c3)


                # Finds the accuracy of the analysis, chooses new limits based on these
                r = 5
                dth = grid[0][0][1]-grid[0][0][0]
                dph = grid[1][1][0]-grid[1][0][0]

                acc += np.array([dth,dph])/(r)

                lims = np.array([av - r*acc, av + r*acc]).T
                
                i+=1
                
                        
        try: 
            self.lon = av[0]
            self.lat = av[1]
            self.errs = acc

            print "Done."
            return av,acc
        
        except:
            print "Unable to triangulate"
            return None
