#!/usr/bin/python3
from pyorbital.orbital import Orbital
from datetime          import datetime
from datetime          import timedelta
from math              import radians
from math              import degrees
from math              import sin, cos, asin, acos

import matplotlib.pyplot as plt
import spacetrack
import requests
import sgp4

class Trajectory(object):


    result_x_coords      = []       #! We need this params to plot a graph
    result_y_coords      = []
    result_z_coords      = []
    thetas               = []
    phies                = []
    times                = []
    lat, lon, height = 0, 0, 0      #! This is lat, lon and height of a place 
                                    #! on Earth from which sat should be seen
    EARTH_RADIUS     = 6378.137

    def __init__(self, coords = (55.93013, 37.51832, 0.2), link='https://celestrak.com/NORAD/elements/active.txt', sat = "NOAA 19                 "):
        """ Initialisation of this class. Tle transcriotion and coords transformation """
        self.lat, self.lon, self.height = coords
        self.object_distance = self.EARTH_RADIUS + self.height

        # Coord transcription
        self.x  = self.object_distance * cos(radians(self.lat))
        self.x *= cos(radians(self.lon))
        self.y  = self.object_distance * sin(radians(self.lon))
        self.y *= cos(radians(self.lat))
        self.z = self.object_distance * sin(self.lat)

        #get satellite coordinats from tle
        with open('tle.txt', 'wb') as f:
            cont = requests.get(link, allow_redirects=True)
            f.write(cont.content)
            f.close()

        with open('tle.txt', 'r') as f:
            l = f.read().split('\n')
            
            for i, value in enumerate(l):
                if (value == sat):
                    print(sat, l[i + 1], l[i + 2], sep='\n')
                    self.tle = [sat, l[i + 1], l[i + 2]]


    
    def get_time_period(self):
        """gets time period of calculstion the orbit"""

        start_time_str = input('Enter a date to start with in format YYYY-MM-DD-HH-MM:    ')

        year, month, day, hour, minute = map(int, start_time_str.split('-'))

        start_time = datetime( year, month, day, hour, minute )

        finish_time_str = input('Enter a date to finish with in format YYYY-MM-DD-HH-MM:   ')

        year, month, day, hour, minute = map(int, finish_time_str.split('-'))

        finish_time = datetime( year, month, day, hour, minute )

        self.start_time  = start_time - timedelta(hours=3)
        self.finish_time = finish_time - timedelta(hours=3) 

        self.current_time = self.start_time

    def calculate_the_orbit(self):
        """Calculates the orbit of spaceship"""

        thetas = []
        phies  = []
        times  = []
        while self.current_time != self.finish_time:
            orb = Orbital("NONO", line1=self.tle[1], line2=self.tle[2])
            lon, lat, alt = orb.get_lonlatalt(self.current_time) 
            lon, lat = radians(lon), radians(lat)
            r = alt + self.EARTH_RADIUS

            x = r * cos(lat) * cos(lon)
            y = r * sin(lon) * cos(lat)
            z = r * sin(lat)

            self.result_x_coords.append(x)
            self.result_y_coords.append(y)
            self.result_z_coords.append(z)
            
            lk = self.x**2 + self.y**2 + self.z**2
            sk = self.x * x + self.y * y + self.z * z

            dist_to_fl  = ( sk - lk ) / ( lk ** 0.5 )
            dist_to_sat = ( (x - self.x) ** 2 + (y - self.y) ** 2 + (z - self.z) ** 2 )
            dist_to_sat = dist_to_sat ** 0.5

            theta = asin(dist_to_fl / dist_to_sat)
            phi   = 0

            if theta >= 0:
                norm_x = norm_y = 0
                norm_z = lk / self.z

                norm_vec = ( norm_x - self.x, norm_y - self.y, norm_z - self.z )

                norm_vec_len = norm_vec[0] ** 2 + norm_vec[1] ** 2 + norm_vec[2] ** 2
                norm_vec_len = norm_vec_len ** 0.5

                k = ( lk + sk ) / lk
                p_x = x + k * self.x
                p_y = y + k * self.y
                p_z = z + k * self.z

                p_vec = ( p_x - self.x, p_y - self.y, p_z - self.z )
                p_vec_len = p_vec[0] ** 2 + p_vec[1] ** 2 + p_vec[2] ** 2
                p_vec_len = p_vec_len ** 0.5

                phi = norm_vec[0] * p_vec[0] + norm_vec[1] * p_vec[1] + norm_vec[2] * p_vec[2]
                phi = acos(phi / (norm_vec_len * p_vec_len))

                phi = degrees( phi )
                
                est_vec = ( norm_vec[1] * self.z - norm_vec[2] * self.y,
                           norm_vec[2] * self.x - norm_vec[0] * self.z,
                           norm_vec[0] * self.y - norm_vec[1] * self.x ) 

                est_vec_len = sum(list(map(lambda x: x**2, est_vec)))
                
                phi1 = p_vec[0] * est_vec[0] + est_vec[1] * p_vec[1] + est_vec[2] * p_vec[2]

                phi1 = acos( phi1 / est_vec_len / p_vec_len )
                phi1 = degrees( phi1 )

                if phi1 > 90:
                    phi = -phi + 360

            thetas.append(theta)
            phies.append(phi)
            times.append(self.current_time + timedelta(hours=1))
            self.needed_elevation_azimuth(thetas, phies, times)
            self.current_time += timedelta(minutes=1)
            

    def needed_elevation_azimuth(self, thetas, phies, times):
        """ Some calculation of what elevations and azimuthes we need """
        _thetas = []
        _phies  = []
        _times  = []

        for i, value in enumerate(thetas):
            if value >= 0:
                _thetas.append(value)
                _phies.append( radians(phies[i]) )
                _times.append( times[i] )
            else:
                if _thetas:
                    self.thetas.append(_thetas)
                    self.phies.append(_phies)
                    self.times.append(_times)
                _thetas = []
                _phies  = []
                _times  = []

    def plot_some_sh___images(self):
        for i, value in enumerate(self.times):
            print(value, end=' ')
            print('Azimuth: ', str(degrees(self.phies[i][0]))[:-12] + 'gr', end=' ')
            print('Max elevation', str(max(self.thetas[i]))[:-12] + 'gr.')
            print()

        fig1 = plt.figure()
        ax   = fig1.add_subplot(111, projection='3d')
        ax.plot(self.result_x_coords, self.result_y_coords, self.result_z_coords)
        ax.scatter(self.x, self.y, self.z, color='red')
        fig1.set_size_inches(7, 7)

        fig2 = plt.figure()
        ax   = fig2.add_subplot(111, projection='polar')
        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)
        ax.set_rlim(bottom=90, top=0)
        for phi, theta in zip(self.phies, self.thetas):
            ax.plot(phi, theta)
        fig2.set_size_inches(7, 7)
        plt.show()

if __name__ == "__main__":
    tr = Trajectory()
    tr.get_time_period()
    tr.calculate_the_orbit()
    tr.plot_some_sh___images()
