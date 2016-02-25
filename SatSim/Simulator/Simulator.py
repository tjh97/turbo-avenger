"""
Created on Apr 21, 2015

@author: tjh97
"""

import numpy as np
import xml.etree.cElementTree as ET
from Satellite.Satellite import Satellite
from model.SatModel import SatModel


class Simulator(object):
    """

    """

    def __init__(self, info_file):
        """

        :param info_file:
        :return:
        """
        
        # Input variable for Satellite object
        sat_info = dict()
        
        # XML parsing
        tree = ET.parse(info_file)
        root = tree.getroot()
        
        # List of scalar parameters
        single_params = ['name', 'shape', 'dry_mass', 'fuel_mass', 
                         'drag_coeff', 'solar_coeff', 'rad_coeff']
        
        for param in single_params:
            try:
                sat_info[param] = float(root.find(param).text)
            except ValueError:
                sat_info[param] = root.find(param).text
            except:
                raise
        
        # List of vector parameters
        vect_params = ['length', 'cog', 'cop', 'moi', 'drag_area', 
                       'solar_area', 'rad_area']
        
        for param in vect_params:
            elist = root.findall(param)
            elist.sort(key=lambda x: x.attrib['dir'])
            vec_list = map(lambda x: float(x.text), elist)
            
            if len(elist) == 3:
                sat_info[param] = np.array(map(lambda x: [x], vec_list))
            elif len(elist) == 9:
                sat_info[param] = np.array([[vec_list[x], 
                                             vec_list[x+1], 
                                             vec_list[x+2] ] for x in range(2)])
            else:
                pass
        
        self.model = SatModel()
        self.satellite = Satellite(sat_info)
        
        
if __name__ == '__main__':
    Simulator('C:\\Users\\tjh97\\Dropbox\\Workspace\\SatelliteSimulator\\SatelliteInfo.xml')
    
    