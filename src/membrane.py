#! /usr/bin/env python3

"""
Module membrane

This module define membrane class with it attributs and associated methods. It also defines residue class.
"""

#Import modules
##########################################################################################################

import os

try:
	import numpy as np
except ImportError as exception:
	print("ERROR: Module Numpy is not installed. Advice: Advice: Use the transs_prot_pred.yaml file to create a powerfull environnment")

#Define class
##########################################################################################################

class Membrane:
    """
    This class define some membrane properties
    """
    def __init__(self, resolution, side, plane_coefficients, mass_center):
        """
        This fucntion assign some instance attribtus for a membrane
        input:
        resolution          float   space between membrane point
        side                float   side length of the membrane
        plane_coefficients  list    list of optimal membrane plane parameters
        center              list    list of mass center portein coordinates
        axis_plane          list    reference axis plane for membrane
        """
        self.resolution = resolution
        self.side = side
        self.plane_coefficients = plane_coefficients
        self.center = mass_center
        self.points_plane_1 = []
        self.points_plane_2 = []
        
    def get_points(self):
        """
        This function allows to add points on the membrane plane
        """
        #Extract mass_center coordinates
        x_center, y_center, z_center = self.center[0], self.center[1], self.center[2]
        #Define x_array and y_array
        x_array = np.arange(x_center - self.side, x_center + self.side, self.resolution)
        y_array = np.arange(y_center - self.side, y_center + self.side, self.resolution)
        #Extract plane coefficient
        a, b, c, d = self.plane_coefficients[0], self.plane_coefficients[1], \
        self.plane_coefficients[2], self.plane_coefficients[3]
        for x in x_array:
            for y in y_array:
                #Compute z value for plane 1
                z = ((-1*a*x) + (-1*b*y) + (-1*d)) / c
                self.points_plane_1.append([x, y, z])
                #Compute z value forplane 2
                d = -1*d
                z = ((-1*a*x) + (-1*b*y) + (-1*d)) / c
                self.points_plane_2.append([x, y, z])


#Define function
###########################################################################################################

def create_pdb_file(input_file, output_file, membrane):
    """
    This function allows to implement membrane in pdb file
    input:
    input_file      str     path to pdb file
    output _file    str     path of the new pdb file
    membrane        object  membrane object to add in the pdb file
    """
    #Create copy 
    os.system(f"cp {input_file} {output_file}")
    #Extract points list
    plane_1 = membrane.points_plane_1
    plane_2 = membrane.points_plane_2
    #Modify copy
    with open(output_file, "a") as pdb_output:
        #Write points
        for i in range(len(membrane.points_plane_1)):
            line = f"{'HETATM':6s}{i + 1:5d} {'DUM':^4s}{'1':1s}             {plane_1[i][0]:8.3f}{plane_1[i][1]:8.3f}{plane_1[i][2]:8.3f}\n"
            pdb_output.write(line)
            f"{'HETATM':6s}{i + 1:5d} {'DUM':^4s}{'2':1s}             {plane_1[i][0]:8.3f}{plane_1[i][1]:8.3f}{plane_1[i][2]:8.3f}\n"
            pdb_output.write(line)                  
