#! /usr/bin/env python3

"""
Module sphere

This module define sphere class with it attributs and associated methods.
It can be used in order to create a sphere with surface points taht embaed a whole protein.
Surface points will be fairly distributed thank to Fibonacci algorithm
"""

#Import module
##########################################################################################################

import math

#Define classes
###########################################################################################################

class Sphere:
    """
    This class define some element of a half sphere which encompass a protein
    """
    def __init__(self, x, y, z, radius):
        """
        This function assigns surface points and center coordinates attributs
        input:
        x       float   x center coordinate
        y       float   y center coordinate
        z       float   z center coordinate
        radius  float   sphere radius
        """
        #Define surface points list
        self.surface_points = []
        #Define center coordinate and radius
        self.x_center, self.y_center, self.z_center = x, y, z
        self.global_radius = radius
        #Define constant which will be useful to generate points
        self.PHI = (1 + math.sqrt(5)) / 2
        self.GOLDEN_ANGLE = (2 * math.pi) / (self.PHI**2)
        #Assign axis
        self.vertical_axis = "y"
        self.cos_axis = "x"
        self.sin_axis = "z"
        #Assign max coordinate value in vertical axis
        self.y_max = self.y_center + self.global_radius

    def get_vertical_incrementation(self, point_number):
        """
        This function compute incrmentation on vertical axis in fucntion of desired
        point number.
        input:
        point_number    int     desired point number on the half sphere 
        """
        self.point_number = point_number
        self.incrementation = self.global_radius / point_number
        
    def add_fibonacci_point(self, incrementation_value):
        """
        This function allows to add a point on the sphere surface using Fibonacci
        algorithm. It must be use in a loop in order to generate many point.
        input:
        incremntation_value  int  a integer providing by the loop (for i in range(point_number)) 
        """
        #Define coordinate on vertical axis
        y = self.y_max - (incrementation_value * self.incrementation)
        #Define local radius which corresponding to a cercle plane
        local_radius = math.sqrt(self.global_radius**2 - (y - self.y_center)**2)
        #Define angle in order to place it on the surface
        local_angle = self.GOLDEN_ANGLE * incrementation_value
        #Define cos_axis and sin_axis coordinate
        x = self.x_center + math.cos(local_angle) * local_radius
        z = self.z_center + math.sin(local_angle) * local_radius
        self.surface_points.append([x, y, z])

    def get_vector(self, surface_point):
        """
        This function compute vector coordinates between sphere center and a given
        surface points.
        input:
        surface_point                   list    A list of spacial coordinate for a surface point
        return:
        x_vector, y_vector, z_vector    tuple   x, y and z vector coordinates
        """
        x_vector = surface_point[0] - self.x_center 
        y_vector = surface_point[1] - self.y_center 
        z_vector = surface_point[2] - self.z_center  
        return x_vector, y_vector, z_vector
