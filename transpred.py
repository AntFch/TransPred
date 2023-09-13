#! /usr/bin/env python3

#Import module
############################################################################################################

import multiprocessing as mp
import sys
import os
import re
import time
import warnings

try:
	import argparse
except ImportError as exception:
	print("ERROR: Module argparse is not installed. Advice: Advice: Use the transs_prot_pred.yaml file to create a powerfull environnment")

try:
	import numpy as np
except ImportError as exception:
	print("ERROR: Module Numpy is not installed. Advice: Advice: Use the transs_prot_pred.yaml file to create a powerfull environnment")

try:
    from src import protein as prot
    from src import sphere as sph
    from src import membrane as memb
except ImportError as exception:
	print("ERROR: src folder doesn't contain protein?py, membrane.py or sphere.py. Please reinstall the tool with github.")
	
#Define argument
############################################################################################################

parser = argparse.ArgumentParser(prog = "transpred.py",
description = "This program allows to predict membrane plane position from a PDB file.")
parser.add_argument("-i", "--input", type = str, help = "PDB file to parse.")
parser.add_argument("-o", "--output", type = str, help = "Output path directory")
parser.add_argument("-t", "--threshold", type = float, 
help = "Threshold for solvent exposed area (default: 0.30)", default = 0.30)
parser.add_argument("-p", "--point",type = int, 
help = "Test point number (default: 1000)", default = 1000)
parser.add_argument("-c", "--cpu", type = int, 
help = "Number of used CPU (default: 1)", default = 1)
args = parser.parse_args()

#Define function
###########################################################################################################	

def parse_args(input_file, output_path, cpu, threshold):
    """
    This function parses given arguments
    
    input:
        input_file      str      path to input PDB file
        output_path     str      output path
        cpu             int      CPU number
        threshold       float    threshold
    return:
        protein_name    str      protein name
    """
    #Parse input file
    if input_file is None:
        sys.exit("ERROR: Input file must be given.")
    
    if not os.path.exists(input_file):
        sys.exit("ERROR: Input file doesn't exist.")
        
    extension = input_file.split(".")[-1].lower()
    if not extension == "pdb":
        sys.exit("ERROR: Input file is not a PDB file.")
        
    protein_name = input_file.split("/")[-1].split(".")[:-1]
    
    #Parse output
    if output_path is not None:
        if not os.path.exists(output_path):
            sys.exit("ERROR: Output doesn't exist.")
    
    #Parse CPU number
    cpu_count = mp.cpu_count()
    if cpu > cpu_count:
        sys.exit("ERROR: CPU number entered is more greater than your \
        computer ability.")
        
    #Parse threshold
    if threshold < 0 or threshold > 1:
        sys.exit("ERROR: Threshold velue must be between 0 and 1.")     
    return protein_name
    
    
def parallel_process(protein, sphere, nb_cpu):
    """
    This function allow to parallelise test points analysis
    input:
        protein         object  protein object to parse
        sphere          object  sphere object which contains test points
        nb_cpu          int     CPU number to use
    return:
        chunk_results   list    list of results per CPU
    """
    #Define chunk_results which will collect results per CPU
    chunk_args = []
    
    #Define point number per CPU
    nb_point = len(sphere.surface_points)
    point_per_cpu = nb_point // nb_cpu
    remainder = nb_point % nb_cpu
    
    #Define start and end index for surface_points list
    for cpu_number in range(nb_cpu):
        start = cpu_number *point_per_cpu
        end = ((cpu_number + 1)*point_per_cpu) - 1
        #Add remainder for the last cpu
        if cpu_number == nb_cpu - 1 and remainder != 0:
            end += remainder
        args = (sphere, protein, start, end)
        chunk_args.append(args)  
          
    #Make parallel process
    with mp.Pool(nb_cpu) as p:
        chunk_results = p.starmap(process_chunk, chunk_args)
    return chunk_results


def process_chunk(sphere, protein, start, end):
    """
    This function allows to process a chunk
    input:
        sphere          object  sphere which contains test points
        protein         object  protein object which contains residues to parse
        start           int     start surface points index
        end             int     start surface points index
    return:
        chunk_result    list    list with point index, distance, plane parameters
                                and hydrophobic residue number
    """
    #Extract mass center coordinate
    x_center, y_center, z_center = sphere.x_center, sphere.y_center, sphere.z_center
    
    #Define start distance
    START_DISTANCE = 7
    
    #Make loop over points
    chunk_result = []
    for i in list(range(start, end + 1)):
        #Extract point
        point = sphere.surface_points[i]
        
        #Extract point coordinate
        x, y, z = point[0], point[1], point[2]
        
        #Compute vector cooridnate
        x_vector, y_vector, z_vector = sphere.get_vector(point)
        
        #Assign plane coefficient (corresponding to orthogonal vector coordinate)
        a, b, c = x_vector, y_vector, z_vector
        
        #Looking for free coefficient start value, end vlaue and incrementation
        free_coeff_center = -1*(a*x_center + b*y_center + c*z_center)
        free_coeff_surface = -1*(a*x + b*y + c*z)
            
        #Giving free coeff variation corresponding to 1 Angström
        incrementation = (free_coeff_surface - free_coeff_center) / sphere.global_radius
        
        #Generate list of free coeff to test
        free_coeff = free_coeff_center +  START_DISTANCE*incrementation #add 7 ANGSRTÖM to start
        list_free_coeff = []
        while abs(free_coeff) < abs(free_coeff_surface):
            list_free_coeff.append(free_coeff)
            free_coeff += incrementation
        
        #Test all free_coeff over a loop        
        for free_coeff in list_free_coeff:
            #Initiliase parameter
            nb_hydrophobic_exposed = 0
            nb_exposed = 0
            #Compute sign value for mass_center
            plane_value_center = a*x_center + b*y_center + c*z_center 
            sign_center_1 = plane_value_center + free_coeff
            sign_center_2 = plane_value_center - free_coeff
                
            #Test all residue    
            for residue in protein.residues:
                #Looking sign value for mass center point
                #It will give the sens to interpret if a point is between the two planes
                plane_value = a*residue.x + b*residue.y + c*residue.z
                sign_1 = plane_value + free_coeff
                sign_2 = plane_value - free_coeff
                #Looking point location
                if sign_center_1 < 0 and sign_center_2 > 0:
                    if sign_1 < 0 and sign_2 > 0:
                        nb_exposed += 1
                        if residue.is_hydrophobic():
                            nb_hydrophobic_exposed += 1
                else:
                    if sign_1 > 0 and sign_2 < 0:
                        nb_exposed += 1
                        if residue.is_hydrophobic():
                            nb_hydrophobic_exposed += 1
            
            if nb_exposed != 0:
                hydrophobic_ratio = round(nb_hydrophobic_exposed / nb_exposed, ndigits = 4)
            else:
                hydrophobic_ratio = 0
            #Report distance for each d value tested                            
            result = [i, (a, b, c, free_coeff), nb_exposed, nb_hydrophobic_exposed, hydrophobic_ratio]
            chunk_result.append(result)                    
    return chunk_result


def parse_chunk_results(chunk_results):
    """
    This function looks chunk_results in order to find
    best plane parameters
    input:
        chunk_results    list   list of chunk result
    return:
        plane_parameters tuple  best plane parameters
    """
    #Create empty optimal result
    optimal_result = []
    #Create global numpy array
    global_result = []
    for chunk_result in chunk_results:
        if len(chunk_result) > 0:
            for result in chunk_result:
                point = int(result[0])
                a, b, c, d = result[1][0], result[1][1], result[1][2], result[1][3] 
                nb_exposed = int(result[2])
                nb_hydrophobic = int(result[3])             
                ratio = result[4]
                to_add = [point, a, b, c, d, nb_exposed, nb_hydrophobic, ratio]
                global_result.append(to_add)

    global_result = np.array(global_result)
    best_value = max(global_result[:,7])
    
    #Reduce global_result according to best value
    to_store = np.in1d(global_result[:,7], np.array([best_value]))
    global_result = global_result[to_store]
    
    #Create best results list
    ##Extract point list
    point_list = list(global_result[:,0])
    ##Remove duplicate
    point_list = list(dict.fromkeys(point_list))
    for point in point_list:
        to_store = np.in1d(global_result[:,0], np.array([point]))
        point_result = global_result[to_store]
        point_result = point_result[0, ]
        optimal_result.append(point_result)
    return optimal_result



def compute_run_time(start_time, end_time):
    """
    This function compute run time of the process
    input:
        start_time  float   start time in seconds
        end_time    float   end time in seconds
    return:
    run_time        str     string which give run time in min:s:cs
    """
    #Get total time in seconds
    tot_time = end_time - start_time
    
    #Compute minute needed
    minute, remainder = int(tot_time // 60), tot_time % 60
    
    #Update tot_time
    tot_time = remainder
    
    #Compute seconds
    second, remainder = int(tot_time // 1), tot_time % 1
    
    tot_time = remainder
    
    #Compute centiseconds
    milisecond = round(tot_time, ndigits = 3)
    
    run_time = f"{minute} min:{second} s::{milisecond}ms" 
    return run_time
    
 #Main program
###########################################################################################################

if __name__ == "__main__":
    
    start_time = time.time()
    warnings.filterwarnings('ignore')
        
    #Parse argument
    protein_name = parse_args(args.input, args.output, args.cpu, args.threshold)
    print(f"Compute membrane plane for {args.input}")
    
    #Parse PDB file
    protein = prot.parse_pdb(args.input, protein_name)
    #max_coordinate, axis = protein.get_max_coordinate()
    print("### PDB file parsed ###")
    
    #Compute mass center of non filtered protein
    x_center, y_center, z_center = protein.compute_mass_center()
    max_distance = protein.get_most_distant_residue([x_center, y_center, z_center])
    print("### Mass center calculated ###")
    
    #Focus on exposed residues
    protein_filter = prot.filter_solvent_accessible_area(args.input, protein, args.threshold)
    print("### Solvent exposed residues got ###")
    
    #Build test points
    sphere = sph.Sphere(x_center, y_center, z_center, max_distance)
    
    point_number = int(args.point / 2)
    sphere.get_vertical_incrementation(point_number)
    
    
    for i in range(point_number + 1):
        sphere.add_fibonacci_point(i)
    
    print("### Test points created ###")
    
    #Process surface points
     
    chunk_results = parallel_process(protein_filter, sphere, args.cpu)
    
    optimal_results = parse_chunk_results(chunk_results)
    
    print("### Test points analysed ###")
    print(f"Found {len(optimal_results)} optimal membrane location for your protein.")
    for i, optimal_result in enumerate(optimal_results):
        plane_coeff = np.array(optimal_result[1:5])
        plane_coeff = np.round(plane_coeff, 3)
        print(f"Membrane {i + 1}:")
        print(f"Plane 1: {plane_coeff[0]} x + {plane_coeff[1]} y + {plane_coeff[2]} z + {plane_coeff[3]} = 0")
        print(f"Plane 2: {plane_coeff[0]} x + {plane_coeff[1]} y + {plane_coeff[2]} z + {-1*plane_coeff[3]} = 0")
    
        #Build membrane
        mass_center = [x_center, y_center, z_center]
        plane_coeff = optimal_results[i][1:5]
        membrane = memb.Membrane(1, sphere.global_radius, plane_coeff, mass_center)
        membrane.get_points()
    
        #Create PDB file
        if args.output is None:
            output_pdb = f"{protein_name[0]}_membrane_{i + 1}.pdb"
        else: 
            output_pdb = os.path.join(args.output, f"{protein_name[0]}_membrane_{i + 1}.pdb")    
        memb.create_pdb_file(args.input, output_pdb, membrane)

    print("### Output generated ###")

    #Compute run time
    end_time = time.time()
    run_time = compute_run_time(start_time, end_time)
    print(f"Run time: {run_time}")