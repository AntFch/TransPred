#! /usr/bin/env python3

"""
Module protein

This module define protein class with it attributs and associated methods. It also defines residue class.
"""
#Import module
##########################################################################################################

import math
import re

try:
    import numpy as np
    from Bio.PDB import PDBParser
    from Bio.PDB.DSSP import DSSP
except ImportError as exception:
	print("ERROR: Numpy or Biopython module is not installed. Advice: Use the transs_prot_pred.yaml file to create a powerfull environnment")

#Define classes
###########################################################################################################

class Protein:
    def __init__(self, prot_name = None):
        """
        This function initialises a protein object by assign a name and 
        creates an empty residues list
        input:
        prot_name   str     protein name
        residues    list    list of selected residues which belong to the protein 
        """
        self.prot_name = prot_name
        self.residues = []

    def compute_mass_center(self):
        """
        This function compute mass center for the protein from residues coordinates. 
        return:
        x_center    float  x mass center coordinate
        y_center    float  y mass center coordinate
        z_center    float  z mass center coordinate
        """
        x_coor = []
        y_coor = []
        z_coor = []
        for residue in self.residues:
            x_coor.append(residue.x)
            y_coor.append(residue.y)
            z_coor.append(residue.z)
        x_center = np.mean(np.array(x_coor))
        y_center = np.mean(np.array(y_coor))
        z_center = np.mean(np.array(z_coor))
        return x_center, y_center, z_center

    def add_residue(self, residue):
        """
        This function allows to add a residue to a protein
        """
        self.residues.append(residue)

    def get_residue(self, chain, position):
        """
        This function allows to extract a specific residue from a given chain and position
        input:
        chain       str     a given protein chain to search
        position    int     a given protein chain position to search
        """
        for i in range(len(self.residues)):
            residue = self.residues[i]
            if residue.chain == chain and residue.position == position:
                break
        return residue

    def get_most_distant_residue(self, mass_center):
        """
        This function computes the most distant residue for a protein
        input:
            mass_center     list   Mass center coordinates (x, y, z)
        output:
            max_distance    float   max observed distance
        """
        distance_list = []
        for residue in self.residues:
            #Compute distance for each residue
            distance = math.sqrt((residue.x - mass_center[0])**2 + (residue.y - mass_center[1])**2 + (residue.z - mass_center[2])**2)
            distance_list.append(distance)
        max_distance = max(distance_list)
        return max_distance


    def get_max_coordinate(self):
        """
        This function allows ro extract max coordinate observed
        """
        max_coordinate = 0
        for residue in self.residues:
            if residue.x > max_coordinate:
                max_coordinate = residue.x
                axis = "x"
            if residue.y > max_coordinate:
                max_coordinate = residue.y
                axis ="y"
            if residue.z > max_coordinate:
                max_coordinate = residue.z
                axis = "z"
        return max_coordinate, axis


class Residue:
    """
    This class define some residual properties
    """
    def __init__(self, x = None, y = None, z = None, atom_type = "CA",
    res_name = None, chain = None, position = None):
        """
        Assign some attributs for a residue
        
        input:
        x           float   x coordinates
        y           float   y coordinates
        z           float   z coordinates
        atom_type   str     atom type (default: CA)
        res_name    str     residue name
        position    int     position in the sequence
        chain       str     chain name
        """
        self.x = x
        self.y = y
        self.z = z
        self.atom_type = atom_type
        self.res_name = res_name
        self.position = int(position)
        self.chain = chain

    def is_hydrophobic(self):
        """
        This function define if a residue is hydrophobic
        """
        if self.res_name in ["PHE", "GLY", "LEU", "ILE", "MET", "VAL", "TRP", "TYR"]:
            return True
        else:
            return False
            
    def is_solvent_exposed(self, aas = None, threshold = None):
        """
        This function define if a residue is exposed to the solvent
        
        input:
        aas         float   relative part of residue area exposed to the solvent
        threshold   float   cutoff value to consider a residue exposed or not
        """
        if aas >= threshold:
            return  True
        else:
            return False

#Define functions
###########################################################################################################

def parse_pdb(pdb_file, protein_name):
    """
    This function allows to parse a given PDB file in order to extract 
    carbon alpha information
    
    input:
        pdb_file    str     path to PDB file to parse
    return
        protein     object  A protein object
    """
    #Define empty residues list
    protein = Protein(protein_name)
    
    #Open pdb file
    with open(pdb_file, "r") as input_pdb:
        for line in input_pdb:
            #Select ATOM line with carbon alpha atom
            if line[0:6].find("ATOM") >= 0:
                if line[12:16].find("CA") >= 0:
                    #Extract coordinates
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    #Extract residue name
                    res_name = line[17:20].strip().upper()
                    #Extract position in the sequence
                    position = float(line[22:26].strip())
                    #extract chain name
                    chain = line[21:22]
                    #Define residue instance
                    residue = Residue(x, y, z, "CA", res_name, chain, position)
                    protein.add_residue(residue)
    return protein


def filter_solvent_accessible_area(pdb_file, protein, threshold):
	"""
	This function computes solvent accessible area from a pdb file using DSSP program 
	with a Biopython API and filters exposed residue according to a threshold
	
	input:
		pdb_file			str		path to PDB file to analyse
		protein             object  A protein object to filter
		threshold           folat   A threshold between 0 and 1 to filter aas values
	return:
		protein_filter	    object	Filtered protein object
	"""
	#Create new filtered protein
	protein_name = protein.prot_name
	protein_filter = Protein(protein_name)
	#Parse PDB file
	p = PDBParser()
	structure = p.get_structure("prot", pdb_file)
	model = structure[0]
	#Use DSSP API to compute solvent accessible area
	dssp = DSSP(model, pdb_file, dssp = "mkdssp")
	#Extract keys to access DSSP results
	residues_keys = list(dssp.keys())
	#Extract for each residue solvent accessible area
	for i in range(len(residues_keys)):
	    #Get relative solvent exposed area and rount it to 4 digits
	    key = residues_keys[i]
	    chain = key[0]
	    position = int(re.findall('\d+', str(key[1]))[0])
	    aas = round(dssp[key][3], ndigits = 4)
	    residue = protein.get_residue(chain, position)
	    is_exposed = residue.is_solvent_exposed(aas, threshold)
	    if is_exposed:
	        protein_filter.add_residue(residue)
	return protein_filter                        
