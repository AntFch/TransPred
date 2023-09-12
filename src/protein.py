"""
Module protein

This module define protein class with its attributs and associated methods. It also defines residue class;
"""



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
