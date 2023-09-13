# TransPred

![TransPred](src/main_image.png)

Welcome to TransPred ! This tool written in Python language can predict the position of membrane for a supposed transmembrane protein. It takes in input a PDB file from experimental strcuture prediction and gives in output optimal membrane plane equations and modified PDB file which contains dummy atoms coordinates that represent predicted membrane planes.

This program was done for a university project, based on the following paper:
> [Reference](https://pubmed.ncbi.nlm.nih.gov/15180935/) :
Tusnády GE, Dosztányi Z, Simon I. Transmembrane proteins in the Protein
Data Bank: identification and classification. Bioinformatics. 2004 Nov
22;20(17):2964-72. Epub 2004 Jun 4. PubMed PMID: 15180935.


## Dependecies
### Python packages:
- argparse (tested with version 1.4.0)
- Biopython (tested with version 1.81)
- Numpy (tested with version 1.25.2)

### External program:
- DSSP (tested with version 3.0.0)

## Installation
This installation requires to have miniconda or anaconda installed. If you don't one of this please visit the following web site to install it: [conda installation](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)

### Download repository

**Option 1: using the HTTPS link**

    git clone https://github.com/AntFch/TransPred.git

**Option 2: using the SSH key**
This option suppose that you already provided a ssh key to github.

    git clone git@github.com:AntFch/TransPred.git

### Setup your environment
src folder contains in TransPred folder, a YAML file that allows to create an appropriate conda environment. Execute the following command to get it.

    cd TransPred
    conda create env -f src/transpred.yaml
    conda activate TransPred

## Usage
To see the help, run the following command:

    python transpred.py --help

You will obtain this output:

    usage: transpred.py [-h] [-i INPUT] [-o OUTPUT] [-t THRESHOLD] [-p POINT] [-c CPU]

    This program allows to predict membrane plane position from a PDB file.

    optional arguments:
      -h, --help            show this help message and exit
      -i INPUT, --input INPUT
                        PDB file to parse.
      -o OUTPUT, --output OUTPUT
                        Output path directory
      -t THRESHOLD, --threshold THRESHOLD
                        Threshold for solvent exposed area (default: 0.30)
      -p POINT, --point POINT
                        Test point number (default: 1000)
      -c CPU, --cpu CPU     Number of used CPU (default: 1)

Differents arguments that can be used are:
- INPUT: a keyword argument which provides path to PDB file to parse.
- OUTPUT: an optional argument which gives path to output folder
- THRESHOLD: an optional value to set threshold to consider a resiude exposed from AAS value providing by DSSP program. This value is between 0 and 1. Default value is set to 0.30.
- POINT: an optional value that defines number of points to test. Default value is set to 1000.
- CPU: an optionnal value that define the number of CPU to use.

## Example usage

Now, you can run your first membrane position prediction. The folder data contains the follwing file: **atp_synthase_chainA_1h68.pdb** extracted from PDBTM database (ID: 1h68). It's about the chain A of ATP synthase located in the mitochondrial membrane. It contains 219 residues. To make the membrane planes prediction run the following command.

    python transpred.py -i /data/atp_synthase_chainA_1h68.pdb

You will get the following outputs:
    
    Compute membrane plane for data/atp_synthase_chainA_1h68.pdb
    ### PDB file parsed ###
    ### Mass center calculated ###
    ### Solvent exposed residues got ###
    ### Test points created ###
    ### Test points analysed ###
    Found 1 optimal membrane location for your protein.
    Membrane 1:
    Plane 1: -1.401 x + 5.556 y + 26.092 z + -197.038 = 0
    Plane 2: -1.401 x + 5.556 y + 26.092 z + 197.038 = 0
    ### Output generated ###
    Run time: 0 min:1 s::0.830ms

Note that run time information depends on your computational power.
The program gives the realised step during execution. It also gives membrane planes equations.

You will find the following file in your current directorie: **atp_synthase_chainA_1h68_membrane_1.pdb**. This file can also be found in data folder after if you don't want to run the command line. It can be visualised in PyMOL. If you don't have it, you can visit the following web site to do the installation: [PyMOL installation](https://pymol.org/dokuwiki/doku.php?id=installation).
With PyMol we get this:
![TransPred](src/main_image.png)

## Acknowledgements
Sincere thanks to Gelly J.C. (Phd and professor associate at Paris-Cité University) and Galochkina T. (Phd at Paris-Cité University) for the lead of this project.