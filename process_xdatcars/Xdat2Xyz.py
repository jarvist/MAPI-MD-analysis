import os
#import asap3
#from asap3.analysis.rdf import RadialDistributionFunction
import ase.io as io
import numpy as np
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-p", "--poscar-file",
                  action="store", type="string", dest="poscarfile", default="POSCAR",
                  help="The POSCAR file to take coords from. [default = POSCAR]")
parser.add_option("-f", "--trajectory-file",
                  action="store", type="string", dest="filename", default="XDATCAR",
                  help="The XDATCAR file to take coords from. [default = XDATCAR]")
parser.add_option("-o", "--output-file",
                  action="store", type="string", dest="outputfile", default="Trajectory.xyz",
                  help="The file to write to. [default = Trajectory.xyz]")
(options, args) = parser.parse_args()


atoms=io.read(options.poscarfile,format="vasp")
natoms = len(atoms)

if isinstance(options.filename, str):
    f = open(options.filename)
else: # Assume it's a file-like object
    f = options.filename

data=f.readlines()

nimages=(len(data)-5)/(natoms+1)
#outfile = open(options.outputfile, "a")
images = []
for i in range(nimages):
    images.append(atoms.copy())
    pos=np.zeros((natoms,3),np.float)
    for j in range(natoms):
        pos[j,0]=float(data[8+i*(natoms+1)+j].split()[0])
        pos[j,1]=float(data[8+i*(natoms+1)+j].split()[1])
        pos[j,2]=float(data[8+i*(natoms+1)+j].split()[2])
    images[i].set_scaled_positions(pos)
    
io.write(options.outputfile,images,format="xyz")
