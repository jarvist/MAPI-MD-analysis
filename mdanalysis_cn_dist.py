# mdanalysis_cn_dist.py - Jarvist Moore Frost April 2014
#  Scripts to analyise CN axis orientation from Ab-Initio MD on MAPI

#When using MDAnalysis in published work, please cite doi:10.1002/jcc.21787
import MDAnalysis
import MDAnalysis.analysis

import matplotlib.pyplot as plt

import numpy
import math
from IPython import embed #iPython magic for interactive session...

#universe item; contains all the mdanalysis stuff
u= MDAnalysis.Universe("1stframe.pdb","Trajectory.xyz")#MAPI_222_equilibr_traj.xyz")

# Above .XYZ file was generated with Keith Butler's Xdat2Xyz.py (depends on ASE)
# The pdb file was generated with Open Babel

#embed()

GenThetas=False    # Theta / Phis to STDOUT for plotting externally
GenXYZ=True        # .XYZ file to STDOUT of CN axes, for Pymol plotting
ExploitSymm=True  # Exploit full symmetry = 42 
Exploit8fold=False # Exploit 8-fold symmetry

thetas=[] # List to collect data for later histogramming
phis=[]

carbons=u.selectAtoms('name C')
nitrogens=u.selectAtoms('name N')

#Hard wired cubic PBCs - the MDAnalysis readers don't seem to pick these up from PDB / XYZ files (?)
mybox=numpy.array([12.5801704656605438,12.5477821243936827,12.5940169967174249],dtype=numpy.float32)
imybox=1/mybox

if GenThetas:
    print "# box: ",mybox, imybox

def spherical_coordinates(cn):
    cn=sorted(abs(cn),reverse=True)
    l=numpy.linalg.norm(cn)
    x=cn[0]
    y=cn[1]
    z=cn[2]
    theta = math.acos(z/l)
    phi   = math.atan2(y,x)
    return (theta,phi)

# Checking values returned by three symmetric alignments
print "[1,0,0] ", spherical_coordinates(numpy.array([1,0,0]))
print "[1,1,0] ", spherical_coordinates(numpy.array([1,1,0]))
print "[1,1,1] ", spherical_coordinates(numpy.array([1,1,1]))
#embed()

for ts in u.trajectory:
    # Of course - this list doesn't actually change between frames; it's part of the topology
    r=MDAnalysis.analysis.distances.distance_array(carbons.coordinates(),nitrogens.coordinates(),mybox)
    # OK, this is now the distance array of all Ns vs. Cs, considering PBCs
#    print r
    if GenXYZ:
        print "16\n" #header for .xyz frames
    for carbon, nitrogenlist in enumerate(r):
        for nitrogen,distance in enumerate(nitrogenlist):
            if distance<1.6: 
                #print carbon, nitrogen, distance
                cn=carbons[carbon].pos - nitrogens[nitrogen].pos
                # Filthily hack in some orthorhombic (only) PBCs as I can't see how MDanalysis does this properly
                for dim in 0,1,2:
                    s=imybox[dim] * cn[dim] # Minimum image convetion, algo from MDAnalysis calc_distances.h
                    cn[dim]=mybox[dim]*(s-round(s)) 
 #               print d

                #OK; cn is now our 3-vector along CN bond
                if Exploit8fold:
                    cn=abs(cn)
                if ExploitSymm:
                    cn=sorted(abs(cn),reverse=True) #Exploit Oh symmetry - see workings in Jarv's notebooks
                # (rear page of Orange book, 16-4-14)
                #print cn,r

#                x=cn[0]
#                y=cn[1]
#                z=cn[2]
#                l=numpy.linalg.norm(cn)
#                theta =  math.acos(z/l)
#                phi   = math.atan2(y,x)
#                print "old inline code...",theta,phi
                theta,phi = spherical_coordinates(numpy.array(cn))
#                print "via spherical_coordinates fn...",theta,phi

                thetas.append(theta) #append this data point to lists
                phis.append(phi)

                #OK, now we fold along theta, phi, to account for symmetry (TODO: Check!)
                if GenThetas:
                    print "%f %f %f %f %f %f %f %f" %(theta,phi,x,y,z,d[0],d[1],d[2])
                if GenXYZ:
                    # quick and dirty .xyz output of animated CN axis
                    cx=carbon%3*0.5*mybox[0] + 1.5*mybox[0]
                    cy=math.floor(carbon/3)*0.5*mybox[0]
                    cz=0.0
                    print "  C %10.5f %10.5f %10.5f" %(cx,cy,cz) #'carbon' as offset
    #                print "N %f %f %f" %(cx+x,cy+y,cz+z) #With +x+y+z --> reduced form
                    print "  N %10.5f %10.5f %10.5f" %(cx-cn[0],cy-cn[1],cz-cn[2]) #With +x+y+z --> reduced form

H,xedges,yedges = numpy.histogram2d(thetas,phis,bins=20)
H.shape, xedges.shape, yedges.shape
extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]
plt.imshow(H,extent=extent,interpolation='nearest')
plt.colorbar()
plt.show()

end

for ts in u.trajectory:
    for c in carbons:
        for n in nitrogens:
            r=c.pos - n.pos
            print r
            #embed()
            print c,n,r
            d=numpy.linalg.norm(r)
            if (d<1.8): #i.e. these are near bonds. TODO: PBCs _NOT_ treated properly.
                #print ts.frame,d,r
                #OK; r is now our 3-vector along CN bond
                cn=sorted(abs(r),reverse=True) #Exploit Oh symmetry - see workings in Jarv's notebooks
                # (rear page of Orange book, 16-4-14)
                #print cn,r

                x=cn[0]
                y=cn[1]
                z=cn[2]

                theta = math.acos(z/d)
                phi   = math.atan2(y,x)
                #OK, now we fold along theta, phi, to account for symmetry (TODO: Check!)
                #theta = theta % math.pi/4
                #phi   = phi % math.pi/4
                print "%f %f %f %f %f" %(theta,phi,x,y,z)
