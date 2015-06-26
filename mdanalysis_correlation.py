# mdanalysis_cn_dist.py - Jarvist Moore Frost April 2014
#  Scripts to analyise CN axis orientation from Ab-Initio MD on MAPI

#When using MDAnalysis in published work, please cite doi:10.1002/jcc.21787
import MDAnalysis
import MDAnalysis.analysis

import matplotlib.pyplot as plt

import numpy
import math
import random

#from IPython import embed #iPython magic for interactive session...
import datetime # current date for log files etc.
now=datetime.datetime.now().strftime("%Y-%m-%d-%Hh%Mm") #String of standardised year-leading time

import sys

if (len(sys.argv)>1):
    trajfilename=sys.argv[1]
else:
    trajfilename="300K.xyz"

print 'Argument List:', str(sys.argv), " Trajectory filename: ", trajfilename

#universe item; contains all the mdanalysis stuff
u= MDAnalysis.Universe("300K.pdb",trajfilename) #Trajectory.xyz")#MAPI_222_equilibr_traj.xyz")

# Above .XYZ file was generated with Keith Butler's Xdat2Xyz.py (depends on ASE)
# The pdb file was generated with Open Babel

#embed()

GenThetas=False    # Theta / Phis to STDOUT for plotting externally
GenXYZ=False        # .XYZ file to STDOUT of CN axes, for Pymol plotting
ExploitSymm=False  # Exploit full symmetry = 48 fold
Exploit8fold=False # Exploit 8-fold symmetry

thetas=[] # List to collect data for later histogramming
phis=[]

dotcount=[0.,0.,0.]

carbons=u.selectAtoms('name C')
nitrogens=u.selectAtoms('name H') # Formadinium - look along C-H axis

# FAPI 2x2x2 supercell
#  1.0000000000000000
#    12.8149204254000004    0.0000000000000000    0.0000000000000000
#    0.0000016686000000   12.5312404631999996    0.0000000000000000
#    0.0000000000000000    0.0000000000000000   12.6751403808000003

#Hard wired cubic PBCs - the MDAnalysis readers don't seem to pick these up from PDB / XYZ files (?)
mybox=numpy.array([12.8149204254000004,12.5312404631999996,12.6751403808000003],dtype=numpy.float32)
imybox=1/mybox

if GenThetas:
    print "# box: ",mybox, imybox

def partition_alignment(cn):
# Dot product of cn vs. [1,0,0],[1,1,0],[1,1,1]
#print "cn: ", cn, "theta: ", theta, "phi:", phi
#print "[1,0,0]", numpy.dot([1,0,0],cn)
#print "[1,1,0]", numpy.dot([math.sqrt(1./2),math.sqrt(1./2),0],cn)
#print "[1,1,1]", numpy.dot([math.sqrt(1./3),math.sqrt(1./3),math.sqrt(1./3)],cn)
    dots=[]
    
    cn=cn/numpy.linalg.norm(cn) #normalise
    sm=sorted(abs(cn),reverse=True)

    for unit in [1,0,0],[1,1,0],[1,1,1]:
        unit=unit/numpy.linalg.norm(unit)
#       print "dot:", unit,numpy.dot(unit,cn),theta,phi
        dots.append(numpy.dot(unit,sm))
    dotcount[numpy.argmax(dots)]=dotcount[numpy.argmax(dots)]+1
#    print dotcount
    print cn[0],cn[1],cn[2],numpy.argmax(dots) #x,y,z, which type (for colour)
#    print "FaceTheta: ", math.acos(min(dots[0],1.0)) #fix to domain errors as dots[] creeping above 1.0
#    print "EdgeTheta: ", math.acos(min(dots[1],1.0))
#    print "DiagTheta: ", math.acos(min(dots[2],1.0))

def test_partition_alignment():
    for i in range(100000):
        sph = numpy.array([0.,0.,0.]) 
        # via http://www.astro.uni-bonn.de/~mmetz/py/rndsphere.php
        #  - though code based on the Mathworld example
        sph[2] = random.uniform(-1.0,1.0)
        z2     = math.sqrt(1.0 - sph[2]*sph[2])
        phi    = (2. * math.pi) * random.random()
        sph[0] = z2 * math.cos(phi)
        sph[1] = z2 * math.sin(phi)
        # apply symm...
        #print sph
        #sph=sorted(abs(sph),reverse=True)
        #print sph
        partition_alignment(sph)

    print dotcount

def spherical_coordinates(cn):
#    cn=sorted(abs(cn),reverse=True) #This forces symmetry exploitation. Used for figuring out what [x,y,z] corresponds to which point in the figure
    l=numpy.linalg.norm(cn)
    x=cn[0]
    y=cn[1]
    z=cn[2]
    theta = math.acos(z/l)
    phi   = math.atan2(y,x)
    return (theta,phi)

print ("Checking symmetry reduced values returned by three symmetric alignments")
print "[1,0,0] ", spherical_coordinates(numpy.array([1,0,0]))
print "[1,1,0] ", spherical_coordinates(numpy.array([1,1,0]))
print "[1,1,1] ", spherical_coordinates(numpy.array([1,1,1]))
#embed()

#test_partition_alignment()
#end
cnsn=numpy.zeros(shape=(u.trajectory.numframes,8,3))

print("Loading Trajectory...")
for ts in u.trajectory:
    # Of course - this list doesn't actually change between frames; it's part of the topology
    r=MDAnalysis.analysis.distances.distance_array(carbons.coordinates(),nitrogens.coordinates(),mybox)
    # OK, this is now the distance array of all Ns vs. Cs, considering PBCs
#    print r
    if GenXYZ:
        print "16\n" #header for .xyz frames; number of atoms per frame
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
                x=cn[0] #stash vector before we do symmetry reduction, for output
                y=cn[1]
                z=cn[2]
#
                #OK; cn is now our 3-vector along CN bond
                if Exploit8fold:
                    cn=abs(cn)
                if ExploitSymm:
                    cn=sorted(abs(cn),reverse=True) #Exploit Oh symmetry - see workings in Jarv's notebooks
                # (rear page of Orange book, 16-4-14)
                
                #print cn,r
                cn=cn/numpy.linalg.norm(cn) #normalise

#                print ts.frame,carbon,cn
                cnsn[ts.frame-1,carbon]=cn #-1 to array index from 0


print cnsn

print("OK, loaded Trajectory and generated alignment vectors. Now calculating correlations as a function of time.")

T=1000 #time steps over which to calculate correlation
correlation=numpy.zeros(T)

for carbon,nitrogenlist in enumerate(r): #iterate over all MA ions
    print("Ions: %d Calculating correlation."%(carbon))
    for i in xrange(u.trajectory.numframes-T): #over number of frames - time to collect data
        for dt in xrange(T): #
#        print i, t
            correlation[dt]=correlation[dt]+numpy.dot(cnsn[i,carbon],cnsn[i+dt,carbon]) 
            #Here comes the science bit, concentrate!
correlation=correlation/((u.trajectory.numframes-T)*8.0) # hard coded number of MA ions

print correlation

fig=plt.figure()
ax=fig.add_subplot(111)

plt.plot(correlation)

ax.set_title("Dot Product correlation of FA-vector")
ax.set_xlabel(r'$\Delta t$, frames')
ax.set_ylabel(r'$r_{T}.r_{T+\Delta t}$')

plt.show()
fig.savefig("%s-mdanalysis_FAPI_correlation_averages.pdf"%now)
fig.savefig("%s-mdanalysis_FAPI_correlation_averages.png"%now)

print("OK; calculating autocorrelation, time to cross zero.")
print ("   alternative analysis, time till cns[i].cns[i+dt] < 0.0 (i.e. falls off to 90 degrees)")

timetozero=[] #list of times to fall to zero

for carbon,nitrogenlist in enumerate(r):
    print("Ion: %d Calculating dot-product correlation."%(carbon))
    for i in xrange(u.trajectory.numframes):
        for dt in xrange (u.trajectory.numframes-i):
#            print i,dt,numpy.dot(cns[i],cns[i+dt])
            #if (numpy.dot(cns[i],cns[i+dt])<0.0):
                #print i,dt
            #    timetozero.append(dt)
            #    break
#            print cnsn[i,carbon],cnsn[i+dt,carbon],numpy.dot(cnsn[i,carbon],cnsn[i+dt,carbon])
            if (numpy.dot(cnsn[i,carbon],cnsn[i+dt,carbon])<0.0):
                timetozero.append(dt)
                break
print "Total frames, ", u.trajectory.numframes
print "Sum time to zero, ", sum(timetozero)
print "Length time to zero, ", len(timetozero)
print "Average number of frames: \n", (sum(timetozero)/len(timetozero))

fig=plt.figure()
ax=fig.add_subplot(111)

plt.hist(timetozero,100,histtype='stepfilled') # we have quite a lot of data here; so 100 bins

ax.set_title("Dot Product correlation of FA-vector, time to cross zero")
ax.set_xlabel(r'$\Delta t$, frames')
ax.set_ylabel(r'$r_{T}.r_{T+\Delta t}<0.0 ?$')

plt.show()
fig.savefig("%s-mdanalysis_FAPI_correlation_timetocrosszero.pdf"%now)
fig.savefig("%s-mdanalysis_FAPI_correlation_timetocrosszero.png"%now)
