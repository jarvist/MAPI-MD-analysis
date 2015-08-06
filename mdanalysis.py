# mdanalysis_cn_dist.py - Jarvist Moore Frost April 2014 -- June 2015
#  Scripts to analyise CN axis orientation from Ab-Initio MD on MAPI
#  Updated & generalised to analyse FAPI June 2015

#When using MDAnalysis in published work, please cite doi:10.1002/jcc.21787
import MDAnalysis
import MDAnalysis.analysis

import matplotlib.pyplot as plt

import numpy
import math
import random

import sys
if (len(sys.argv)>1):
    trajfilename=sys.argv[1]
else:
    trajfilename="300K.xyz"

print 'Argument List:', str(sys.argv), " Trajectory filename: ", trajfilename

#from IPython import embed #iPython magic for interactive session...
import datetime # current date for log files etc.
now=datetime.datetime.now().strftime("%Y-%m-%d-%Hh%Mm") #String of standardised year-leading time
fileprefix= now + "-" + trajfilename.rsplit( ".", 1 )[ 0 ] + "-mdanalysis"  # regexp to get filename.ext --> filename

print 'Filename prefix: ', fileprefix

#universe item; contains all the mdanalysis stuff
u= MDAnalysis.Universe("300K.pdb",trajfilename) #Trajectory.xyz")#MAPI_222_equilibr_traj.xyz")

# Above .XYZ file was generated with Keith Butler's Xdat2Xyz.py (depends on ASE)
# The pdb file was generated with Open Babel

#embed()

GenThetas=False    # Theta / Phis to STDOUT for plotting externally
GenXYZ=False        # .XYZ file to STDOUT of CN axes, for Pymol plotting
#ExploitSymm=False  # Exploit full symmetry = 48 fold
#Exploit8fold=False # Exploit 8-fold symmetry
DisplayFigures=False # interrupt program with matplotlib, or just silent write to file?

MDTimestep=0.025 # Number of ps per _FRAME_ of data supplied. 
                 # For me, this is MD with dt=0.5 fs ; saving every 50th frame --> 25 fs --> 0.025 ps
                 # Just used to scale outputs / histograms

dotcount=[0.,0.,0.]

carbons=u.selectAtoms('name C')
nitrogens=u.selectAtoms('name H') 
# Formamidinium - look along C-H axis
# Methylammonium - only axis is C-N

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
#    print cn[0],cn[1],cn[2],numpy.argmax(dots) #x,y,z, which type (for colour)
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

    #print dotcount

def spherical_coordinates(cn):
#    cn=sorted(abs(cn),reverse=True) #This forces symmetry exploitation. Used for figuring out what [x,y,z] corresponds to which point in the figure
    l=numpy.linalg.norm(cn)
    x=cn[0]
    y=cn[1]
    z=cn[2]
    theta = math.acos(z/l)
    phi   = math.atan2(y,x)
    theta   = theta - math.pi/2 #to agree with Matplotlib view of angles...
    return (theta,phi)

print ("Checking symmetry reduced values returned by three symmetric alignments")
print "[1,0,0] ", spherical_coordinates(numpy.array([1,0,0]))
print "[1,1,0] ", spherical_coordinates(numpy.array([1,1,0]))
print "[1,1,1] ", spherical_coordinates(numpy.array([1,1,1]))
#embed()

#test_partition_alignment()
#end
cnsn=numpy.zeros(shape=(u.trajectory.numframes,8,3))

print("Loading Trajectory... ['.' = 100 frames]")
for ts in u.trajectory:
    if(ts.frame%100==0):
        sys.stdout.write('.') 
        sys.stdout.flush()
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
               
                #print cn,r
                cn=cn/numpy.linalg.norm(cn) #normalise; needed for orientation work
#                print ts.frame,carbon,cn
                cnsn[ts.frame-1,carbon]=cn #-1 to array index from 0; THIS STORES VALUES FOR CORRELATION
                
                partition_alignment(cn) # Can anyone remember what this does?!

                #Optional output routines...
                if GenThetas:
                    print "%f %f %f %f %f %f %f %f" %(theta,phi,cn[0],cn[1],cn[2],x,y,z)
                if GenXYZ:
                    # quick and dirty .xyz output of animated CN axis
                    cx=carbon%3*0.5*mybox[0] + 1.5*mybox[0]
                    cy=math.floor(carbon/3)*0.5*mybox[0]
                    cz=0.0
                    print "  C %10.5f %10.5f %10.5f" %(cx,cy,cz) #'carbon' as offset
    #                print "N %f %f %f" %(cx+x,cy+y,cz+z) #With +x+y+z --> reduced form
                    print "  N %10.5f %10.5f %10.5f" %(cx-cn[0],cy-cn[1],cz-cn[2]) #With +x+y+z --> reduced form

print(" OK...")
print "Found ", len(cnsn), " alignment vectors." # Print out what we found...

print("Loaded Trajectory and generated alignment vectors. ")
sys.stdout.flush()

### ANALYSIS code below here; from merged files

def correlation():
    print("Now calculating correlations as a function of time.")
    T=1000 #dataframes over which to calculate correlation
    correlation=numpy.zeros(T)

    for carbon,nitrogenlist in enumerate(r): #iterate over all MA ions
        print("Ions: %d Calculating correlation."%(carbon))
        sys.stdout.flush()
        for i in xrange(u.trajectory.numframes-T): #over number of frames - time to collect data
            for dt in xrange(T): #
#           print i, t
                correlation[dt]=correlation[dt]+numpy.dot(cnsn[i,carbon],cnsn[i+dt,carbon]) 
                #Here comes the science bit, concentrate!
    correlation=correlation/((u.trajectory.numframes-T)*8.0) # hard coded number of MA ions

    #print correlation

    fig=plt.figure()
    ax=fig.add_subplot(111)

    plt.plot(correlation)

    ax.set_title("Dot Product correlation of FA-vector")
    ax.set_xlabel(r'$\Delta t$, frames')
    ax.set_ylabel(r'$r_{T}.r_{T+\Delta t}$')

    if (DisplayFigures):
        plt.show()
    fig.savefig("%s-correlation_averages.pdf"%fileprefix)
    fig.savefig("%s-correlation_averages.png"%fileprefix)
    
    # Save data for future plotting; in units rescaled to ps
    f = open("%s-correlation_averages.txt"%fileprefix, "w")
    for frame,count in enumerate(correlation):
        f.write( str(frame*MDTimestep) + ' ' + str(count) + '\n'  )
    f.close()

    print("OK; calculating autocorrelation, time to cross zero.")
    print ("   alternative analysis, time till cns[i].cns[i+dt] < 0.0 (i.e. falls off to 90 degrees)")
    sys.stdout.flush()

    timetozero=[] #list of times to fall to zero

    for carbon,nitrogenlist in enumerate(r):
        print("Ion: %d Calculating dot-product correlation."%(carbon))
        sys.stdout.flush()
        for i in xrange(u.trajectory.numframes):
            for dt in xrange (u.trajectory.numframes-i):
#                print i,dt,numpy.dot(cns[i],cns[i+dt])
#            print cnsn[i,carbon],cnsn[i+dt,carbon],numpy.dot(cnsn[i,carbon],cnsn[i+dt,carbon])
                if (numpy.dot(cnsn[i,carbon],cnsn[i+dt,carbon])<0.0):
                    timetozero.append(dt)
                    break
    print "Total frames, ", u.trajectory.numframes
    print "Sum time to zero, ", sum(timetozero)
    print "Length time to zero, ", len(timetozero)
    print "Average number of frames: \n", (sum(timetozero)/len(timetozero))
    sys.stdout.flush()

    fig=plt.figure()

    ax=fig.add_subplot(111)

    plt.hist(timetozero,100,histtype='stepfilled') # we have quite a lot of data here; so 100 bins

    ax.set_title("Dot Product correlation of FA-vector, time to cross zero")
    ax.set_xlabel(r'$\Delta t$, frames')
    ax.set_ylabel(r'$r_{T}.r_{T+\Delta t}<0.0 ?$')

    if (DisplayFigures):
        plt.show()
    fig.savefig("%s-correlation_timetocrosszero.pdf"%fileprefix)
    fig.savefig("%s-correlation_timetocrosszero.png"%fileprefix)

    # Save data for future plotting; in units rescaled to ps
    f = open("%s-correlation_timetocrosszero.txt"%fileprefix, "w")
    for t in timetozero:
        f.write(str(t*MDTimestep) + '\n')
    f.close()

### OTHER FILE ###

def orientation_density():
    
    thetas=[] # List to collect data for later histogramming
    phis=[]

    thetasOh=[]
    phisOh=[]

    for frame in cnsn[:,:]:
        for cn in frame:

            theta,phi = spherical_coordinates(numpy.array(cn)) # Values used for ORIENTATION 
            thetas.append(theta) #append this data point to lists
            phis.append(phi)
 
    # Now apply symmetry operations for looking at absolute orientation...
    # Exploits full Oh symmetry with two lines - see workings in Jarv's notebooks
    # (rear page of Orange book, 16-4-14)
            cn=abs(cn)
            cn=sorted(abs(cn),reverse=True) #Exploit Oh symmetry - see workings in Jarv's notebooks
            
            thetaOh,phiOh=spherical_coordinates(numpy.array(cn))
            thetasOh.append(thetaOh)
            phisOh.append(phiOh)

    # 2D density plot of the theta/phi information
    fig=plt.figure()
    ax=fig.add_subplot(111)

# Use Cubehelix colourmap, as it is wicked. 
# Read why: http://www.ifweassume.com/2013/05/cubehelix-or-how-i-learned-to-love.html
    plt.hexbin(phis,thetas,gridsize=36,marginals=False,cmap=plt.cm.cubehelix_r) #PuRd) #cmap=plt.cm.jet)
    plt.colorbar()
    pi=numpy.pi

# Full spherical coordinate axes
    plt.xticks( [-pi,-pi/2,0,pi/2,pi],
                [r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],
                fontsize=14)
    plt.yticks( [-pi/2,0,pi/2],
                [r'$-\pi/2$',r'$0$',r'$\pi/2$'],
                fontsize=14)

    if (DisplayFigures):
        plt.show()
    fig.savefig("%s-orientation_density_nosymm.png"%fileprefix,bbox_inches='tight', pad_inches=0,dpi=200)
    fig.savefig("%s-orientation_density_nosymm.pdf"%fileprefix,bbox_inches='tight', pad_inches=0)
    fig.savefig("%s-orientation_density_nosymm.eps"%fileprefix,bbox_inches='tight', pad_inches=0.2)
    
    # 2D density plot of the theta/phi information - MOLLWEIDE projection
    fig=plt.figure()
    ax=fig.add_subplot(111,projection = 'mollweide')

    plt.hexbin(phis,thetas,gridsize=36,marginals=False,cmap=plt.cm.cubehelix_r) #PuRd) #cmap=plt.cm.jet)

    plt.colorbar()
    pi=numpy.pi

# Full spherical coordinate axes
#    plt.xticks( [-pi,-pi/2,0,pi/2,pi],
#                [r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],
#                fontsize=14)
#    plt.yticks( [0,pi/2,pi],
#                [r'$0$',r'$\pi/2$',r'$\pi$'],
#                fontsize=14)
    plt.xticks([],[])

    if (DisplayFigures):
        plt.show()
    fig.savefig("%s-orientation_density_nosymm_mollweide.png"%fileprefix,bbox_inches='tight', pad_inches=0,dpi=200)
    fig.savefig("%s-orientation_density_nosymm_mollweide.pdf"%fileprefix,bbox_inches='tight', pad_inches=0)
    fig.savefig("%s-orientation_density_nosymm_mollweide.eps"%fileprefix,bbox_inches='tight', pad_inches=0.2)
    
    fig=plt.figure()
    ax=fig.add_subplot(111)

    plt.hexbin(phisOh,thetasOh,gridsize=36,marginals=False,cmap=plt.cm.cubehelix_r) #PuRd) #cmap=plt.cm.jet)
    plt.colorbar()
    pi=numpy.pi
   
    # Full symm axes
    plt.xticks( [0.01,pi/4], 
                [r'$0$',r'$\pi/4$'],
                fontsize=14)
    # Now _THAT_ is a magic number. 
    #  (OK, actually it's a Trig identity, just a bit of a complex one)
    plt.yticks( [-0.6154797,-0.01],
                [r'$-0.62$',r'$0$'],
                fontsize=14)
    
    if (DisplayFigures):
        plt.show()
    fig.savefig("%s-orientation_density_Oh_symm.png"%fileprefix,bbox_inches='tight', pad_inches=0,dpi=200)
    fig.savefig("%s-orientation_density_Oh_symm.pdf"%fileprefix,bbox_inches='tight', pad_inches=0)
    fig.savefig("%s-orientation_density_Oh_symm.eps"%fileprefix,bbox_inches='tight', pad_inches=0.2)
 


print("Calculating orientation density...")
sys.stdout.flush()
orientation_density()
print("Calculating correlation...")
sys.stdout.flush()
correlation()

