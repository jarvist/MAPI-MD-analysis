# mdanalysis_cn_dist.py - Jarvist Moore Frost April 2014
#  Scripts to analyise CN axis orientation from Ab-Initio MD on MAPI

#When using MDAnalysis in published work, please cite doi:10.1002/jcc.21787
import MDAnalysis
import MDAnalysis.analysis

import matplotlib.pyplot as plt

import numpy
import math
import random

from IPython import embed #iPython magic for interactive session...

#universe item; contains all the mdanalysis stuff
u= MDAnalysis.Universe("1stframe.pdb","process_xdatcars/aggregate_222.xyz") #Trajectory.xyz")#MAPI_222_equilibr_traj.xyz")

# Above .XYZ file was generated with Keith Butler's Xdat2Xyz.py (depends on ASE)
# The pdb file was generated with Open Babel

#embed()

GenThetas=False    # Theta / Phis to STDOUT for plotting externally
GenXYZ=False        # .XYZ file to STDOUT of CN axes, for Pymol plotting
ExploitSymm=True  # Exploit full symmetry = 48 fold
Exploit8fold=False # Exploit 8-fold symmetry

thetas=[] # List to collect data for later histogramming
phis=[]

dotcount=[0.,0.,0.]

carbons=u.selectAtoms('name C')
nitrogens=u.selectAtoms('name N')

#Hard wired cubic PBCs - the MDAnalysis readers don't seem to pick these up from PDB / XYZ files (?)
mybox=numpy.array([12.5801704656605438,12.5477821243936827,12.5940169967174249],dtype=numpy.float32)
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

    for unit in [1,0,0],[1,1,0],[1,1,1]:
        unit=unit/numpy.linalg.norm(unit)
#       print "dot:", unit,numpy.dot(unit,cn),theta,phi
        dots.append(numpy.dot(unit,cn))
    dotcount[numpy.argmax(dots)]=dotcount[numpy.argmax(dots)]+1
#    print dotcount
    print cn[0],cn[1],cn[2],numpy.argmax(dots) #x,y,z, which type (for colour)
    print "FaceTheta: ", math.acos(min(dots[0],1.0)) #fix to domain errors as dots[] creeping above 1.0
    print "EdgeTheta: ", math.acos(min(dots[1],1.0))
    print "DiagTheta: ", math.acos(min(dots[2],1.0))

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

                theta,phi = spherical_coordinates(numpy.array(cn))

                thetas.append(theta) #append this data point to lists
                phis.append(phi)

#                partition_alignment(cn)

                #OK, now we fold along theta, phi, to account for symmetry (TODO: Check!)
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

print dotcount

# 2D density plot of the theta/phi information
fig=plt.figure()
ax=fig.add_subplot(111)


plt.hexbin(phis,thetas,gridsize=36,marginals=False,cmap=plt.cm.jet)
plt.colorbar()
pi=numpy.pi

# Full spherical coordinate axes
#plt.xticks( [-pi,-pi/2,0,pi/2,pi],
#            [r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],
#            fontsize=14)
#plt.yticks( [0,pi/2,pi],
#            [r'$0$',r'$\pi/2$',r'$\pi$'],
#            fontsize=14)

# Full symm axes
plt.xticks( [0,pi/4],
            [r'$0$',r'$\pi/4$'],
            fontsize=14)
plt.yticks( [0.9553166181245,pi/2],
            [r'$arcos(\frac{1}{\sqrt{3}})$',r'$\pi/2$'],
            fontsize=14)



plt.show()

fig.savefig("mdanalysis_cn_dist.png",bbox_inches='tight', pad_inches=0)
fig.savefig("mdanalysis_cn_dist.pdf",bbox_inches='tight', pad_inches=0)


end

phi_bins=[]
for i in range(24):
    phi_bins.append( 0.3*math.pi+math.asin(i/72.))
#theta_bins=numpy.arange(0.,math.pi/4,math.pi/72)

theta_bins  =numpy.arange(0.,math.pi/4,math.pi/72)

print theta_bins, phi_bins

H,xedges,yedges = numpy.histogram2d(thetas,phis,bins=36)
H2, _, _ = numpy.histogram2d(thetas,phis,bins=[phi_bins,theta_bins])

H.shape, xedges.shape, yedges.shape
extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]

#Contours - via http://micropore.wordpress.com/2011/10/01/2d-density-plot-or-2d-histogram/
# - Data are too noisy for this to be useful. Also, they're upside down?! Weird!
#fig.subplots_adjust(bottom=0.15,left=0.15)
#levels = (5.0e1, 4.0e1, 3.0e1, 2.0e1)
#cset = plt.contour(H, levels, origin='lower',colors=['black','green','blue','red'],linewidths=(1.9, 1.6, 1.5, 1.4),extent=extent)
#plt.clabel(cset, inline=1, fontsize=10, fmt='%1.0i')
#for c in cset.collections:
#    c.set_linestyle('solid')

plt.imshow(H,extent=extent,interpolation='nearest')
plt.colorbar()
plt.show()

plt.imshow(H2,extent=extent,interpolation='nearest')
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
