# mdanalysis_cn_dist.py - Jarvist Moore Frost April 2014
#  Scripts to analyise CN axis orientation from Ab-Initio MD on MAPI

#When using MDAnalysis in published work, please cite doi:10.1002/jcc.21787
import MDAnalysis
import MDAnalysis.analysis

import numpy
import math
from IPython import embed #iPython magic for interactive session...

#universe item; contains all the mdanalysis stuff
u= MDAnalysis.Universe("1stframe.pdb","MAPI_222_equilibr_traj.xyz")

#embed()

carbons=u.selectAtoms('name C')
nitrogens=u.selectAtoms('name N')

mybox=numpy.array([12.5801704656605438,12.5477821243936827,12.5940169967174249],dtype=numpy.float32)
imybox=1/mybox
print mybox, imybox

for ts in u.trajectory:
    r=MDAnalysis.analysis.distances.distance_array(carbons.coordinates(),nitrogens.coordinates(),mybox)
    # OK, this is now the distance array of all Ns vs. Cs, considering PBCs
    print r
    for carbon, nitrogenlist in enumerate(r):
        for nitrogen,distance in enumerate(nitrogenlist):
            if distance<1.6: 
                print carbon, nitrogen, distance
                d=carbons[carbon].pos - nitrogens[nitrogen].pos
                # Filthily hack in some orthorhombic (only) PBCs as I can't see how MDanalysis does this properly
                for dim in 0,1,2:
                    s=imybox[dim] * d[dim] # Minimum image convetion, algo from MDAnalysis calc_distances.h
                    d[dim]=mybox[dim]*(s-round(s)) 
                print d
                

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
