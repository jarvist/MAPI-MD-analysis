#!/usr/bin/awk -f
BEGIN {
    unitcell=12.57398/2.0 #average vector of pseudo cubic unit cell
    # divided by two as we're a 2x2x2 supercell (MAGIC NUMBER)
}
{
    pos[1]=$2
    pos[2]=$3
    pos[3]=$4

    for (i=1;i<=3;i++)
    {
        while (pos[i] >= unitcell) # while any lattice vector > unitcell, decrement by unitcell
            pos[i]=pos[i]-unitcell
    }

#    print $1,pos[1],pos[2],pos[3]   # XYZ output angstrom
    print $1,pos[1]/unitcell,pos[2]/unitcell,pos[3]/unitcell # Fractional coordinates
}
