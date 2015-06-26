grep Pb aggregate_222.xyz > Pb.dat
grep I aggregate_222.xyz > I.dat
grep ^C aggregate_222.xyz > C.dat

cat Pb.dat | ./Pb_I_distribution_PBCS.awk > Pb_symm.dat
cat I.dat | ./Pb_I_distribution_PBCS.awk > I_symm.dat
cat C.dat | ./Pb_I_distribution_PBCS.awk > C_symm.dat

