for i
do
 python mdanalysis.py "${i}" 2>&1 | tee -a ${i}.log 
done
