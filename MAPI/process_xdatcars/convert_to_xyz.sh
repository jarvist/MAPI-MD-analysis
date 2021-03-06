scp ultra:/home/pr1u1109/pr1u1109/pr1ujarv/work/MAPI-MD/testmd2-nonselective/XDATCAR ./testmd2-nonselective_XDATCAR
scp ultra:/home/pr1u1109/pr1u1109/pr1ujarv/work/MAPI-MD/testmd3/XDATCAR ./testmd3_XDATCAR
scp ultra:/home/pr1u1109/pr1u1109/pr1ujarv/work/MAPI-MD/testmd4/XDATCAR ./testmd4_XDATCAR
scp ultra:/home/pr1u1109/pr1u1109/pr1ujarv/work/MAPI-MD/testmd5/XDATCAR ./testmd5_XDATCAR

source ~/Virtualenvs/python-ase-3.8.1.3440/bin/activate

python ./Xdat2Xyz.py -f ./testmd2-nonselective_XDATCAR -o ./testmd2-nonselective.xyz
python ./Xdat2Xyz.py -f ./testmd3_XDATCAR -o ./testmd3.xyz
python ./Xdat2Xyz.py -f ./testmd4_XDATCAR -o ./testmd4.xyz
python ./Xdat2Xyz.py -f ./testmd4_XDATCAR -o ./testmd5.xyz

cat ./testmd2-nonselective.xyz ./testmd3.xyz ./testmd4.xyz ./testmd5.xyz > aggregate_222.xyz

