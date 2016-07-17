csh build_AMPS.sh


cd AMPS
cp ../amps .
cp ../Makefile.test .

export LD_LIBRARY_PATH=/home2/ccmc/kameleon/lib/:/home2/ccmc/kameleon/lib/ccmc/
make -f Makefile.test test_CCMC-Individual_Trajectories--KAMELEON



