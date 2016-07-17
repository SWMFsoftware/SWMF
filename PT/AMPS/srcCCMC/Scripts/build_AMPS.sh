module load mpi/openmpi-x86_64

cd AMPS

./Config.pl -install -compiler=gfortran
./Config.pl -application=ccmc-individual_trajectories--kameleon 
./Config.pl -boost-path=/home2/ccmc/kameleon/include/boost 
./Config.pl -kameleon-path=/home2/ccmc/kameleon/include 
./Config.pl -cplr-data-path=.
./Config.pl -model-data-path=.
./Config.pl -link-option=-L,/home2/ccmc/kameleon/lib,-L,/home2/ccmc/kameleon/lib/ccmc,-lccmc,-lcdf,-lhdf5,-lboost_thread,-lpython2.7,-lstdc++,-lmpi_cxx,-lboost_wave,-lboost_python,-lhdf5_hl_cpp,-lboost_program_options,-lboost_filesystem,-lhdf5_cpp,-lhdf5_hl,-lboost_date_time,-lboost_chrono,-lboost_system,/home2/ccmc/kameleon/lib/ccmc/pyKameleon.so

make -j

cp amps ..
cd ..
