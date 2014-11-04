DEFAULT_TARGET : picsolver 
CC=icpc    
#CC=/usr/local/mpi/bin/mpicxx
#CC=mpicxx
#CC=openmpicxx
#CC=/opt/intel/bin/icpc 

#CC=/opt/local/bin//openmpic++
#CC=/left/Compiler/Intell/11.1/bin/mpicxx
 


#SOURCES=/Users/vtenishe/Debugger/eclipse-workspace/SOURCES
#CWD=/Users/vtenishe/Debugger/eclipse-workspace/mars
CWD=/nobackupp8/ylee5/SIMULATION/HotOxygen/IntensityCal_MGITM_test/ALICE/4BG_MGITM_4AltPts-03hr45min-AEQUMAX-UT14
#CWD=/Users/yunilee/1D_MAVEN
#SPICE=/Users/vtenishe/SPICE/Toolkit/cspice

Lib=  -lm -lmpi 
SHELL=/bin/bash

clean:
	cd ${CWD}/srcTemp/general; rm -f *.o *.a
	cd ${CWD}/srcTemp/pic; rm -f *.o *.a
	cd ${CWD}/srcTemp/meshAMR; rm -f *.o *.a
	rm -f *.o *.a a.out

tar:
	cd ../pic-tower/sources/general; rm -f *.o *.a
	cd ../pic-tower/sources/dsmc; rm -f *.o *.a
	tar -cvf sources.tar sources

#Flags=-g -fno-inline -ffinite-math-only  -ftrapping-math -fsignaling-nans -wd383 -wd981 -wd869 -wd1418 -LANG:std -Wall -DMPI_ON -DParticleSolver=dsmc -DCompilationTarget=0
#Flags=-O3  -fasm-blocks -use-asm  -fprefetch-loop-arrays -funroll-loops -unroll=3  -mmmx  -wd383 -wd981 -wd869 -wd1418 -LANG:std -Wall -DMPI_ON -DParticleSolver=dsmc -DCompilationTarget=0
#Flags=-g   -use-asm   -fprefetch-loop-arrays -funroll-loops -unroll=3    -LANG:std -Wall -DMPI_ON -DParticleSolver=dsmc -DCompilationTarget=0

#Flags=-O3 -fasm-blocks  -wd383 -wd981 -wd869 -wd1418 -LANG:std -Wall -DMPI_ON -DParticleSolver=dsmc -DCompilationTarget=0

Flags=-O3  -wd383 -wd981 -wd869 -wd1418 -LANG:std -Wall -DMPI_ON -DParticleSolver=dsmc -DCompilationTarget=0
IncludeList=-I${CWD}/srcTemp/pic -I${CWD}/srcTemp/main  -I${CWD}/srcTemp/meshAMR -I${CWD}/srcTemp/general  -I${CWD}/srcTemp/species -I${CWD}/srcTemp/models/exosphere 

.SUFFIXES: .cpp .o 
.cpp.o: 
	${CC} -c ${Flags} ${IncludeList} $< 
 
picsolver :

	cd ${CWD}/srcTemp/general; rm -f *.o *.a
	cd ${CWD}/srcTemp/pic; rm -f *.o *.a
	cd ${CWD}/srcTemp/meshAMR; rm -f *.o *.a
	cd ${CWD}/srcTemp/main; rm -f *.o *.a 
	cd ${CWD}/srcTemp/species; rm -f *.o *.a
	rm -f *.o *.a a.out

	cd ${CWD}/srcTemp/general; make -f ${CWD}/makefile general.a
	cd ${CWD}/srcTemp/meshAMR; make -f ${CWD}/makefile mesh.a
	cd ${CWD}/srcTemp/pic; make -f ${CWD}/makefile pic.a
	cd ${CWD}/srcTemp/species; make -f ${CWD}/makefile species.a
	cd ${CWD}/srcTemp/main; make -f ${CWD}/makefile mars.a

	${CC} -o mars srcTemp/main/mars.a srcTemp/species/species.a srcTemp/pic/pic.a srcTemp/general/general.a srcTemp/meshAMR/mesh.a  ${Lib}

PICOBJ = pic.o pic_parallel.o pic_bc.o pic_init_const.o pic_mesh.o pic_mol.o pic_mol_parser.o pic_mover.o pic_parser.o pic_pbuffer.o pic_weight_time.o pic_internal_boundary_sphere.o pic_sample_distribution_function.o pic_sample_flux_distribution.o pic_vel_distribution.o pic_ices.o pic_photolytic_reactions.o pic_volume_particle_injection.o pic_background_atmosphere.o pic_column_integration.o  
pic.a : $(PICOBJ) 
	ar -src pic.a $(PICOBJ) 

MESHOBJ = meshAMRgeneric.o
mesh.a : $(MESHOBJ)
	ar -src mesh.a $(MESHOBJ)

GENERALOBJ = specfunc.o rnd.o 
general.a : $(GENERALOBJ)
	ar -src general.a $(GENERALOBJ)

SPECIESOBJ = Na.o
species.a : $(SPECIESOBJ)
	ar -src species.a $(SPECIESOBJ)

EXOSPHEREOBJ = Exosphere.o Exosphere_OrbitalMotion.o Exosphere_SourceProcesses.o Exosphere_Parallel.o Exosphere_ReferenceGroundBasedObservations.o
exosphere.a : $(EXOSPHEREOBJ)
	ar -src exosphere.a $(EXOSPHEREOBJ)

MARSOBJ = mars.o newMars.o MarsBackgroundAtmosphereFox.o
mars.a : $(MARSOBJ)
	ar -src mars.a $(MARSOBJ)

testmesh:
	make clean
	cd mesh; make -f ../makefile mesh.a
	cd general; make -f ../makefile general.a
	${CC} -c ${Flags} -Igeneral -Imesh -Ipic -I.. -I. 3dMesh.cpp  
	${CC} -o a.out 3dMesh.o mesh/mesh.a general/general.a ${Lib}
