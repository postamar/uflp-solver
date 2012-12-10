# (C) 2010-2011, by Marius Posta. 
# Check LICENSE.txt for the legal blah-blah.
#
#
#
# compilation flags
#
FLAGS = -m64 -O3 -DNDEBUG -ffast-math #-march=native 
COMPILER = g++

all: Bundle.o MinQuad.o BMinQuad.o QPBundle.o bundle_solver.o
	ar rcs libbtt.a Bundle.o MinQuad.o BMinQuad.o QPBundle.o bundle_solver.o

#
# Antonio Frangioni's [B]TT
#

Bundle.o: Bundle.C 
	$(COMPILER) -c $(FLAGS) $*.C -o $@ 

MinQuad.o: MinQuad.C
	$(COMPILER) -c $(FLAGS) $*.C -o $@ 

BMinQuad.o: BMinQuad.C 
	$(COMPILER) -c $(FLAGS) $*.C -o $@ 

QPBundle.o: QPBundle.C 
	$(COMPILER) -c $(FLAGS) $*.C -o $@


#
# Paul-Virak Khuong's C99 interface to the object-oriented mess that is [B]TT
#
bundle_solver.o: bundle_solver.C 
	$(COMPILER) -c $(FLAGS) $*.C -o $@ 


#
# clean rule: deletes temporary files.
#
clean :
	/bin/rm -rf *.o *~ *.class *.a
	/bin/rm -rf *.mps *.ord *.sos *.lp *.sav *.net *.msg *.log


