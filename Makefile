# (C) 2011, Marius Posta (postamar@iro.umontreal.ca)
# Check LICENSE.txt for the legal blah-blah.

CFLAGS = -m64 -W -Wall -O3 -ffast-math
LFLAGS = -lm -lblas -L./BTT -lbtt 
COMPILER = clang 


solver: dual.o node.o node_aux.o tighten.o primal.o instance.o main.o RngStream.o
	cd BTT; make all
	g++ $(CFLAGS) $^ -o solver $(LFLAGS)

main.o: main.c instance.h RngStream.h primal.h dual.h node.h
	$(COMPILER) -c $(CFLAGS) $*.c -o $@

instance.o: instance.c instance.h RngStream.h
	$(COMPILER) -c $(CFLAGS) $*.c -o $@

RngStream.o: RngStream.c RngStream.h 
	$(COMPILER) -c $(CFLAGS) $*.c -o $@

primal.o: primal.c primal.h instance.h RngStream.h
	$(COMPILER) -c $(CFLAGS) $*.c -o $@

dual.o: dual.c dual.h node.h instance.h RngStream.h
	$(COMPILER) -c $(CFLAGS) $*.c -o $@ 

node.o: node.c node.h instance.h RngStream.h
	$(COMPILER) -c $(CFLAGS) $*.c -o $@

node_aux.o: node_aux.c node.h instance.h RngStream.h
	$(COMPILER) -c $(CFLAGS) $*.c -o $@

tighten.o: tighten.c dual.h node.h instance.h RngStream.h 
	$(COMPILER) -c $(CFLAGS) $*.c -o $@ 


clean :
	/bin/rm -rf *.o *~ *.class
	/bin/rm -rf *.mps *.ord *.sos *.lp *.sav *.net *.msg *.log




