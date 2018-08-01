E_vs_it:	E_vs_it.o nrutil.o solvde.o difeq.o bksub.o red.o pinvs.o qromb.o trapzd.o polint.o
	gcc -o E_vs_it E_vs_it.o nrutil.o solvde.o difeq.o bksub.o red.o pinvs.o qromb.o trapzd.o polint.o -lm

E_vs_it.o:	E_vs_it.c nrutil.h headerfile.h
	gcc -c E_vs_it.c

nrutil.o:	nrutil.c nrutil.h
	gcc -c nrutil.c

solvde.o:	solvde.c headerfile.h
	gcc -c solvde.c

difeq.o:	difeq.c headerfile.h
	gcc -c difeq.c

bksub.o:	bksub.c headerfile.h
	gcc -c bksub.c

red.o:	red.c headerfile.h
	gcc -c red.c

pinvs.o:	pinvs.c headerfile.h
	gcc -c pinvs.c

qromb.o:	qromb.c	headerfile.h
	gcc -c qromb.c

trapzd.o:	trapzd.c	headerfile.h
	gcc -c trapzd.c

polint.o:	polint.c	headerfile.h
	gcc -c polint.c
