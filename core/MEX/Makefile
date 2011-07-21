MLAB = /usr/share/matlab

MXSFX = mexlx
MEX = ${MLAB}/bin/mex

MFLAGS = -V5
INC = -I ${MLAB}/extern/include
LIB = -L ${MLAB}/extern/lib

CC = gcc -Wall -pedantic
C_OPTIMIZE_SWITCH = -O2    ## For GCC
CFLAGS = ${C_OPTIMIZE_SWITCH} ${INC} ${LIB}

all: pointOp.${MXSFX} demod2d_mod.${MXSFX} demod2d_c.${MXSFX}

clean:
	/bin/rm *.o

pointOp.${MXSFX}: pointOp.o
	${MEX} ${MFLAGS} pointOp.o

demod2d_c.${MXSFX}: demod2d_c.o
	${MEX} ${MFLAGS} demod2d_c.o

demod2d_mod.${MXSFX}: demod2d_mod.o
	${MEX} ${MFLAGS} demod2d_mod.o

%.o : %.c
	${CC} -c ${CFLAGS} $<		

