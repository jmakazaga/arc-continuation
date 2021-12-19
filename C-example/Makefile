
# MHDIR must be the directory where the file arc-continuation.c is placed
MHDIR = ..


INCDIR = -I./ -I${MHDIR}

CC = /usr/bin/gcc
CXX = /usr/bin/c++

BINARIES = mainArcContinuation
# in the example we do not have users functions in separate file, so we don't need user_fcns.o
#OBJECTS = mainArcContinuation.o ${MHDIR}/arc-continuation.o user_fcns.o
OBJECTS = mainArcContinuation.o ${MHDIR}/arc-continuation.o 

all : $(BINARIES)

mainArcContinuation: ${OBJECTS}
	${CC} ${EXTRA_CFLAGS} -I${INCDIR}   ${OBJECTS}  -lm -llapack -o $@

mainArcContinuation.o: mainArcContinuation.c 
	${CC} -c -I${INCDIR}   mainArcContinuation.c -o $@

${MHDIR}/arc-continuation.o: ${MHDIR}/arc-continuation.c 
	${CC} -c -I${INCDIR}   ${MHDIR}/arc-continuation.c -o $@

#if users functions (rest and jacobian were in separate file uncomment next:
#user_fcns.o : user_fcns.c
#	${CC} ${EXTRA_CFLAGS} -I${INCDIR}  -c  user_fcns.c  -o $@
