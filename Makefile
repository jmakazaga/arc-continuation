# This makefile if for mathematica 

#VERSION=10.0
#SYS = Linux-x86-64
#MLLIB =ML64i3

# MHDIR must be the directory where the files arc-continuation.c and terminalArcContinuation.c are placed
# MHDIR =/home/joseba/OwnCloud/ik/arc-continuation
#MHDIR =/home/eli/ik/homotopia_double

ifeq (${SYS}, Linux-x86-64)
	MLLIB = ML64i3
else
	MLLIB = ML32i3
endif

MLINKDIR = /usr/local/Wolfram/Mathematica/${VERSION}/SystemFiles/Links/MathLink/DeveloperKit/${SYS}/CompilerAdditions
CADDSDIR = ${MLINKDIR}
EXTRA_CFLAGS= -O3

INCDIR = ${CADDSDIR} -I./ -I${MHDIR}
LIBDIR = ${CADDSDIR}

# If you have the command mprep simbolicaly linkeked in /usr/bin uncomment next line
# MPREP = mprep
# if not in /usr/bin then uncomment next one
MPREP = ${CADDSDIR}/mprep
RM = rm

CC = /usr/bin/gcc
CXX = /usr/bin/c++

BINARIES = mathArcContinuation

all : $(BINARIES)


mathArcContinuation : ${MHDIR}/mathArcContinuationtm.o ${MHDIR}/mathArcContinuation.o ${MHDIR}/arc-continuation.o user_fcns.o
	${CC} ${EXTRA_CFLAGS} -I${INCDIR} -L${MLINKDIR}  ${MHDIR}/arc-continuation.o ${MHDIR}/mathArcContinuationtm.o user_fcns.o  -l${MLLIB} -lm -luuid -lpthread -lrt -lstdc++ -lm -llapack -o $@


terminalArcContinuation: ${MHDIR}/terminalArcContinuation.o ${MHDIR}/arc-continuation.o user_fcns.o
	${CC} ${EXTRA_CFLAGS} -I${INCDIR}   ${MHDIR}/terminalArcContinuation.o  user_fcns.o -lm -llapack -o $@


${MHDIR}/terminalArcContinuation.o: ${MHDIR}/terminalArcContinuation.c 
	${CC} -c -I${INCDIR}   ${MHDIR}/terminalArcContinuation.c -o $@


${MHDIR}/mathArcContinuation.o: ${MHDIR}/mathArcContinuation.c 
	${CC} -c -I${INCDIR}   ${MHDIR}/terminalArcContinuation.c -o $@


${MHDIR}/arc-continuation.o : ${MHDIR}/arc-continuation.c
	${CC} ${EXTRA_CFLAGS} -I${INCDIR}  -c ${MHDIR}/arc-continuation.c -o $@ 

${MHDIR}/mathArcContinuationtm.o : ${MHDIR}/mathArcContinuationtm.c
	${CC} ${EXTRA_CFLAGS} -I${INCDIR}  -c ${MHDIR}/mathArcContinuationtm.c -o $@

user_fcns.o : user_fcns.c
	${CC} ${EXTRA_CFLAGS} -I${INCDIR}  -c  user_fcns.c  -o $@

${MHDIR}/mathArcContinuationtm.c : ${MHDIR}/mathArcContinuation.tm
	${MPREP} ${MHDIR}/mathArcContinuation.tm -o ${MHDIR}/mathArcContinuationtm.c


clean :
	@ ${RM} -rf *.o ${MHDIR}/*.o ${MHDIR}/*tm.c $(BINARIES)
