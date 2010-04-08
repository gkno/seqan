CC = g++
CFLAGS = -O0 -g3 -Wall -c -fmessage-length=0
LIB_SEQAN = -IF:\Development\workspace_cpp_test\seqan\projects\library
DIR_BIN = binaries/

####################################################
#Target CPPKurs
journal: main.o
	${CC} main.o -o ${DIR_BIN}journal.exe 
	@echo 'ready'
	
main.o: main.cpp
	${CC} ${CFLAGS} $(LIB_SEQAN) main.cpp
	
####################################################
#Target clean

clean:
	rm -rf *o ${DIR_BIN}*.exe *.d