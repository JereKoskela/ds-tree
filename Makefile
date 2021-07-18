CFLAGS=-Wall -Wextra -O3
LDFLAGS= -lgsl -lgslcblas -lconfig++
CC = g++


main: abc.cc 
	${CC} ${CFLAGS} -o abc abc.cc ${LDFLAGS} 

sim: simulate.cc 
	${CC} ${CFLAGS} -o simulate simulate.cc ${LDFLAGS} 
