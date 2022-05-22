CFLAGS=-Wall -Wextra -O3
LDFLAGS= -lgsl -lgslcblas -lconfig++
CC = g++

sim: simulate.cc 
	${CC} ${CFLAGS} -o simulate simulate.cc ${LDFLAGS} 
