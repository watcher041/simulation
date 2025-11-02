

OBJECT = $(shell echo *.c  | sed 's/\.c/\.o/g'  )
TEMP=$(shell grep -o 'M [0-9]*' 00_PI.h | grep -o '[0-9]*' )
OPTION = -O3 -Wall -std=c++11

test : ${OBJECT}
	@g++ -o test${TEMP} ${OPTION} ${OBJECT} -lm
	@make clean
	@reset
	@echo "make test${TEMP}"
.c.o:
	@g++ -c ${OPTION} $<
clean:
	@rm -f  *.o core 
clear:
	@rm *.c~ *.h~ *.sh~ Makefile~ test* *.dat

