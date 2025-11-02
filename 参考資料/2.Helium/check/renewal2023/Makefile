

# 分割コンパイルを行うスクリプト #

COMPILER = g++
OBJECT = $(shell echo *.c  | sed 's/\.c/\.o/g'  )
THREAD=$(shell grep -o 'define M [0-9]*' *.h | grep -o '[0-9]*' )
OPTION = -O3 -Wall -std=c++11

test : ${OBJECT}
	@${COMPILER} -o test${THREAD} ${OPTION} ${OBJECT} -lm
	@make clean
	@reset
	@echo "make test${THREAD}"
.c.o:
	@${COMPILER} -c ${OPTION} $<
clean:
	@rm -f  *.o core 
clear:
	@rm *.c~ *.h~ *.sh~ Makefile~ test* *.dat

