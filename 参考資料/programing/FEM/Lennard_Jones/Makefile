
#＜clapackを用いたプログラムをコンパイルするスクリプト＞#

OPTION = -Wall -O3 -o
LIB = -llapack -lf2c -lblas -lm 

model:
	@gcc ${OPTION} test LJ.c ${LIB}
	
