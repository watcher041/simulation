

// ＜ 並列計算を行うためのプログラム１ ＞ //

/* 各コアごとに「Hello world!」を出力するもの */

#include <stdio.h>

int main(void)
{

  #pragma omp parallel
  {
      printf("Hello world!\n");
  }
  
  return 0;

}
