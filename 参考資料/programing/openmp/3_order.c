
// コアの順番通りに実行

#include<stdio.h>
#include<omp.h>

int main(void)
{

  int i,a[100];

  #pragma omp parallel for ordered
  for(i=0;i<100;i++){
      a[i]=0;
      #pragma omp ordered
      printf("%d a[%d]=%d\n",omp_get_thread_num(),i,a[i]);
  }

}
