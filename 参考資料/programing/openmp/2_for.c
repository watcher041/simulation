
#include<stdio.h>

int main(void)
{

  int i,a[100];

  #pragma omp parallel
  {
    #pragma omp for
    for(i=0;i<100;i++){
        a[i]=0; 
        printf("%d\n",a[i]);
    }
  }

  #pragma omp parallel for
  for(i=0;i<100;i++){
      a[i]=1;
      printf("%d\n",a[i]);
  }

}
