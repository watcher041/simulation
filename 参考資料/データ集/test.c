#include <stdio.h>

 int main(void)
{
   int M;
   double T;

  for(M=1;M<=100;M++){
      T=1.0/(0.010*(double)M);
      printf("M=%d T=%f\n",M,T);
  }

   return 0;
}

