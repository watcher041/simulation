
 #include "split.h"

 int main(void)
{
   #pragma omp parallel num_threads(3)
   output();

   int i,sum,dammy;

   sum=0;
   #pragma omp parallel for reduction(+:sum) private(dammy) num_threads(3)
   for(i=1;i<=100;i++){
       dammy=i;
       sum+=dammy;
   }

   printf("sum=%d\n",sum);

   return 0; 
}
