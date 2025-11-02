
#include"00_PI.h"

// 経路全体でフリップを行う関数
void globalflip(double *r)
{
   int n,k;
   double r_copy[3*M];
   double xi,ds,w,phase,phi;
   double ag = sqrt(TAU);
   
   k=(int)round(3.0*Random()); phi=Random();     
   xi=1.0-2.0*Random();
   for(n=0;n<3*M;n++){
       phase=2.0*M_PI*(double)k*(( (double)n/(double)(3*M) )+phi);
       r_copy[n]=r[n]+ag*xi*cos(phase);
   }

   ds=0.0; 
   for(n=0;n<3*M;n++){
        ds+=Tn(r_copy,n)-Tn(r,n)+
            Un(r_copy,n)-Un(r,n)+
            Fn(r_copy,n)-Fn(r,n);
   }

   if(ds<0.0){
       w=1.0;
   }else{
       w=exp(-TAU*ds);
   }
   if(Random()<=w){
     for(n=0;n<3*M;n++){
         r[n]=r_copy[n];
     }
   }
}

