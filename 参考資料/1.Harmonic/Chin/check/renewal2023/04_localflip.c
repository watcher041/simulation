
#include"00_PI.h"
#include <string.h>

// スレッドごとにフリップ判定を行う関数
void localflip(double *r)
{
   int n,under;
   double xi,ds,w;
   double al = sqrt(TAU);
   double r_copy[3*M];
   
   for(n=0;n<3*M;n++){
      r_copy[n] = r[n];
   }

   /* 各スレッドごとにフリップ判定を行う */
   for(n=0;n<3*M;n++){

     /* 各成分ごとにフリップさせる */
     xi=1.0-2.0*Random();
     r_copy[n]=r[n]+al*xi;

     /* ひとつ下のスレッドを指す番号 */
     under=(n+3*M-1)%(3*M);

     /* フリップさせた後の相対確率を計算 */
     ds=Tn(r_copy,under)-Tn(r,under)+
        Tn(r_copy,n)-Tn(r,n)+
        Un(r_copy,n)-Un(r,n)+
        Fn(r_copy,n)-Fn(r,n);

     /* 経路を変更するか判定する。 */
     if(ds<0.0){
         w=1.0;
     }else{
         w=exp(-TAU*ds);
     }
     if(Random()<=w){
         r[n]=r_copy[n];
     }else{
         r_copy[n]=r[n];
     }
   }
}