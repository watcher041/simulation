

// 「メトロポリス=ヘイスティング法によるフリップ判定」 //
 #include "PI.h"

// スレッドごとにフリップ判定を行う関数
 int Localflip(double *r,int num)
{
   int n,under;
   double old,xi,Sn,So,dS,w,p,judge;
   static long long int count=0,all=0;
   static double a=AL;

   /* フリップする割合を計算するときに必要 */
   while(num==2){

       judge=(double)count/(double)all;
       count=0; all=0;

       if( 0.46 > judge ){
           printf("　al=%f judge=%f\n",a,judge);
           a*=0.99;
       }else if( 0.54 < judge ){
           printf("　al=%f judge=%f\n",a,judge);
           a*=1.01;
       }else{
           printf("　alの値 a=%f\n",a);
           return 0;
       }

     return 1;
   }

   /* 各スレッドごとにフリップ判定を行う */
   for(n=0;n<3*M;n++){

     /* ひとつ下のスレッドを指す番号 */
     under=(n+3*M-1)%(3*M);

     /* フリップする前の相対確率を計算 */
     So=ST(r,under)+ST(r,n)+SU(r,n)+SF(r,n);

     /* 各成分ごとにフリップさせる */
     xi=1.0-2.0*Random(); old=r[n];
     r[n]=old+a*xi;

     /* フリップさせた後の相対確率を計算 */
     Sn=ST(r,under)+ST(r,n)+SU(r,n)+SF(r,n) ;

     /* 経路を変更するか判定する。 */
     dS=Sn-So;
     if(dS<0.0){
         w=1.0;
     }else{
         w=exp(-tau*dS);
     }
     p = Random();
     if(p>w){
         r[n]=old;
     }else{
         count+=num;
     }
     all+=num;
   }

   return 0;
}

// 経路全体でフリップを行う関数
 int Pathflip(double *r,int num)
{
   int n,k;
   double R[3*M],error;
   double xi,dS,w,p,phase,phi,judge;
   static long long int count=0,all=0;
   static double a=AG;

   while(num==2){
       
       judge=(double)count/(double)all;
       count=0; all=0;

       if( 0.46 > judge ){
           printf("　ag=%f judge=%f\n",a,judge);
           a*=0.99;
       }else if( 0.54 < judge ){
           printf("　ag=%f judge=%f\n",a,judge);
           a*=1.01;
       }else{
           printf("　agの値 a=%f\n",a);
           return 0;
       }

     return 1;
   }
   
   k=(int)round(3.0*Random()); phi=Random();     
   xi=1.0-2.0*Random();
   for(n=0;n<3*M;n++){
       phase=2.0*M_PI*(double)k*(( (double)n/(double)(3*M) )+phi);
       R[n]=r[n]+a*xi*cos(phase);
   }

   dS=0.0; error=0.0;
   for(n=0;n<3*M;n++){
       Kahan( &dS,&error,-(ST(r,n)+SU(r,n)+SF(r,n)) );
       Kahan( &dS,&error,ST(R,n)+SU(R,n)+SF(R,n) );
   }

   if(dS<0.0){
       w=1.0;
   }else{
       w=exp(-tau*dS);
   }

   p = Random();

   if(p<=w){
     for(n=0*M;n<3*M;n++){
         r[n]=R[n];
     }
     count+=num;
   }
   all+=num;
        
   return 0;
  
}

