

// 「シミュレーションの本体」 //
 #include "PI.h"

// シミュレーションを行う関数
 void Simulation(double *r)
{
   int step,judge1,judge2;
   
   printf("\n　シミュレーションを実行中…\n\n");

   /* 熱化後、フリップの最大幅が適正かチェックする */
   while(1){

     /* 経路を初期化 */
     Initialize(r);

     /* 熱化を実施 */
     for(step=1;step<=HEAT_MAX;step++){
         Flipfunc(r,step,1);
     }

     /* 各フリップが5割になっているかをチェック */
     judge1=Localflip(r,2); judge2=Pathflip(r,2);
     if( judge1==0 && judge2==0 ) break;
s
  }

   /* 熱力学的量の期待値を計算 */
   Sampling(r);

}

// フリップ判定を行う関数
 void Flipfunc(double *r,int step,int num)
{
   /* 各スレッドにある分子をフリップさせるか判定する */
   Localflip(r,num);
        
   /* 経路全体をフリップさせるか判定する */
   if(step%10==0){
       Pathflip(r,num);
   }
}

// 総和を計算する関数（Kahanのアルゴリズム）
 void Kahan(double *sum,double *error,double value)
{
   double t,y;

   y=value-(*error);         /* 誤差の反映 */
   t=(*sum)+y;               /* 試しに加算 */
   (*error)=( t-(*sum) )-y;  /* 誤差の算出 */
   (*sum)=t;                 /* 加算した値を代入 */
   
} 
 

