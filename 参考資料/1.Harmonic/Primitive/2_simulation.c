


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

  }

  /* シミュレーションに関わる定数を出力 */
   printf("\n");
   printf("（シミュレーションに関わる定数）\n");
   printf("　温度 T=%f\n",TEMP);
   printf("　分割数 τ=%f\n",M);

   /* 理論値を出力 */
   double beta=1.0/TEMP;
   double x=0.5*beta;
   double E=0.5/tanh(x);
   double C=(beta*beta)/( 4.0*sinh(x)*sinh(x) );
   printf("　理論値 E=%e C=%e\n\n",E,C);

   /* 熱力学的量の期待値を計算 */
   Sampling(r);

   /* 見やすくするために、改行を入れておく。 */
   printf("\n");
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
 

