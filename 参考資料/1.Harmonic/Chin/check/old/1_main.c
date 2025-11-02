

// 「プログラムの初期化」 //
 #include "PI.h"
 #include <stdlib.h>
 #include <time.h>

 double tau;

// メイン関数（初期設定を行う関数）
 int main(void)
{
   /* 逆温度βとスレッド数Mの計算 */
   double beta=1.0/TEMP; tau=beta/(double)M;

   /* 経路を形を記録する配列を定義 */
   double r[3*M];

   /* 経路の初期化 */
   Initialize(r);

   /* 乱数の初期化 */
   srand48((unsigned) time(NULL));

     /* シミュレーションに関わる定数を出力 */
   printf("\n");
   printf("（シミュレーションに関わる定数）\n");
   printf("　温度 T=%f\n",TEMP);
   printf("　分割幅 τ=%f\n",tau);

   /* 理論値を出力 */
   double x=0.5*beta;
   double E=0.5/tanh(x);
   double C=(beta*beta)/( 4.0*sinh(x)*sinh(x) );
   printf("　理論値 E=%e C=%e\n",E,C);

   /* シミュレーションを実行 */
   Simulation(r);

   /* シミュレーションを終了 */
   return 0;
}

// 経路の初期化を行う関数
 void Initialize(double *r)
{
   int n;

   /* 位置の初期化 */
   for(n=0;n<3*M;n++){
       r[n]=R0;
   }
}

// 一様乱数を作成する関数
 double Random()
{
   double answer;

   /* 低温でBoltzmann因子と比較できるようにするためdrand48を用いる */
   answer=drand48();

   return answer;
}

