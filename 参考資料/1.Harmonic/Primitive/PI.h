


//  ＜1体DIM次元調和振動子におけるシミュレーション＞  //

// 座標ごとにフリップを行うと良い。

// ヘッダーファイルの指定
 #include <stdio.h>
 #include <math.h>
 
// 観測に関わる定数 
 #define M 1
 #define TEMP 0.2            /* 温度ごとにつけた番号 */
 #define MCS_MAX 1.0e+08     /* サンプル数 */
 #define HEAT_MAX 1.0e+04    /* 熱化の実行回数（実際は指定した値の二倍行う）	 */
 #define AL 5.00e-01         /* 局所更新の初期最大幅 */
 #define AG 5.00e-01         /* 対局更新の初期最大幅 */
 #define N_SAMP 10           /* サンプリングの試行回数 */  
 #define R0 0.0              /* スレッドごとの経路の初期配置（ポテンシャルに応じて変更） */ 
 
// 期待値を計算する必要がある量のマクロ関数
 #define ST(X1,X2) ( 0.5*Tn(X1,X2)/(tau*tau) )
 #define SU(X1,X2) ( Un(X1,X2) )
 
// グローバル変数の定義
 extern double tau;

// プロトタイプ宣言
 void Initialize(double *);
 void Simulation(double *);
 void Sampling(double *);
 void Flipfunc(double *,int,int);
 void Kahan(double *,double *,double);
 int Localflip(double *,int);
 int Pathflip(double *,int);
 double Tn(double *,int);
 double Un(double *,int);
 double U(double);
 double Random();
 
 
