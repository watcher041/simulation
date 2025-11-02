

 /* ヘッダーファイルの指定 */
  #include <stdio.h>
  #include <math.h>
  #include <stdlib.h>
  #include <random>

 /* シミュレーションに関わる定数 */
  #define M 50                   /* スレッド数ごとにつけた番号 */
  #define TAU 0.010              /* 分割幅[K^(-1)] */
  #define MCS_MAX 4.0e+05        /* モンテカルロステップ */
  #define HEAT_MAX 1.0e+04       /* 熱化の実行回数	 */
  #define AL 7.00e-02            /* 局所更新の初期最大幅 */
  #define AG 7.00e-02            /* 大域更新の初期最大幅 */
  #define SAMP_N 5               /* サンプリングの試行回数 */
  #define MAX 1000               /* シンプソン積分での刻み幅 */
  #define DIGIT 1.0e-04          /* ポテンシャルの計測間隔 */
  #define THREADS 10             /* 交換を行う幅（13が限界） */
  #define ADJUST_MAX 10          /* 熱化を行う回数 */
  #define RANGE 1.0

 /* 観測する物質（ヘリウム）に関わる定数 */
  #define NL 4                   /* 一辺あたりのセルの個数 */
  #define NUM (NL*NL*NL)    /* 粒子数に相当する */
  #define MASS 4.00              /* m [g/mol] */
  #define RHO 0.021843        /* n [個/Å^3] */
  #define RM 2.963               /* rm [Å] */
  #define EPSILON 10.948         /* ε [K] */
  #define A 1.8443101e+05
  #define D 1.4826
  #define ALPHA 10.43329537
  #define BETA -2.27965105
  #define C6 1.36745214
  #define C8 0.42123807
  #define C10 0.17473318

 /* Chin-Actionに必要な定数 */
  #define A1 0.33
  #define T0 0.082

 /* マクロ関数の定義 */
  #define ST(X1,X2,X3)     ( 0.5*Tn(X1,X2,X3)/(eta*tau*tau) )
  #define SU(X1,X2,X3,X4)  ( Un(X1,X2,X3,X4) )
  #define SF(X1,X2,X3,X4)  ( tau*tau*eta*Fn(X1,X2,X3,X4) )

 /* グローバル変数の定義 */
  extern int N;
  extern double eta,tau,L;

 /* 経路の情報を含む構造体の定義 */
  struct path {
    int num;
    double al;
    double ag;
    struct thread {
      struct part1 {
        int link;              /* 上のスレッドでの居場所 */
        double r[3];
        double f[3];
      } pe[NUM];
    } td[3*M];
    struct part2 {
      struct interaction {
        int nb;
      } pj[NUM];
    } pi[NUM];
  };

 /* プロトタイプ宣言 */
  void Initialize(path *);
  void Amplitude(path *);
  void Sampling(path *);
  void Flipfunc(path *,int *);
  void Periodic(double *,int *);
  void Adjustment(path *,int *,int *);
  void Kahan(double *,double *,double *);
  int Localflip(path *);
  int Globalflip(path *);
  int Exchange(path *);
  double Winding(path *);
  double Tn(path *,int *,int *);
  double Un(path *,int *,int *,int *);
  double Fn(path *,int *,int *,int *);
  double U(double *);
  double dU(double *);
  double Correction();
  double Random();


