
// ヘッダーファイルの指定
  #include <stdio.h>
  #include <math.h>
  #include <stdlib.h>

// シミュレーションに関する定数
  #define TAU 0.1              /* 分割幅 */
  #define M 50                 /* 分割数 */
  #define RANK 3               /* 次元数 */
  #define MCS_MAX 1.0e+04       /* モンテカルロステップ */
  #define HEAT_MAX 1.0e+04      /* 熱化の実行回数	 */
  #define SAMP_N 10              /* サンプリングの試行回数 */

// ポテンシャル（Aziz型）に関する定数
  #define A 1.8443101e+05
  #define D 1.4826
  #define ALPHA 10.43329537
  #define BETA -2.27965105
  #define C6 1.36745214
  #define C8 0.42123807
  #define C10 0.17473318

// Chin-Actionに関する定数
  #define A1 0.33
  #define T0 0.082

// 観測する物質に関する定数
  #define NX 3             /* 一辺当たりの粒子数 */
  #define RHO 0.02186      /* 数密度 [個/Å^3] */
  #define MASS 4.00        /* 質量 [g/mol] */
  #define RM 2.963         /* Azizパラメータ rm [Å] */
  #define EPS 10.948       /* Azizパラメータ ε [K] */
  #define NUM (NX*NX*NX)   /* struct内でN配列を利用するために定義 */

// 物理定数
  #define KB 1.38064852e-23  /* ボルツマン定数 kB[J/K] */
  #define AMU 1.6605402e-27  /* 原子質量単位 amu[kg/個] */
  #define HBAR 1.0545717e-34 /* ディラック定数[J・s] */

// マクロ関数の定義 
  #define Random()      ( drand48() )

// 計算して求める定数は変数にする
  extern int N;            /* ループに必要 */
  extern double L;         /* 周期境界条件適用などで必要 */
  extern double ETA,T;     /* 重み計算に必要 */

// 経路の情報を記録するための構造体
// まとめて代入するために中に構造体を作っている
  struct Path{
    struct {
      int link;
      struct {
          double r[RANK];
          double f[RANK];
      } tinfo[3*M];
    } pinfo[NUM];
  };

// プロトタイプ宣言
  void reset(Path *);
  void sampling(Path *);
  void periodic(double*,int);
  void force(Path *,int,int);
  void fij(Path *,int,int,int,int);
  void flipfunc(Path *,int);
  void localflip(Path *);
  void globalflip(Path *);
  double Tn(Path *,int,int);
  double Un(Path *,int,int);
  double Fn(Path *,int,int);
  double U(double);
  double dU(double);