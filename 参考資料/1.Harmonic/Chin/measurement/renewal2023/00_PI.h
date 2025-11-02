
// ヘッダーファイルの指定
#include <stdio.h>
#include <math.h>

// シミュレーションに関わる定数
#define TAU 1.0             /* 分割幅 */
#define M 5                 /* 分割数 */
#define HEAT_MAX 1.0e+04    /* 熱化の実行回数（10^4だけ行わないと、フリップ割合が収束しない） */
#define MCS_MAX 1.0e+06     /* サンプル数 */
#define N_SAMP 10           /* サンプリングの試行回数 */

// Chin-Actionのパラメータ
#define A1 0.33
#define T0 0.12

// 期待値などを計算するときに用いるマクロ関数
#define Random() drand48()

// グローバル変数
extern double T;

// プロトタイプ宣言
 void reset(double *);
 void simulation(double *);
 void flipfunc(double *,int);
 void localflip(double *);
 void globalflip(double *);
 double Tn(double *,int);
 double Un(double *,int);
 double Fn(double *,int);
 double U(double);
