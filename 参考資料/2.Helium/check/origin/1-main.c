

 /* ヘッダーファイルの指定 */
  #include "pimc.h"
  #include <time.h>

 /* グローバル変数の定義 */
  int N;
  double L,eta,tau;

// メイン関数（シミュレーションの全体像）
  int main(void)
 {

    /* ----------------- シミュレーションの設定 -------------------------- */

    /* 物理定数 */
    const double kB=1.38064852e-23;    /* ボルツマン定数 kB[J/K] */
    const double amu=1.6605402e-27;    /* 原子質量単位 amu[kg/個] */
    const double hbar=1.0545717e-34;   /* ディラック定数 [J・s] */

    /* 変数の定義 */
    double eps,m,rm;

    /* 変数の初期化 */
    eps=EPSILON*kB; m=MASS*amu; rm=RM*1.0e-10;

    /* グローバル変数の初期化 */
    N=NUM; L=pow( (NUM/RHO) , (1.0/3.0) )/RM;
    tau=(EPSILON/TEMP)*( 1.0/(double)M ); eta=(hbar*hbar)/(m*rm*rm*eps);

    /* 乱数の初期化 */
    srand48((unsigned) time(NULL));

    /* 設定の確認 */
    printf("\n");
    printf("　温度 %f [K]\n",TEMP );
    printf("　分割幅 %f [K-1]\n",(tau/EPSILON) );
    printf("　数密度 %f [Å-3]\n",RHO);
    printf("　総粒子数 %d\n",N);
    printf("　スレッド数 %d\n",3*M);
    printf("　量子パラメータ %f\n",eta);
    printf("　セル一辺の長さ %f\n",L);
    printf("\n");

    /* ----------------------------------------------------------------- */

    /* 経路の情報を保存する構造体の定義 */
    struct path ps;

    /* 経路の情報を初期化 */
    Initialize(&ps);

    /* フリップ幅が5割になるように最大幅を設定 */
    Amplitude(&ps);

    /* 期待値を求める */
    Sampling(&ps);

    return 0;
 }



