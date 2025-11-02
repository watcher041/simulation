
#include "00_header.h"

// グローバル定数を初期化
  int N = NUM;   
  double L = pow( ((double)N/RHO) , (1.0/3.0) )/RM;
  double TAU = EPS / (T*(double)M);
  double ETA = (HBAR*HBAR)/(MASS*AMU*(RM*1.0e-10)*(RM*1.0e-10)*EPS*KB);

// メイン関数（プログラムの全体像）
  int main(void)
{
    /* 経路情報を記録する構造体を定義 */
    Path path;

    /* 乱数の初期化 */
    srand48((unsigned) time_t(NULL));

    /* 設定の確認 */
    printf("\n");
    printf("　温度 %f [K]\n",T);
    printf("　分割幅 %f [K^-1]\n",(TAU/EPS) );
    printf("　数密度 %f [Å^-3]\n",RHO);
    printf("　総粒子数 %d\n",N);
    printf("　スレッド数 %d\n",3*M);
    printf("　量子パラメータ %e\n",ETA);
    printf("　セル一辺の長さ %f\n",L);
    printf("\n");

    /* 経路の熱化後にサンプリングを行う */
    sampling(&path);

    /* 正常終了 */
    return 0;
}