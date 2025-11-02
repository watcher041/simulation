
/* ヘッダーファイルの指定 */
#include"00_PI.h"

/* グローバル変数 */
double TAU = 1.0/( T*(double)M );

/* プログラムの全体像 */
int main(void)
{
    /* 経路を形を記録する配列を定義 */
    double r[3*M];

    /* 乱数の初期化 */
    srand48((unsigned) time_t(NULL));

    /* 理論値を計算 */
    double x=0.5/T;
    double E=0.5/tanh(x);
    double C=1.0/( 4.0*T*T*sinh(x)*sinh(x) );

    /* シミュレーションに関わる値を表示 */
    printf("\n");
    printf("（シミュレーションに関わる定数）\n");
    printf("　温度 T=%f\n",T);
    printf("　分割幅 τ=%f\n",TAU);
    printf("　理論値 E=%e C=%e\n",E,C);

    /* シミュレーションを実行 */
    simulation(r);

    /* 正常終了 */
    return 0;
}