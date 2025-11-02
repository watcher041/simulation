
#include"00_PI.h"

// 経路の初期化を行う関数
void reset(double *r)
{
    /* ループに利用する変数 */
    int n;

    /* 位置の初期化 */
    for(n=0;n<3*M;n++){
        r[n]=0.0;
    }
}