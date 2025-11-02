
// ＜CLAPACKのサンプルプログラム＞ //
//
//    ( 12 0 4 )     ( 6 0 2 )
//  A=(  0 6 0 )   B=( 0 2 0 )  
//    (  4 0 4 )     ( 2 0 2 )
//
// に対して Ax=λBx の固有値 λ と固有ベクトル x を求める。
 
// ヘッダーの定義
#include <stdio.h>
#include <f2c.h>
#include <clapack.h>

// マクロの定義
#define N 3

// メイン関数
int main(void)
{

    /* 変数の定義 */
    char job,form;
    double A[N*N],B[N*N],E[N],work[3*N],M1,M2,M3;
    static long int i,type,n,Ac,Bc,lwork=3*N,info,lc;

    /* N×N行列のAに値を代入（999が入っている配列の値は利用されない） */
    A[0]= 12.0; A[3]=  0.0; A[6]=4.0; 
    A[1]=999.0; A[4]=  6.0; A[7]=0.0;
    A[2]=999.0; A[5]=999.0; A[8]=4.0;

    /* N×N行列のBに値を代入（999が入っている配列の値は利用されない） */
    B[0]=  6.0; B[3]=  0.0; B[6]=2.0;
    B[1]=999.0; B[4]=  2.0; B[7]=0.0;
    B[2]=999.0; B[5]=999.0; B[8]=2.0;

    /* 固有値方程式を計算 */
    type=1; job='V'; form='U'; n=N; Ac=n; Bc=n; 
    dsygv_(&type,&job,&form,&n,A,&Ac,B,&Bc,E,work,&lwork,&info);

    /* 計算結果を出力 */
    printf("\n");
    for(i=0;i<N;++i){
        lc=i*N; M1=A[lc]; M2=A[lc+1]; M3=A[lc+2];
        printf("  λ=%lf : vector=[ %9.6f %9.6f %9.6f ]\n",E[i],M1,M2,M3);
    }
    printf("\n");

    return 0;

}

