

// ＜有限要素法（FEM法）によるシミュレーション＞ //

/* 一般化固有値方程式は、CLAPACK（シーレパック）により解く。 */

// 連続的に広範囲で見れると良い。

/* 原子ごとのパラメータ一覧（LJポテンシャルの場合）

   　　 σ[nm]   ε[K]    m[amu]
　・Ne… 0.274  36.2149  20.1797
　・Ar… 0.340  120.958  39.948
　・Kr… 0.365  162.967  83.798
　・Xe… 0.398  231.775  131.293

*/

// ヘッダーファイル
 #include <stdio.h>
 #include <math.h>
 #include <stdlib.h>
 #include <f2c.h>
 #include <clapack.h>

// 観測に関する物理量
 #define N 500          /* 分割数N */
 #define Ns 2           /* 始点とする地点（Xs=Ns*h） */
 #define Xf 30.0        /* 終点の位置（σでスケーリング） */
 #define T_MIN 0.01     /* 観測する温度の最低値 */
 #define T_MAX 1.0      /* 観測する温度の最大値 */
 #define DIGIT 0.01     /* 観測する温度の間隔*/
 #define Z 2            /* 物質の番号 */

 double eta;

// プロトタイプ宣言
 void Matrix(double *,double *);
 double V(double *);

// メイン関数
 int main(void)
 {
    int i;
    char job,triangle;
    long int type,n,line,column,lw,info;
    double A[N*N],B[N*N],Ei[N],w[3*N];
    double T,beta,numerator,denominator,eigen,E,EE,C;
    FILE *EIGEN_FILE,*DATA_FILE;
    
    /* 物理定数の定義(kittel参照、圧力は不明) */
   const double kB = 1.38064852e-23;    /* ボルツマン定数 kB[J/K] */
   const double amu = 1.6605402e-27;    /* 原子質量単位 amu[kg/個] */
   const double hbar = 1.0545717e-34;   /* ディラック定数 [J・s] */
   const double mass[5] = {4.00,20.2,39.9,83.8,131};
   const double sigma[5] = {2.56,2.74,3.40,3.65,3.98};
   const double epsilon[5] = {14,50,167,225,320};
   const double m = mass[Z-1]*amu;            /* m [kg/個] */
   const double sig = sigma[Z-1]*1.0e-10;     /* σ [Å] */
   const double eps = epsilon[Z-1]*1.0e-23;   /* ε [J] */

   /* 無次元量ηの値を求める。（どれだけ、量子的なものかを表す） */
   eta=(hbar*hbar)/(0.5*m*sig*sig*eps);  /* 慣性質量でスケーリングするため */

    /* 変数の初期化 */
    type = 1; job = 'V'; triangle = 'U';
    n = N; line= n; column = n;  /* 正方行列の次元数 */
    lw = 3*n;   

    /* 固有値方程式を解く */
    Matrix(A,B); 
    dsygv_(&type,&job,&triangle,&n,A,&line,B,&column,Ei,w,&lw,&info);

    /* 出力ファイルが残っていれば削除する */
    remove("fem.dat"); remove("eigen.dat");

    /* エネルギー固有値を出力 */
    EIGEN_FILE=fopen("eigen.dat","w");
    for(i=0;i<N;i++){
        fprintf(EIGEN_FILE,"%3d %3e\n",i,Ei[i]*(eps/kB));
    }
    fclose(EIGEN_FILE);
  
    printf( "\n η  = %e \n", eta);
    printf( " ε  = %e K\n",eps/kB );
    printf( " E0 = %e K\n\n",Ei[0]*(eps/kB) );
  
   for(T=T_MIN;T<=T_MAX;T+=DIGIT){

      numerator=0.0; denominator=0.0; beta=1.0/T;              
      for(i=0;i<N;i++){  
          eigen=Ei[i]-Ei[0];
          numerator+=eigen*exp(-beta*eigen);
          denominator+=exp(-beta*eigen);     
      }
      E=numerator/denominator;
    
      numerator=0.0; 
      for(i=0;i<N;i++){
          eigen=Ei[i]-Ei[0];
          numerator+=eigen*eigen*exp(-beta*eigen);
      }
      EE=numerator/denominator; 
      C=beta*beta*(EE-E*E);
    
      DATA_FILE=fopen("fem.dat","a"); 
      fprintf(DATA_FILE,"%f %e %e\n",T*(eps/kB),(E+Ei[0])*(eps/kB),C);
      fclose(DATA_FILE);
    
     printf( " T = %f K\n",T*(eps/kB) );
     printf( " E = %e K\n",(E+Ei[0])*(eps/kB));
     printf( " C = %e \n\n",C );

    }
    return 0;
 }

// シュレーディンガー方程式の演算子にあたる行列を作成
 void Matrix(double *A,double *B)
 {

    int i,j,k;
    double h,inv_h,x_i,x_m,x_p,x_j,f1,f2,f3;

    h=(double)Xf/(double)N;
    inv_h=1.0/h;

    for(i=0;i<N;i++){
        x_i=h*(double)(i+Ns); /* 始点は、ポテンシャルが有限なるところから */
      for(j=0;j<N;j++){
        if(j==i){
            x_m=h*(double)(i-1+Ns);
            x_p=h*(double)(i+1+Ns);
            f1=2.0;
            f2=V(&x_m)+6.0*V(&x_i)+V(&x_p);
            f3=4.0;
        }else if(j==i-1 || j==i+1){
            x_j=h*(double)(j+Ns);
            f1=-1.0;
            f2=V(&x_i)+V(&x_j);
            f3=1.0;      
        }else{
            f1=0.0;
            f2=0.0;
            f3=0.0;
        }
        k=i*N+j;
        A[k]=0.5*eta*f1*(inv_h)+(h*f2/12.0);
        B[k]=h*f3/6.0;      
      }
    }

 }

// ポテンシャルエネルギー（*x：分子間の距離）
 double V(double *x)
 {
    double answer,r2,r6,R;    
    
    r2=(*x)*(*x);
    r6=r2*r2*r2; R=1.0/r6;
    answer=4.0*R*(R-1.0);
    
    return answer;
 }


