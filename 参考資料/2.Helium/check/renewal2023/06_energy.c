
 /* ヘッダーファイルの指定 */
  #include "00_header.h"

// 変位の二乗を求める関数（Cn:現在のスレッドでの番号、Cd:下のスレッド）
  double Tn(Path *path,int number,int thread)
 {
    /* 変数の定義 */
    int d,on;
    double answer,dr;

    /* Chin-Actionの係数 */
    constexpr double t1=0.5-T0;
    constexpr double ti[3]={1.0/t1,1.0/t1,0.5/T0};

    /*Tnの値を計算  */
    answer = 0.0; on=(thread+1)%(3*M);
    for(d=0;d<RANK;d++){
        dr = (path->pinfo[number].tinfo[on].r[d])-(path->pinfo[number].tinfo[thread].r[d]);
        periodic(&dr,1); answer += dr*dr;
    }
    answer *= 0.5*ti[thread%3]/(ETA*TAU*TAU);

    return answer;
 }

// 相互作用を計算する関数
  double Un(Path *path,int number,int thread)
 {
    /* 変数 */
    int i,d;
    double x,rij,answer;

    /* Chin-Actionの係数 */
    constexpr double t1=1.0-2.0*T0,v1=1.0/(6.0*t1*t1),v2=1.0-2.0*v1;
    constexpr double vi[3] = {v1,v2,v1};

    /* ポテンシャルの計算を行う */
    answer=0.0;
    for(i=0;i<N;i++){
        if(i == number) continue;
        rij=0.0;
        for(d=0;d<3;d++){
            x = path->pinfo[number].tinfo[thread].r[d]-path->pinfo[i].tinfo[thread].r[d];
            periodic(&x,1); rij += x*x;
        }
        
        rij = sqrt(rij); answer += U(rij);
    } 
    answer *= vi[thread%3];

    return answer;
 }

// 力に関する重みを計算する関数
 double Fn(Path *path,int number,int thread)
{
    /* 変数の定義 */
    int d;
    double answer;

    /* Chin-Actionで用いる定数 */
    constexpr double t1 = 1.0-2.0*T0;
    constexpr double v1 = 1.0/(6.0*t1*t1);
    constexpr double u0 = ( 1.0-(1.0/t1)+(v1/t1) )/12.0;
    constexpr double ai[3] = {A1,1.0-2.0*A1,A1};

    /* 初期化を行う */
    answer=0.0;
    for(d=0;d<RANK;d++){
        answer += (path->pinfo[number].tinfo[thread].f[d])*(path->pinfo[number].tinfo[thread].f[d]);
    }
    answer *= u0*TAU*TAU*ETA*ai[thread%3];

    return answer;
}

// ポテンシャルを計算する関数
  double U(double x)
 {
    /* ポテンシャルの計算 */
    double x2,x6,x8,x10,fx,u1,u2,u3,u4,value;

    x2=x*x; x6=x2*x2*x2; x8=x6*x2; x10=x8*x2;
    u1=A*exp(-ALPHA*x+BETA*x2); u2=C6/x6; u3=C8/x8; u4=C10/x10;

    if(x<D){
        value = ( (D/x)-1.0 ); 
        fx = exp(-value*value);
    }
    else{
        fx = 1.0;
    }

    value = u1-(u2+u3+u4)*fx;

    return value;
}
