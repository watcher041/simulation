
// 「ポテンシャルに関係する関数一覧」 //
#include"00_PI.h"

// 分子の変位の二乗が含まれる
 double Tn(double *r,int tn)
{
   int on=(tn+1)%(3*M);
   double answer,dr;
   
   /* Chin-Actionに用いる係数 */
   constexpr double t1 = 0.5-T0;
   constexpr double ti[3]={(1.0/t1),(1.0/t1),(0.5/T0)};

   dr=(r[on]-r[tn]);
   answer=0.5*ti[tn%3]*dr*dr/(TAU*TAU);

   return answer;
}

 double Un(double *r,int tn)
{
   double answer;
   
   /* Chin-Actionに用いる係数 */
   constexpr double t1 = 1.0-2.0*T0;
   constexpr double v1 = 1.0/(6.0*t1*t1), v2 = 1.0-2.0*v1;
   constexpr double vi[3] = {v1,v2,v1};

   answer=vi[tn%3]*U(r[tn]);

   return answer;
}

 double Fn(double *r,int tn)
{
   double answer;
   
   /* Chin-Actionに用いる係数 */
   constexpr double t1 = 1.0-2.0*T0, v1 = 1.0/(6.0*t1*t1);
   constexpr double ai[3] = {A1,(1.0-2.0*A1),A1};
   constexpr double u0 = ( 1.0-(1.0/t1)+(v1/t1) )/12.0;

   answer=u0*TAU*TAU*ai[tn%3]*(r[tn]*r[tn]);

   return answer;
}

// ポテンシャルエネルギーの形
 double U(double x)
{
   double answer;

   answer=0.5*(x*x);

   return answer;
}


