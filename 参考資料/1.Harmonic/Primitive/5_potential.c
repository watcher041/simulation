

// 「ポテンシャルに関係する関数一覧」 //
 #include "PI.h"

// 分子の変位の二乗が含まれる
 double Tn(double *r,int tn)
{
   int on=(tn+1)%M;
   double answer,dr;

   dr=( r[on]-r[tn] );
   answer=dr*dr;

   return answer;
}

 double Un(double *r,int tn)
{
   double answer;

   answer=U(r[tn]);

   return answer;
}

// ポテンシャルエネルギーの形
 double U(double x)
{
   double answer;

   answer=0.5*(x*x);

   return answer;
}


