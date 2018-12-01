#include "up.h"
#include <math.h>

long double sinc(long double x) {if(x==0)return(1); else return(sin(x)/x);}

//обычная реализация атомарных функций и спектров

long double fup(long double u, word m)
 {
  long double f=sinc(u/2);
  if(m>1023){m=1023;}
  for(word i=2; i<=m; i++) {f=f*sinc(u/pow(2,i));}
  return f;
 }

long double up(long double x, dword s, word m)
 {
  if((x<=-1)||(x>=1))return(0);
  if(x==0)return(1);
  long double y=0.5;
  for(dword k=1; k<=s; k++) {y=y+(fup(M_PI*double(2*k-1),m)*cos(M_PI*double(2*k-1)*x));}
  return y;
 }

long double fupn(long double u, dword n, word m)
 {
  long double f=pow(sinc(u/2),n);
  if(m>1023){m=1023;}
  for(word i=1; i<=m; i++) {f=f*sinc(u/pow(2,i));}
  return f;
 }

long double upn(long double x, dword n, dword s, word m)
 {
  long double l=(double(n)/2+1);
  if((x<=-l)||(x>=l))return(0);
  long double y=0.5;
  for(dword k=1; k<=s; k++) {y=y+fupn(2*M_PI*double(k)/double(2+n),n,m)*cos(2*M_PI*double(k)*x/double(2+n));}
  y=y*2/double(n+2);
  return y;
 }

//реализация объекта с upn и таблицей

atomfunc :: atomfunc(const dword n1, word np1)                                              //конструктор списка простых сомножителей.
 {
  n=n1; //максимальный порядок функции
  np=(np1<4096)?np1:4096;//количество множителей
  nc=2*n; //количество строк таблицы
  ns=(double(n*2)+1/4<16)?n*2:16;   //..количество слагаемых в каждой строке
  fupna=new long double[nc*ns];  //массив коэффициентов
  for(dword i=1;i<=nc;i++)  //для каждой строки, означающей порядок функции
    for(dword j=1;j<=ns;j++) //для каждого элемента строки - слагаемого
      fupna[(i-1)*ns+j-1]=fupn((2*M_PI*double(j))/double(2+i),i);
 }

//fupn(2*M_PI*double(k)/double(2+n),n)

long double atomfunc::fupn(long double u, dword n)
 {
  long double f=pow(sinc(u/2),n);
  for(word i=1; i<=np; i++) {f=f*sinc(u/pow(2,i));}
  return f;
 }

//fupna[(n-1)*ns+k-1]

long double atomfunc::upn(long double x, dword n)
 {
  long double l=(double(n)/2+1); //границы носителя функции
  if((x<=-l)||(x>=l))return(0); //если за пределами, там 0
  long double y=0.5; //нулевой член ряда
  for(dword k=1; k<=ns; k++) //для каждого элемента n-й строки (соотв. порядку)
    y=y+fupna[(n-1)*ns+k-1]*cos(2*M_PI*double(k)*x/double(2+n));
  y=y*2/double(n+2);
  return y;
 }


/*long double up(double x, dword m)
 {
  if(x>=1)return(0);
  if(x<=-1)return(0);
  if(x>=0.5)return(up_table[int(99-x*100)]);
  if(x<=-0.5)return(up_table[int(99+x*100)]);
  if(x>0)return(1-up_table[int(x*100)]);
  if(x<0)return(1-up_table[int(-x*100)]);
  if(x==0)return(1);
 }*/
 
/*long double ma_up(dword n)
 {
  long double a=0;
  long double a1=0;
  if (n==0) return 1;
  if (n&1!=0) return 0;
  a1=factorial(n); a1=a1/((2<<(n-1))-1);
  for (dword k=1; k<n/2+1; k++) 
   {a=a+ma_up(n-2*k)/(factorial(n-2*k)*factorial(2*k+1));}
  a=a1*a;
  return a;
 }

long double mb_up(dword n)
 {
  long double b;
  if (n&1==0) b=ma_up(n)/2;
  else switch(n)
   {
    case 1: b=0.1388888888889; break;
    case 3: b=0.02648148148148; break;
    case 5: b=0.0005043426835755; break;
    case 7: b=0.0006210698819994; break;
    default: return 0;
   }
  return b;
 }
*/
