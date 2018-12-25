// сия программа для отработки интерполяции аф по интегралам.
// интерполяция по интегралам изображения 9х9 отработана.
// 12.05.16

#include <iostream>
//#include <complex>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include "./image/image.h"
#include "./imageop/imageop.h"
#include "./up/up.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector.h>

//#include "defs.h"

long double testimage[]={
1.0000000025607425, 1.0000034575202832, 1.0002296211519204, 1.0018904349931375, 1.0023142471714943, 1.0004992417451533, 1.000014508934933,  1.0000000180353237, 1.,
1.0000034575202827, 1.0016632502727547, 1.039626961378202,  1.159235825072223,  1.183710906381325,  1.0671559030462219, 1.0045356086362005, 1.0000233550591353, 1.0000000008312357,
1.00022962115192,   1.039626961378202,  1.3916518416991812, 1.7760000905949596, 1.8187883375855227, 1.5188570417426293, 1.0913768176133072, 1.0014071283674104, 1.00000037564291,
1.0018904349931375, 1.1592358250722232, 1.7760000905949596, 1.991071162368386,  1.9963394706910373, 1.8923778295707876, 1.3199030875159177, 1.0110094784114048, 1.000008783990875,
1.0023142471714943, 1.183710906381325,  1.8187883375855232, 1.9963394706910373, 1.9993574803514214, 1.935098728892264,  1.3869059162726183, 1.0158564132099182, 1.000014927443585,
1.0004992417451533, 1.067155903046222, 1.5188570417426293,  1.8923778295707876, 1.935098728892264,  1.7106066847410506, 1.176184482317732,  1.0038892399418635, 1.0000017940912411,
1.0000145089349328, 1.0045356086362007, 1.091376817613307,  1.3199030875159177, 1.3869059162726183, 1.1761844823177319, 1.0188013458363614, 1.000163833123728,  1.0000000151678552,
1.0000000180353237, 1.0000233550591353, 1.0014071283674104, 1.0110094784114048, 1.0158564132099184, 1.0038892399418635, 1.0001638331237281, 1.000000467881383,  1.,
0.9999999999999999, 1.0000000008312357, 1.0000003756429099, 1.000008783990875, 1.000014927443585, 1.0000017940912411,   1.0000000151678552, 1.,                 1.
};//значения из математического пакета для сверки

using namespace std;
// размер хвостов
dword Cx_(dword Nx){return(Nx-1);}         
// количество членов полинома Mx=2*Nx-1
dword Mx_(dword Nx){return(Nx+Cx_(Nx));}        
//шаг базисных функций Hx=2*(Nx-1)/
long double Hx_(dword Nx){return((Nx+Cx_(Nx)-1)/(Mx_(Nx)-1));} 
//один хвост
dword cx_(dword Nx){return(Cx_(Nx)/2);}             
//порядок базисной функции
dword NUpx_(dword Nx){return((Nx-1)/2);}       

dword factorial(dword n)                //факториал
 {
  dword i=1;
  for(dword j=1; j<=n; j++) i=i*j;
  return i;
 }

dword binomial(dword n, dword m)        //вычисление биноминальных коэффициентов
  {
   dword c,i=1,j;
   if(m<(n-m)){c=m;m=n-m;}else c=n-m;   //m-большее, c-меньшее
   for(j=m+1;j<=n;j++) i=i*j;           //числитель
   c=factorial(c);                      //знаменатель
   c=i/c;
   return c;
  }

long double upintegrate(atomfunc upn1, dword n, long double a, long double b, long double h) //интеграл атомарной функции
//Передаём: upn1 - атомарная функция; n - порядок функции; a, b - пределы интегрирования; h - шаг интегрирования
  {
   long double s=0; //интегральная сумма сначала равна нулю
   for (long double i=a; i<b; i=i+h) s=s+upn1.upn(i+h/2,n); //метод прямоугольников 
   s=s*h; //умножаем сумму значений на шаг
   return s;
  }

/* не работает
//формирование массива кусков интегралов для текущей решётки дискрет
void iupinit(atomfunc upn1, dword n, dword Nx, long double h, long double* iup)
 {
   long double Hx=Hx_(Nx);
   for (dword i =0; i<(Nx+1)/(2*Hx)+1; i++)
    iup[i]=upintegrate(upn1,n,double(i)*Hx,(double(i)+1)*Hx,h);
 }*/

int UnitStep(long double x){return ((x<0)?0:1);} //функция ступеньки

long double intup(atomfunc upn1, dword n, dword Nx, long double ax,long double bx,dword i,long double h) //интеграл атомарной функции в полиноме
//Передаём: upn1 - атомарная функция; n - порядок функции; Nx - размер интерполируемого блока
//ax, bx - пределы интегрирования; i - номер базисной функции в полиноме (определяет сдвиг); h - шаг интегрирования
  {
     long double Hx=Hx_(Nx); //шаг базисных функций
     dword cx=cx_(Nx); //хвосты базисных функций
     long double sign,ival;
     if (ax<bx) sign=1; else {sign=ax;ax=bx;bx=sign;sign=-1;} //если пределы перепутаны, знак результата будет противоположным
     ival=upintegrate(upn1,n,ax-((i-1)*Hx-cx),bx-((i-1)*Hx-cx),h); //интегрируем сдвинутую функцию

/* не работает
// предварительное вычисление таблицы кусков интегралов, поэтому интегрирование
// атомарных функций выполняется без предварительной таблицы кусков и медленно.
     if (ax<=((i-1)*Hx-cx-(Nx+1)/2)&&(bx>=(i-1)*Hx-cx+(Nx+1)/2)) return 1; //если пределы интегрирования покрывают носитель
     if (bx<=((i-1)*Hx-cx-(Nx+1)/2)||(ax>=(i-1)*Hx-cx+(Nx+1)/2)) return 0; //если пределы интегрирования вне носителя функции
     for (dword k=dword(ceil(fmax((ax-((i-1)*Hx-cx))/Hx,-(Nx+1)/(2*Hx))))+1; k<=dword(floor(fmin((bx-((i-1)*Hx-cx)),(Nx+1)/(2*Hx)))); k++)
       { ival=ival+iup[abs(int(k)-UnitStep(-k))];} //суммируем входящие куски
     ival=ival+upintegrate(n,ax-((i-1)*Hx-cx),ceil(fmax((ax-((i-1)*Hx-cx))/Hx,-(Nx+1)/(2*Hx))),h); //добавляем хвосты, не составляющие целого куска слева
     ival=ival+upintegrate(n,floor(min((bx-((i-1)*Hx-cx)),(Nx+1)/(2*Hx))),bx-((i-1)*Hx-cx),h); //и справа 
*/
     return sign*ival; //возвращаем результат со знаком
  }

void upsxinit(atomfunc upn1, dword Nx, long double* upsx)// Инициализация матрицы A
//передаём upn1 - атомарная функция; Nx - размер блока для интерполяции; 
//получаем upsx - матрица A.
  {
    dword i,j,k,Mx=Mx_(Nx),NUpx=NUpx_(Nx);
//    cout << "binomials...";
    for(k=0; k<=(Mx-Nx)/2-1; k++)//k - номер строки, первые (Mx-Nx)/2 строк
      for (i=0; i<Mx; i++) //для каждого i-го элемента строки
        upsx[k*Mx+i]=binomial(k+NUpx+1,k+NUpx+1-i)*pow(float(-1),int(i));  //биномиальные коэффициенты
//    cout << "binomials finite differences...";
    for(k=0;k<=(Mx-Nx)/2-1;k++) //k - номер строки (порядок конечной разности)
      for (i=(Mx-Nx)/2-1;i>k;i--)//i - номер строки, первые (Mx-Nx)/2 строк
        for (j=0;j<Mx;j++)//для каждого j-го элемента строки
          upsx[i*Mx+j]=upsx[(i-1)*Mx+j]-upsx[i*Mx+j];//конечные разности биномиальных коэффициентов
    for(k=0; k<=(Mx-Nx)/2-1; k++)//k - номер строки, первые (Mx-Nx)/2 строк
      for (i=0; i<Mx; i++) //для каждого i-го элемента строки
        upsx[(k+(Mx+Nx)/2)*Mx+i]=upsx[k*Mx+i];
//    cout << "frequencies...";
    for(k=0;k<Mx-Nx-float(Mx-Nx+1)/2;k++) //k+(Mx-Nx)/2 - номер строки от (Mx-Nx)/2+1 до (Mx-Nx)
      for(i=0;i<Mx;i++) //для каждого i-го элемента строки
        upsx[(k+Nx)*Mx+i]=cos(i*M_PI*Nx*(1/float(Mx))*(2+k/float(Mx-Nx-float(Mx-Nx+1)/2)))+sin(i*M_PI*Nx*(1/float(Mx))*(2+k/float(Mx-Nx-float(Mx-Nx+1)/2)));//на оставшихся позициях прижимаем гармоники
//    cout << "atomic...\n";    
    for(k=0;k<Nx;k++)  //последние Nx строк
      for(i=0;i<Mx;i++) //заполняем значениями интегралов интерполяционного полинома
        upsx[k*Mx+i]=intup(upn1,NUpx,Nx,k,k+1,i+1,0.0001);
  }

void ridotsinit(Image im0, dword Nx, long double* ridots) //Инициализация вектора B
//Передаём im0 - прочитанное изображение; Nx - размер блока для интерполяции
//получаем ridots - вещественная матрица со значениями в пределах от 1 до 2.
  {
    dword i,j,Mx=Mx_(Nx);
    for(i=0;i<Mx*Nx;i++)ridots[i]=0;//Обнуляем значения ridots
    for(i=0;i<Nx;i++)
      for(j=0;j<Nx;j++)
        ridots[i*Nx+j]=double(im0.getrpixel(i,j))/255+1; //читаем изображение (монохромное) и нормируем в диапазон от 1 до 2.
//        ridots[i*Nx+j]=testimage[i*Nx+j];
  }

void findDDD1(long double* upsy, long double* idots, dword Nx, long double* DDDY) //Вычисление коэффициентов полинома по строкам.
//Передаём upsy - матрица A; idots - вектор B; Nx - размер блока для интерполяции
//Получаем DDDY - матрица с коэффициентами полинома по строкам.
//Для вычисление используем gsl
  {
    dword i,j,k,Mx=Mx_(Nx);
    int s;
    gsl_matrix* A = gsl_matrix_alloc(Mx,Mx); //Матрица A
    gsl_vector* B = gsl_vector_alloc(Mx); //Вектор B
    gsl_vector* X = gsl_vector_alloc(Mx); //Вектор X
    gsl_permutation* P = gsl_permutation_alloc(Mx); //Объект для LU декомпозиции
    for(k=0;k<Nx;k++)                            //Для каждой строки изображения
      {
        for(i=0;i<Mx;i++)for(j=0;j<Mx;j++)gsl_matrix_set(A,i,j,upsy[i*Mx+j]); //Задаём A
        for(i=0;i<Mx;i++)gsl_vector_set(B,i,idots[i*Nx+k]); //Задаём B
        gsl_linalg_LU_decomp(A,P,&s); //Выполняем LU разложение матрицы A
        gsl_linalg_LU_solve(A,P,B,X); //получаем решение СЛАУ
        for(i=0;i<Mx;i++)DDDY[k*Mx+i]=gsl_vector_get(X,i); //Записываем результат в строчку DDDY
      }
    gsl_permutation_free (P);
    gsl_vector_free (X);
    gsl_vector_free (B);
    gsl_matrix_free (A);
  }

void findDDD2(long double* upsx, long double* DDDY, dword Nx, long double* DDDX) //Вычисление коэффициентов 2D полинома.
//Передаём upsx - матрица A; DDDY - вектор B; Nx - размер блока для интерполяции
//Получаем DDDX - матрица с коэффициентами 2D полинома
//Для вычисление используем gsl
  {
    dword i,j,k,Mx=Mx_(Nx);
    int s;
    gsl_matrix* A = gsl_matrix_alloc(Mx,Mx); //Матрица A
    gsl_vector* B = gsl_vector_alloc(Mx); //Вектор B
    gsl_vector* X = gsl_vector_alloc(Mx); //Вектор X
    gsl_permutation* P = gsl_permutation_alloc(Mx); //Объект для LU декомпозиции
    for(k=0;k<Mx;k++)  //для каждого столбца коэффициентов, полученных для строк
      {
        for(i=0;i<Mx;i++)for(j=0;j<Mx;j++)gsl_matrix_set(A,i,j,upsx[i*Mx+j]); //Задаём A
        for(i=0;i<Nx;i++)gsl_vector_set(B,i,DDDY[i*Mx+k]); //Задаём B из коэффициентов, полученных для строк
        for(i=Nx;i<Mx;i++)gsl_vector_set(B,i,0); //Задаём B, недостающие элементы нулевые (условия по производным и гармоникам)
        gsl_linalg_LU_decomp(A,P,&s); //Выполняем LU разложение матрицы A
        gsl_linalg_LU_solve(A,P,B,X); //получаем решение СЛАУ
        for(i=0;i<Mx;i++)DDDX[k*Mx+i]=gsl_vector_get(X,i);// Записываем результат в столбец DDDX
      }
    gsl_permutation_free (P);
    gsl_vector_free (X);
    gsl_vector_free (B);
    gsl_matrix_free (A);
  }

long double ifupf1 (atomfunc upn1, long double x, long double y, dword Nx, long double* DDDX) //Вычисление значения 2D полинома
//Передаём upn1 - атомарная функция; x,y - значения аргументов, Nx - размер блока для интерполяции;
//DDDX - матрица с коэффициентами 2D полинома
  {
    dword i,j,Mx=Mx_(Nx),NUpx=NUpx_(Nx);
    long double Hx=Hx_(Nx),cx=cx_(Nx);
    long double sx=0,sy,u;
    for(i=0;i<Mx;i++) //в этом цикле вычисляется значение двумерной свёртки для точки (x,y)
      {
        sy=0;
        for(j=0;j<Mx;j++) sy=sy+DDDX[j*Mx+i]*upn1.upn(y+cx-double(j)*Hx,NUpx);
        sx=sx+sy*upn1.upn(x+cx-double(i)*Hx,NUpx);
      }
    return sx;
  }

void showiup(dword Nx, long double* iup) //показывает значения кусков интегралов для диагностики
  {
    long double Hx=Hx_(Nx);
    for (dword i =0; i<(Nx+1)/(2*Hx)+1; i++) cout << iup[i]<<",";
    cout <<"\n";
  }

void showupsx(dword Nx, long double* upsx) //показывает матрицу A для диагностики
  {
     dword Mx=Mx_(Nx);
     for(int i=0;i<Mx;i++)
       {
         for (int j=0;j<Mx;j++) cout << upsx[i*Mx+j] <<",";
         cout <<"\n";
       }
  }

void showidots(dword Nx, long double* idots) //Показывает вектор B для диагностики
  {
     dword Mx=Mx_(Nx);
     for(int i=0;i<Mx;i++)
       {
         for (int j=0;j<Nx;j++) cout << idots[i*Nx+j] <<",";
         cout <<"\n";
       }
  }

void showDDD1(dword Nx, long double* DDD1) //Показывает матрицу коэффициентов по строкам для диагностики
  {
     dword Mx=Mx_(Nx);
     for(int i=0;i<Nx;i++)
       {
         for (int j=0;j<Mx;j++) cout << DDD1[i*Mx+j] <<",";
         cout <<"\n";
       }
  }

void showDDD2(dword Nx, long double* DDD2) //Показывает матрицу коэффициентов многочлена для диагностики
  {
     dword Mx=Mx_(Nx);
     for(int i=0;i<Mx;i++)
       {
         for (int j=0;j<Mx;j++) cout << DDD2[i*Mx+j] <<",";
         cout <<"\n";
       }
  }

void showf1(atomfunc upn1, dword Nx, long double* DDD2) //Вычисляет таблицу значений интерполяционного многочлена для диагностики
  {
     for(int i=0;i<Nx;i++)
       {
         for (int j=0;j<Nx;j++) cout << ifupf1(upn1, double(i), double(j), Nx, DDD2) <<",";
         cout <<"\n";
       }
  }

void savef1(atomfunc upn1, dword Nx, char* filnam, long double* DDD2) //Сохранение значений полинома в файл
//Передаём upn1 - атомарная функция; Nx - размер блока для интерполяции; filnam - имя файла; DDD2 - матрица с коэффициентами 2D полинома
  {
    Image i2(Nx,Nx);        //здесь формируем изображение
    for(int i=0;i<Nx;i++)   //для каждой строки
      for(int j=0;j<Nx;j++)  //для всех элементов строки
        i2.putpixel(i,j,0,0,round((ifupf1(upn1, double(i), double(j), Nx, DDD2)-0.5)*128)); //выводим изображение синим, при этом немного нормируем
    i2.saveimage(filnam);    //Сохраняем изображение в файл
  }

void saveridots(dword Nx, char* filnam, long double* img) //Сохранение вектора B в файл для диагностики
  {
    Image i2(Nx,Nx);
    for(int i=0;i<Nx;i++)
      for(int j=0;j<Nx;j++)
        i2.putpixel(i,j,0,round((img[i*Nx+j]-0.5)*128),0);
    i2.saveimage(filnam);      
  }

void saveupsx(dword Nx, long double* upsx) //Сохранение матрицы A в файл для диагностики
  {
    dword Mx=Mx_(Nx);
    Image i2(Mx,Mx);
    long double umax = upsx[0];
    long double umin = upsx[0];
    long double k;
    for(int i=1;i<Mx*Mx;i++){
      if(umax<upsx[i]) umax=fabs(upsx[i]);                  //поиск максимума
      if(umin>upsx[i]) umin=fabs(upsx[i]);                  //поиск минимума
      }
    cout << " max=" << umax << "; min=" << umin << "\n";                     //минимум и максимум
    for(int i=0;i<Mx;i++)
      for(int j=0;j<Mx;j++){
        k=upsx[j*Mx+i]*1024/(umax-umin);          //вычисляем нормированное значение
        if(k>0) 
          if(k<256) i2.putpixel(i,j,0,int(k),0); //чёрный зеленеет
          else if(k<512) i2.putpixel(i,j,int(k-255),255,0); //зелёный желтеет
               else if(k<768) i2.putpixel(i,j,255,int(767-k),0); //жёлтый краснеет
                    else i2.putpixel(i,j,255,int(k-767),int(k-767)); //красный СВЕТЛЕЕТ
        else
          if(fabs(k)<256) i2.putpixel(i,j,0,int(fabs(k)),int(fabs(k))); //чёрный голубеет
          else if(fabs(k)<512) i2.putpixel(i,j,0,int(511-fabs(k)),255); //голубой синеет
               else if(fabs(k)<768) i2.putpixel(i,j,int(fabs(k)-511),0,255); //синий розовеет
                    else i2.putpixel(i,j,int(1023-fabs(k)),0,int(1023-fabs(k))); //розовый темнеет
               
      }
    i2.saveimage("UPsx.bmp");      
  }

void usage ()                                                                   //help screen
 {
  cout << "atomic interval interpolation by Vlad Efremof, v2@bk.ru" << "\n";   //
  cout << "usage: <this processor> ???" << "\n"; //
 }

int main(int argc, char *argv[])
 {
  int t1=clock();cout<<t1<<" main started;\n";                  //время старта
//читаем файл
  if (argc !=3) { usage(); return 100; };
  char* inbmpname = argv[1];
  char* outbmpname = argv[2];
  Image Im0(inbmpname);
//  Image Im0("w.bmp");
//основные структуры данных
  dword iw=Im0.getwidth();//ширина
  dword ih=Im0.getheight();//высота
  dword Nx=(iw<ih)?iw:ih ;// размер блока для интерполяции
  long double Hx=Hx_(Nx);//вычисляем шаг базисных функций
  dword Mx=Mx_(Nx);//Количество базисных функций в составе полинома
  cout <<"iw="<<iw<<",ih="<<ih<<",Nx="<<Nx<<",Hx="<<Hx<<",Mx="<<Mx<<"\n";
  atomfunc upn1(Nx,1000);  //базисная функция
//  long double *iup;          //куски интегралов
  long double *upsx;         //матрица А
  long double *ridots;        //массив векторов B из изображения.
  long double *DDD1;          //коэффициенты после первого прохода
  long double *DDD2;          //коэффициенты двумерного полинома
//  iup=new long double[dword((Nx+1)/(2*Hx))];
  upsx=new long double[Mx*Mx];
  ridots=new long double[Mx*Nx];
  DDD1=new long double[Nx*Mx];
  DDD2=new long double[Mx*Mx];

/* не работает
// предварительное вычисление таблицы кусков интегралов, поэтому интегрирование
// атомарных функций выполняется без предварительной таблицы кусков и медленно.
  iupinit(Nx,0.001);
  cout<<clock()<<" table iup of integrals done;\n";
  showiup();*/

  upsxinit(upn1,Nx,upsx); //в upsx матрица A атомарных функций для СЛАУ для построения атомарного многчлена
//  showupsx(Nx,upsx); saveupsx(Nx,upsx);
  cout<<clock()<<" table upsx of coefficients done;\n";
  ridotsinit(Im0,Nx,ridots);//в ridots вектор B для СЛАУ
//  showidots(Nx,ridots); saveridots(Nx,"ridots.bmp",ridots);
  cout<<clock()<<" table idots of right parts done;\n";
  findDDD1(upsx,ridots,Nx,DDD1); //построение многочлена решением СЛАУ по строкам
//  showDDD1(Nx,DDD1);
  cout<<clock()<<" table DDDY done;\n";
  findDDD2(upsx,DDD1,Nx,DDD2);  //построение 2D многочлена решением СЛАУ по столбцам
//  showDDD2(Nx,DDD2);
  cout<<clock()<<" table DDDX done;\n";
  atomfunc upn2(Nx,1000);  //базисная функция (почему-то upn1 после вызова findDDD1 портится)
//  showf1(upn2,Nx,DDD2); //таблица точек, вычисленных с помощью интерполяционного многочлена
  savef1(upn2,Nx,outbmpname,DDD2); //сохранение изображения, вычисленного с помощью интерполяционного многочлена
  int t2=clock(); t1=t2-t1; cout << t2 << " main ended, work TIME = " << t1 << "\n";        //время работы
  delete [] upsx; //убираем за собой
  delete [] ridots;
  delete [] DDD1;
  delete [] DDD2;
  return 1;
 }
