#include "../image/defs.h"
#include "../image/image.h"
#include "imageop.h"
#include <iostream>                                                             //подключаем библиотеку для работы с потоками
#include <stdlib.h>

using namespace std;

Image grayscale(Image i1, double rc=0.11, double gc=0.59, double bc=0.3)        //сделать черно-белым.
 {
  byte scale; dword i,j, iheight = i1.getheight(), iwidth = i1.getwidth(); 
  Image i2(i1);
  for(i=0; i<iwidth; i++) for(j=0; j<iheight; j++)
     {
      scale = byte(i1.getbpixel(i,j)*0.11+i1.getrpixel(i,j)*0.59+i1.getgpixel(i,j)*0.3);
      i2.putpixel(i,j,scale,scale,scale);
     }
  return i2;
 }

Image I_conv(Image i1, Image conv)                                      //дискретная свертка
 {
  dword iheight = i1.getheight(), iwidth = i1.getwidth();                  //высота и ширина изображения
  dword cheight = conv.getheight(), cwidth = conv.getwidth();              //высота и ширина сверточной последовательности
  dword theight = iheight + cheight, twidth = iwidth + cwidth;             //
  Image t1(twidth, theight), i2(i1);
  dword sumr,sumg,sumb,i,j,k,l,rweight=0,gweight=0,bweight=0;
  for(i=0; i<cwidth; i++) for(j=0; j<cheight; j++)
   { rweight=rweight+conv.getrpixel(i,j); gweight=gweight+conv.getgpixel(i,j); bweight=bweight+conv.getbpixel(i,j); }
  if(rweight==0) rweight=1; if(gweight==0) gweight=1; if(bweight==0) bweight=1;
  for(i=0; i<twidth; i++) for(j=0; j<theight; j++) t1.putpixel(i,j,0,0,0);
  for(i=0; i<iwidth; i++) for(j=0; j<iheight; j++)
    t1.putpixel(i+cwidth/2,j+cheight/2,i1.getrpixel(i,j),i1.getgpixel(i,j),i1.getbpixel(i,j));
  for(i=0; i<iwidth; i++) for(j=0; j<iheight; j++)
     {
      sumr=0,sumg=0,sumb=0;
      for(k=0; k<cwidth; k++) for(l=0; l<cheight; l++)
          {
           sumr=sumr+t1.getrpixel(i+k,j+l)*conv.getrpixel(k,l);
           sumg=sumg+t1.getgpixel(i+k,j+l)*conv.getgpixel(k,l);
           sumb=sumb+t1.getbpixel(i+k,j+l)*conv.getbpixel(k,l);
          }
      i2.putpixel(i,j,sumr/rweight,sumg/gweight,sumb/bweight);
     }
  return i2;
 }

Image I_mul(Image i1, Image mult)     //умножение
 {
  dword i,j, iheight = i1.getheight() <? mult.getheight(), iwidth = i1.getwidth() <? mult.getwidth();
  Image i2(i1);
  for (i=0; i<iwidth; i++) for (j=0; j<iheight; j++)
    i2.putpixel(i,j,i1.getrpixel(i,j)*mult.getrpixel(i,j)/255,i1.getgpixel(i,j)*mult.getgpixel(i,j)/255,i1.getbpixel(i,j)*mult.getbpixel(i,j)/255);
  return i2;
 }

Image I_sub(Image i1, Image sub)       //вычитание
 {
  dword i,j, iheight = i1.getheight() <? sub.getheight(), iwidth = i1.getwidth() <? sub.getwidth();
  Image i2(i1);
  for (i=0; i<iwidth; i++) for (j=0; j<iheight; j++)
    i2.putpixel(i,j,i1.getrpixel(i,j)-sub.getrpixel(i,j),i1.getgpixel(i,j)-sub.getgpixel(i,j),i1.getbpixel(i,j)-sub.getbpixel(i,j));
  return i2;
 }

Image I_dif(Image i1, Image diff)       //нормированная разница изображений (по модулю)
 {
  dword i,j, iheight = i1.getheight() <? diff.getheight(), iwidth = i1.getwidth() <? diff.getwidth();
  Image i2(i1);
  for (i=0; i<iwidth; i++)  for (j=0; j<iheight; j++)
     i2.putpixel(i,j,abs(i1.getrpixel(i,j)-diff.getrpixel(i,j)),abs(i1.getgpixel(i,j)-diff.getgpixel(i,j)),abs(i1.getbpixel(i,j)-diff.getbpixel(i,j)));
  return i2;
 }

Image I_norm(Image i1, dword bound)
 {
  dword iheight = i1.getheight(), iwidth = i1.getwidth(), i, j; double n=0;
  Image i2(i1);
  for(i=0;i<iwidth;i++) for(j=0;j<iheight;j++)                    //коэффициент нормирования
     { n= n >? i2.getrpixel(i,j); n=n >? i2.getgpixel(i,j); n=n >? i2.getbpixel(i,j); }
  n = n/bound;
  for (i=0; i<iwidth; i++)  for (j=0; j<iheight; j++)
     i2.putpixel(i,j,byte(i2.getrpixel(i,j)/n),byte(i2.getgpixel(i,j)/n),byte(i2.getbpixel(i,j)/n));
  return i2;
 }

Image hystogram(Image i1, dword mh)
 {
  dword i,j,n=0,ih=i1.getheight(),iw=i1.getwidth();                             //индексы, ширина, высота
  dword *_h;_h=new dword[768];for(i=0;i<768;i++)_h[i]=0;                        //массив счетчиков гистограммы
  for(i=0;i<iw;i++)for(j=0;j<ih;j++)                                            //по всем элементам изображения
   {
    _h[i1.getrpixel(i,j)*3]++;                                                  //красный
    _h[i1.getgpixel(i,j)*3+1]++;                                                //зеленый
    _h[i1.getbpixel(i,j)*3+2]++;                                                //синий
   }
  for(i=0;i<768;i++)n=n>?_h[i];                                                 //максимальное количество элементов
  if(n>mh){for(i=0;i<768;i++)_h[i]=_h[i]*mh/n;n=mh;}                            //если выходит за границы - нормируем
  Image i2(1024,mh);                                                            //изображение на выход
  for(i=0;i<1024;i++)for(j=0;j<mh;j++)i2.putpixel(i,j,192,192,192);             //заполняем серым тоном
  for(i=0;i<1024;i++)                                                           //по всем столбцам
   {for(j=0;j<n;j+=n/4)i2.putpixel(i,mh-j-1,0,0,0);i2.putpixel(i,mh-n,0,0,0);}  //горизонтальные полосы сетки
  for(i=0;i<256;i++)                                                            //по всем столбцам по 4 элемента
   {
    for(j=0;j<_h[i*3];j++)i2.putpixel(i*4,mh-j-1,255,0,0);                      //красный
    for(j=0;j<_h[i*3+1];j++)i2.putpixel(i*4+1,mh-j-1,0,255,0);                  //зеленый
    for(j=0;j<_h[i*3+2];j++)i2.putpixel(i*4+2,mh-j-1,0,0,255);                  //синий
    for(j=0;j<mh;j++)i2.putpixel(i*4+3,mh-j-1,0,0,0);                           //вертикальные линии сетки
   }
  delete [] _h;                                                                 //прибиваем массив
  return i2;                                                                    //возвращаем изображение
 }
