#ifndef _IMAGEOP_H_                                                               //если не определено
#define _IMAGEOP_H_                                                               //определяем

Image grayscale(Image i1);                         //сделать черно-белым.
Image I_conv(Image i1, Image conv);         //дискретная свертка
Image I_mul(Image i1, Image mult);        //умножение
Image I_sub(Image i1, Image sub);           //вычитание
Image I_dif(Image i1, Image diff);             //разница изображений (по модулю)
Image I_norm(Image i1, dword bound=255);        //нормирование
Image hystogram(Image i1, dword mh=768);        //гистограмма

#endif                                                                          //если определено
