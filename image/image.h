/*  image.h                          библиотека  обработки  изображений v1.4. */
/*  базовый заголовок - содержит определение  класса изображения изображения. */
/*  Автор   -   Vlad  Efremov.     v2@free.kursknet.ru     fido 2:5035/10.34. */


#ifndef _IMAGE_H_                                                               //если не определено
#define _IMAGE_H_                                                               //определяем


#include "defs.h"                                                               //подключаем файл с определениями


class Image                                                                     //описание класса изображения.
 {
  private:                                                                      //внутренние переменные
    dword width, height; byte *mem;                                             //ширина, высота, изображение.

  public:
    Image() : width(0), height(0), mem(0) {}                                    //конструктор без параметров
    Image(const char *fname);                                                   //конструктор от файла
    Image(const Image &image);                                                  //конструктор от изображения
    Image(const dword w,const dword h):width(w),height(h),mem(new byte[w*h*3]){}//конструктор по заданным параметрам.

    ~Image() { delete [] mem; }                                                 //деструктор

    Image& operator= (const Image &image);                                      //оператор присваивания

    dword getwidth() const { return width; }                                    //возвращает ширину
    dword getheight() const { return height; }                                  //возвращает высоту

    byte getrpixel(dword x, dword y) { return mem[y*width*3+x*3+2]; }           //красная составляющая.
    byte getgpixel(dword x, dword y) { return mem[y*width*3+x*3+1]; }           //зеленая составляющая.
    byte getbpixel(dword x, dword y) { return mem[y*width*3+x*3]; }             //синяя составляющая.

    void putpixel(dword x, dword y, byte r, byte g, byte b);                    //кладем пиксел.

    byte saveimage(const char *fname) const;                                    //записываем изображение в файл

private:                                                                        //внутренние функции
    byte loadimage(const char *fname);                                          //читаем изображение из файла
 };

#endif                                                                          //если определено
