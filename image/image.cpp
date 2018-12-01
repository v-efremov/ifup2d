/*  image.cpp                        библиотека  обработки  изображений v1.4. */
/*  содержит конструкторы, операторы, функции доступа к элементу изображения. */
/*  Автор   -   Vlad  Efremov.     v2@free.kursknet.ru     fido 2:5035/10.34. */


#include <memory.h>                                                             //подключаем библиотеку для работы с памятью
#include <iostream>                                                             //подключаем библиотеку для работы с потоками
#include "image.h"                                                              //подключаем к основному заголовку библиотеки
using namespace std;
Image :: Image(const char *fname) { loadimage(fname); }                         //конструктор от файла

Image :: Image(const Image &image)                                              //изображение от изображения
 {
  width = image.width; height = image.height; dword size = width*height*3;      //определяем ширину высоту габариты
  mem = new byte [size]; memcpy(mem, image.mem, size);                          //пытаемся организовать и копируем изображение
 }

Image& Image :: operator= (const Image &image)                                  //оператор присваивания
 {
  if(mem) delete [] mem;                                                        //если в выходном объекте что-то есть прибиваем его
  width = image.width; height = image.height; dword size = width*height*3;      //определяем ширину высоту габариты
  mem = new byte [size]; memcpy(mem, image.mem, size);                          //пытаемся организовать и копируем изображение
  return *this;                                                                 //возвращаем указатель на этот объект
 }

void Image :: putpixel(dword x, dword y, byte r=0, byte g=0, byte b=0)          //положи пиксел
 { dword n=y*width*3+x*3; mem[n]=b; mem[++n]=g; mem[++n]=r; }                   //кладем пиксел
