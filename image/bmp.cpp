/*  bmp.cpp                          библиотека  обработки  изображений v1.4. */
/*  содержит  функции  чтения  и  записи  изображений  в  файл  формата  bmp. */
/*  Автор   -   Vlad  Efremov.     v2@free.kursknet.ru     fido 2:5035/10.34. */
/*  exit   code   errors:    1 - open;  2 - read;  3 - wrong img;  4 - write. */

#include <stdio.h>                                                              //подключаем стандартную библиотеку ввода-вывода
#include "image.h"                                                              //подключаем к основному заголовку библиотеки
#include "defs.h"                                                               //подключаем файл с определениями


byte Image::loadimage(const char *fname)                                        //читаем изображение из файла - только bmp
 {
  BMPHEADER header; dword i;                                                    //переменная для заголовка дополнительная переменная
  FILE *f; if((f=fopen(fname,"rb"))==NULL)return 1;                             //файл, попытка открыть, если ошибка - выход
  if(fread(&header,sizeof(BMPHEADER),1,f)<1){fclose(f);return 2;}               //попытка прочитать заголовок, если ошибка - выход
  if(header.fh.Type!=0x4D42||header.ih.Compression){fclose(f);return 3;}        //wrong image type
  width=header.ih.Width; height=header.ih.Height;                               //ширина высота
  if(header.ih.BitCount==8){fclose(f);return 3;}                                //wrong image type
  dword len=width*3,size=height*len;                                            //расчет строки в байтах, размера
  mem=new byte[size];byte *memory=mem+size;                                     //выделяем память под изображение указатель на конец
  byte temp[4];dword plus=(4-len&3)&3;                                          //массив в 4 байта для выравнивания
  for(i=height;i>0;i--)                                                         //по количеству строк задом наперед поехали
   {
    memory-=len;                                                                //уменьшаем указатель конца на длину строки    
    if(fread(memory,len,1,f)<0)break;                                           //считываем в конец изображения одну строку из файла
    if(plus&&fread(temp,1,plus,f)<plus)break;                                   //считываем из файла выравнивающие байты
   }
  fclose(f);                                                                    //закрываем файл
  return (i>0)?2:0;                                                             //возвращаем read error, или 0, если все в порядке.
 }


byte Image::saveimage(const char *fname)const                                   //запись изображения в файл, BMP only
 {
  dword i,len=3*width,plus=(4-len&3)&3,size=(len+plus)*height;                  //счетчик, длина строки, выравнивание, размеры
  BMPHEADER header; FILE *f; if((f=fopen(fname,"wb"))==NULL)return 1;           //заголовок, файл, Попытка открыть файл, open error
  header.fh.Type = 0x4D42;                                                      // "BM" or 0x4D42
  header.fh.OffBits = sizeof(header.fh)+sizeof(header.ih);              // Offset in file where the bits begin
  header.fh.Reserved = 0;                                                       //зарезервировано
  header.fh.Size = header.fh.OffBits+size;                                      // Size of file in bytes
  header.ih.Planes = 1;                                                         //number of plenes
  header.ih.Compression = 0;                                                    //no compression only
  header.ih.BitCount = 24;                                                      //bits per pixel
  header.ih.SizeImage = size;                                                   //size of image
  header.ih.ClrUsed = 1L<<24;                                                   //Используемых цветов
  header.ih.ClrImportant = 0xFFFFFFFF;                                          //Важных цветов
  header.ih.Size = sizeof(header.ih);                                       //size of infoheader
  header.ih.Width = header.ih.XPelsPerMeter = width;                            //ширина
  header.ih.Height = header.ih.YPelsPerMeter = height;                          //высота
  if(fwrite(&header,sizeof(BMPHEADER),1,f)<1){fclose(f);return 4;}              //Попытка записать заголовок, write error

  char temp[4]={0,0,0,0};byte *memory=mem+(len*height);                         //4 байта временно, указатель на конец изображения
  for(i=height;i>0;i--)                                                         //записываем строки задом наперед
   {
    memory-=len;                                                                //указатель уменьшаем на строку
    if(fwrite(memory,len,1,f)<1)break;                                          //записываем 1 строку изображения из памяти в файл
    if(plus&&fwrite(temp,1,plus,f)<plus)break;                                  //записываем выравнивающие байты
   }
   fclose(f);                                                                   //закрыли файл
   return (i > 0) ? 4 : 0;                                                      //возвращаем write error, или 0, если все в порядке
 }
