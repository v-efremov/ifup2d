/*  defs.h                           библиотека  обработки  изображений v1.4. */
/*  содержит описания  типов  элементов и применяемых  структур  изображения. */
/*  Автор   -   Vlad  Efremov.     v2@free.kursknet.ru     fido 2:5035/10.34. */


#ifndef _DEFS_H                      //если не определено
#define _DEFS_H                      //определяем


typedef unsigned short word;         //слово
typedef unsigned long dword;         //двойное слово
typedef unsigned char byte;          //байт
typedef byte BYTE;                   //то же самое большими буквами.
typedef word WORD;                   //то же самое большими буквами.
typedef dword DWORD;                 //то же самое большими буквами.

#pragma pack(1)
typedef struct tagBITMAPFILEHEADER   //структура, заголовок bmp файла
{
   WORD  Type;                       // "BM" or 0x4D42
   DWORD Size;                       // Size of file in bytes
   DWORD Reserved;             	     // Set to 0
   DWORD OffBits;                    // Offset in file where the bits begin
} BMPFILEHEADER;

#pragma pack(1)
typedef struct tagBITMAPINFOHEADER   //структура, информационный заголовок файла
{
   DWORD Size;                       // Size of the structure
   DWORD Width;                      // Width in pixels
   DWORD Height;                     // Height in pixels
   WORD  Planes;                     // # of color Planes: Set to 1
   WORD  BitCount;                   // Color bits per pixel
   DWORD Compression;                // Compression Scheme
   DWORD SizeImage;                  // Number of bitmap bytes
   DWORD XPelsPerMeter;              // Horizontal Resolution
   DWORD YPelsPerMeter;              // Vertical Resolution
   DWORD ClrUsed;                    // Number of colors used
   DWORD ClrImportant;               // Important colors
} BMPINFOHEADER;

#pragma pack(1)
typedef struct tagBMPHEADER          //структура - заголовок bmp
{
    struct tagBITMAPFILEHEADER fh;   //структура, заголовок bmp файла
    struct tagBITMAPINFOHEADER ih;   //структура, информационный заголовок файла
} BMPHEADER;

#endif                               //если определено
