/*  defs.h                           ����������  ���������  ����������� v1.4. */
/*  �������� ��������  �����  ��������� � �����������  ��������  �����������. */
/*  �����   -   Vlad  Efremov.     v2@free.kursknet.ru     fido 2:5035/10.34. */

#include <windows.h>

#ifndef _DEFS_H                      //���� �� ����������
#define _DEFS_H                      //����������


typedef unsigned short word;         //�����
typedef unsigned long dword;         //������� �����
typedef unsigned char byte;          //����
typedef byte BYTE;                   //�� �� ����� �������� �������.
typedef word WORD;                   //�� �� ����� �������� �������.
typedef dword DWORD;                 //�� �� ����� �������� �������.

#pragma pack(1)
typedef struct tagBITMAPFILEHEADER1   //���������, ��������� bmp �����
{
   WORD  Type;                       // "BM" or 0x4D42
   DWORD Size;                       // Size of file in bytes
   DWORD Reserved;             	     // Set to 0
   DWORD OffBits;                    // Offset in file where the bits begin
} BMPFILEHEADER1;

#pragma pack(1)
typedef struct tagBITMAPINFOHEADER1   //���������, �������������� ��������� �����
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
} BMPINFOHEADER1;


#pragma pack(1)
typedef struct tagBMPHEADER          //��������� - ��������� bmp
{
    struct tagBITMAPFILEHEADER1 fh;   //���������, ��������� bmp �����
    struct tagBITMAPINFOHEADER1 ih;   //���������, �������������� ��������� �����
} BMPHEADER;

#endif                               //���� ����������
