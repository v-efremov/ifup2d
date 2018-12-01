/*  defs.h                           ������⥪�  ��ࠡ�⪨  ����ࠦ���� v1.4. */
/*  ᮤ�ন� ���ᠭ��  ⨯��  ����⮢ � �ਬ��塞��  �������  ����ࠦ����. */
/*  ����   -   Vlad  Efremov.     v2@free.kursknet.ru     fido 2:5035/10.34. */


#ifndef _DEFS_H                      //�᫨ �� ��।�����
#define _DEFS_H                      //��।��塞


typedef unsigned short word;         //᫮��
typedef unsigned long dword;         //������� ᫮��
typedef unsigned char byte;          //����
typedef byte BYTE;                   //� �� ᠬ�� ����訬� �㪢���.
typedef word WORD;                   //� �� ᠬ�� ����訬� �㪢���.
typedef dword DWORD;                 //� �� ᠬ�� ����訬� �㪢���.

#pragma pack(1)
typedef struct tagBITMAPFILEHEADER   //�������, ��������� bmp 䠩��
{
   WORD  Type;                       // "BM" or 0x4D42
   DWORD Size;                       // Size of file in bytes
   DWORD Reserved;             	     // Set to 0
   DWORD OffBits;                    // Offset in file where the bits begin
} BMPFILEHEADER;

#pragma pack(1)
typedef struct tagBITMAPINFOHEADER   //�������, ���ଠ樮��� ��������� 䠩��
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
typedef struct tagBMPHEADER          //������� - ��������� bmp
{
    struct tagBITMAPFILEHEADER fh;   //�������, ��������� bmp 䠩��
    struct tagBITMAPINFOHEADER ih;   //�������, ���ଠ樮��� ��������� 䠩��
} BMPHEADER;

#endif                               //�᫨ ��।�����
