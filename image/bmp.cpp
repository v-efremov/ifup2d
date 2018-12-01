/*  bmp.cpp                          ����������  ���������  ����������� v1.4. */
/*  ��������  �������  ������  �  ������  �����������  �  ����  �������  bmp. */
/*  �����   -   Vlad  Efremov.     v2@free.kursknet.ru     fido 2:5035/10.34. */
/*  exit   code   errors:    1 - open;  2 - read;  3 - wrong img;  4 - write. */

#include <stdio.h>                                                              //���������� ����������� ���������� �����-������
#include "image.h"                                                              //���������� � ��������� ��������� ����������
#include "defs.h"                                                               //���������� ���� � �������������


byte Image::loadimage(const char *fname)                                        //������ ����������� �� ����� - ������ bmp
 {
  BMPHEADER header; dword i;                                                    //���������� ��� ��������� �������������� ����������
  FILE *f; if((f=fopen(fname,"rb"))==NULL)return 1;                             //����, ������� �������, ���� ������ - �����
  if(fread(&header,sizeof(BMPHEADER),1,f)<1){fclose(f);return 2;}               //������� ��������� ���������, ���� ������ - �����
  if(header.fh.Type!=0x4D42||header.ih.Compression){fclose(f);return 3;}        //wrong image type
  width=header.ih.Width; height=header.ih.Height;                               //������ ������
  if(header.ih.BitCount==8){fclose(f);return 3;}                                //wrong image type
  dword len=width*3,size=height*len;                                            //������ ������ � ������, �������
  mem=new byte[size];byte *memory=mem+size;                                     //�������� ������ ��� ����������� ��������� �� �����
  byte temp[4];dword plus=(4-len&3)&3;                                          //������ � 4 ����� ��� ������������
  for(i=height;i>0;i--)                                                         //�� ���������� ����� ����� ������� �������
   {
    memory-=len;                                                                //��������� ��������� ����� �� ����� ������    
    if(fread(memory,len,1,f)<0)break;                                           //��������� � ����� ����������� ���� ������ �� �����
    if(plus&&fread(temp,1,plus,f)<plus)break;                                   //��������� �� ����� ������������� �����
   }
  fclose(f);                                                                    //��������� ����
  return (i>0)?2:0;                                                             //���������� read error, ��� 0, ���� ��� � �������.
 }


byte Image::saveimage(const char *fname)const                                   //������ ����������� � ����, BMP only
 {
  dword i,len=3*width,plus=(4-len&3)&3,size=(len+plus)*height;                  //�������, ����� ������, ������������, �������
  BMPHEADER header; FILE *f; if((f=fopen(fname,"wb"))==NULL)return 1;           //���������, ����, ������� ������� ����, open error
  header.fh.Type = 0x4D42;                                                      // "BM" or 0x4D42
  header.fh.OffBits = sizeof(header.fh)+sizeof(header.ih);              // Offset in file where the bits begin
  header.fh.Reserved = 0;                                                       //���������������
  header.fh.Size = header.fh.OffBits+size;                                      // Size of file in bytes
  header.ih.Planes = 1;                                                         //number of plenes
  header.ih.Compression = 0;                                                    //no compression only
  header.ih.BitCount = 24;                                                      //bits per pixel
  header.ih.SizeImage = size;                                                   //size of image
  header.ih.ClrUsed = 1L<<24;                                                   //������������ ������
  header.ih.ClrImportant = 0xFFFFFFFF;                                          //������ ������
  header.ih.Size = sizeof(header.ih);                                       //size of infoheader
  header.ih.Width = header.ih.XPelsPerMeter = width;                            //������
  header.ih.Height = header.ih.YPelsPerMeter = height;                          //������
  if(fwrite(&header,sizeof(BMPHEADER),1,f)<1){fclose(f);return 4;}              //������� �������� ���������, write error

  char temp[4]={0,0,0,0};byte *memory=mem+(len*height);                         //4 ����� ��������, ��������� �� ����� �����������
  for(i=height;i>0;i--)                                                         //���������� ������ ����� �������
   {
    memory-=len;                                                                //��������� ��������� �� ������
    if(fwrite(memory,len,1,f)<1)break;                                          //���������� 1 ������ ����������� �� ������ � ����
    if(plus&&fwrite(temp,1,plus,f)<plus)break;                                  //���������� ������������� �����
   }
   fclose(f);                                                                   //������� ����
   return (i > 0) ? 4 : 0;                                                      //���������� write error, ��� 0, ���� ��� � �������
 }
