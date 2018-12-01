/*  image.cpp                        ����������  ���������  ����������� v1.4. */
/*  �������� ������������, ���������, ������� ������� � �������� �����������. */
/*  �����   -   Vlad  Efremov.     v2@free.kursknet.ru     fido 2:5035/10.34. */


#include <memory.h>                                                             //���������� ���������� ��� ������ � �������
#include <iostream>                                                             //���������� ���������� ��� ������ � ��������
#include "image.h"                                                              //���������� � ��������� ��������� ����������
using namespace std;
Image :: Image(const char *fname) { loadimage(fname); }                         //����������� �� �����

Image :: Image(const Image &image)                                              //����������� �� �����������
 {
  width = image.width; height = image.height; dword size = width*height*3;      //���������� ������ ������ ��������
  mem = new byte [size]; memcpy(mem, image.mem, size);                          //�������� ������������ � �������� �����������
 }

Image& Image :: operator= (const Image &image)                                  //�������� ������������
 {
  if(mem) delete [] mem;                                                        //���� � �������� ������� ���-�� ���� ��������� ���
  width = image.width; height = image.height; dword size = width*height*3;      //���������� ������ ������ ��������
  mem = new byte [size]; memcpy(mem, image.mem, size);                          //�������� ������������ � �������� �����������
  return *this;                                                                 //���������� ��������� �� ���� ������
 }

void Image :: putpixel(dword x, dword y, byte r=0, byte g=0, byte b=0)          //������ ������
 { dword n=y*width*3+x*3; mem[n]=b; mem[++n]=g; mem[++n]=r; }                   //������ ������
