/*  image.h                          ����������  ���������  ����������� v1.4. */
/*  ������� ��������� - �������� �����������  ������ ����������� �����������. */
/*  �����   -   Vlad  Efremov.     v2@free.kursknet.ru     fido 2:5035/10.34. */


#ifndef _IMAGE_H_                                                               //���� �� ����������
#define _IMAGE_H_                                                               //����������


#include "defs.h"                                                               //���������� ���� � �������������


class Image                                                                     //�������� ������ �����������.
 {
  private:                                                                      //���������� ����������
    dword width, height; byte *mem;                                             //������, ������, �����������.

  public:
    Image() : width(0), height(0), mem(0) {}                                    //����������� ��� ����������
    Image(const char *fname);                                                   //����������� �� �����
    Image(const Image &image);                                                  //����������� �� �����������
    Image(const dword w,const dword h):width(w),height(h),mem(new byte[w*h*3]){}//����������� �� �������� ����������.

    ~Image() { delete [] mem; }                                                 //����������

    Image& operator= (const Image &image);                                      //�������� ������������

    dword getwidth() const { return width; }                                    //���������� ������
    dword getheight() const { return height; }                                  //���������� ������

    byte getrpixel(dword x, dword y) { return mem[y*width*3+x*3+2]; }           //������� ������������.
    byte getgpixel(dword x, dword y) { return mem[y*width*3+x*3+1]; }           //������� ������������.
    byte getbpixel(dword x, dword y) { return mem[y*width*3+x*3]; }             //����� ������������.

    void putpixel(dword x, dword y, byte r, byte g, byte b);                    //������ ������.

    byte saveimage(const char *fname) const;                                    //���������� ����������� � ����

private:                                                                        //���������� �������
    byte loadimage(const char *fname);                                          //������ ����������� �� �����
 };

#endif                                                                          //���� ����������
