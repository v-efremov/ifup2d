#ifndef _IMAGEOP_H_                                                               //���� �� ����������
#define _IMAGEOP_H_                                                               //����������

Image grayscale(Image i1);                         //������� �����-�����.
Image I_conv(Image i1, Image conv);         //���������� �������
Image I_mul(Image i1, Image mult);        //���������
Image I_sub(Image i1, Image sub);           //���������
Image I_dif(Image i1, Image diff);             //������� ����������� (�� ������)
Image I_norm(Image i1, dword bound=255);        //������������
Image hystogram(Image i1, dword mh=768);        //�����������

#endif                                                                          //���� ����������
