function [area1,area2,mix_area,prf1,prf2,mix_prf]=MixPrf(comp1,comp2,X)
%��1������comp1�����1�����ƣ�һ��Ϊ����
%��2������comp2�����2�����ƣ�һ��Ϊ�ܼ�
%��3������X����ֵ����������X(1)��Ӧ���1����ɣ�X(2)��Ӧ���2����ɡ���
[area1,prf1]=PurePrf(comp1);
[area2,prf2]=PurePrf(comp2);
mix_area=area1*X(1)+area2*X(2);%����������[�������������������������������]
mix_prf=prf1*X(1)+prf2*X(2);

