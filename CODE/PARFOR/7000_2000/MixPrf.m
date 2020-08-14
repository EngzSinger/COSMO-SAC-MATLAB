function [area1,area2,mix_area,prf1,prf2,mix_prf]=MixPrf(comp1,comp2,X)
%第1个参数comp1是组分1的名称，一般为溶质
%第2个参数comp2是组分2的名称，一般为溶剂
%第3个参数X是组分的组成向量，X(1)对应组分1的组成，X(2)对应组分2的组成……
[area1,prf1]=PurePrf(comp1);
[area2,prf2]=PurePrf(comp2);
mix_area=area1*X(1)+area2*X(2);%这里的面积：[总面积，非氢键部分面积，氢键部分面积]
mix_prf=prf1*X(1)+prf2*X(2);

