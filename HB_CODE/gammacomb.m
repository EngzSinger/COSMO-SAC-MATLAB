function [gamma,lngamma]=gammacomb(Vcosmo,Acosmo,x,flag)
%第一个参数Vcosmo为各组分的COSMO体积向量，默认为行向量
%第二个参数Acosmo为各组分的COSMO面积向量，默认为行向量
%第三个参数x为各组分的组成向量，默认为行向量
%第四个参数flag为计算组合项活度系数的选择，flag=0：选用STAVERMAN-GUGGENHEIM方程，flag=1：选用UNIFAC方程
%
%canstants for STAVERMAN-GUGGENHEIM equation
vnorm=66.69;  %volume normalization constant (A^3)
anorm=79.63;  %area normalization constant (A^2)
z=10.0;         %the cooridination number
%canstants for UNIFAC equation
% stdr=25.19933555; %standard volume segment
% stdq=41.5282392;  %standard surface area segment
stdr=66.693927000000002; %standard volume segment,来自Wang Shu C++程序
stdq=79.531953999999999;  %standard surface area segment,来自Wang Shu C++程序
numPoint=length(x(:,1)); %number of data points
comp=max(size(Vcosmo));   %the numbers of compunds
lngamma=zeros(numPoint,comp);
gamma=ones(numPoint,comp);
%calculate the combinatorial activity of each data point

Vcosmo
Acosmo
%calculate the combinatorial activity coefficients by STAVERMAN-GUGGENHEIM equation
if flag==0
    
    %+calculate the normal volume and area for STAVERMAN-GUGGENHEIM equation
    %+allocate space
    rnorm=zeros(1,comp);      %allocate the space of rnorm
    qnorm=zeros(1,comp);     %allocate the space of qnorm
    theta=zeros(numPoint,comp);
    phi=zeros(numPoint,comp);
    l=zeros(1,comp);
    %++calculate
    rnorm=Vcosmo/vnorm;       %calculate the rnorm
    qnorm=Acosmo/anorm;       %calculate the qnorm
rnorm
qnorm
    l=(z/2.0)*(rnorm-qnorm)-(rnorm-1.0);   %与文献不一致
    l
    for i=1:numPoint
    theta(i,:)=x(i,:).*qnorm/sum(x(i,:).*qnorm);
    phi(i,:)=x(i,:).*rnorm/sum(x(i,:).*rnorm);
    lngamma(i,:)=log(phi(i,:)./x(i,:))+z/2.0*qnorm.*log(theta(i,:)./phi(i,:))+l-(phi(i,:)./x(i,:))*sum(x(i,:).*l);
    end
%calculate the combinatorial activity coefficients by UNIFAC equation    
else if flag==1
        V=zeros(1,comp);
        F=zeros(1,comp);
        q=zeros(1,comp);
        V=Vcosmo/stdr;
        F=Acosmo/stdq;
        q=F;
        for i=1:numPoint
        sumV=sum(x(i,:).*V);
        sumF=sum(x(i,:).*F);
        V=V/sumV;
        F=F/sumF;
        lngamma(i,:)=1-V+log(V)-5*q.*(1-V./F+log(V./F));
        end
    end
end
%calculate the combinatorial activity coefficients

gamma=exp(lngamma);
end


    
