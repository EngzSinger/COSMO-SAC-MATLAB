function [gamma]=gammaOW(comp3,x)
% This function aims to calculate the Octanol-Water partition coefficient
% of comp3 at 298.15 K
% x(1,1),x(1,2),x(1,3) are the mole fraction of octano,water,comp3 in octanol-reach phase,respectively
% x(2,1),x(2,2),x(2,3) are the mole fraction of octano,water,comp3 in
% octanol-reach phase,respectively
%+选择计算组合项活度系数的公式 combFlag=1;
%+体系类型选择，ILflag（=1：离子液体体系；=0：非离子液体体系）ILflag=1;
%+计算交互作用的参数选择,parameterFlag(=1：COSMO-SAC(2007)中的参数；=2：Wang Shu C++程序中的参数；=3：杨犁博士论文的参数) parameterFlag=2;
%+选择cosmo文件来源cosmoFlag(=1：COSMOTherm软件；=2：VT-2005；=3：G09生成) cosmoFlag=2;
global combFlag parameterFlag cosmoFlag Aeff;
%constant for calculation
%Aeff=7.25;              %有效标准片段面积，与片段平均化中数值保持一致,from the paper of COSMO-SAC(2007)
rgas=0.001987;          %气体常数，kcal/mol/K
tol=0.000000001;            %收敛标准
numComp=3;      %组分数
numPoint=2;     %数据点
T=298.15;
comp1='octanol';
comp2='water';
[area1,volume1,prf1]=PurePrf(comp1,cosmoFlag);
[area2,volume2,prf2]=PurePrf(comp2,cosmoFlag);
[area3,volume3,prf3]=PurePrf(comp3,cosmoFlag);
Acosmo=[area1(1),area2(1),area3(1)];
Vcosmo=[volume1,volume2,volume3];
sumA=x*(Acosmo');       %混合后总面积,列向量，对应每个点的混合面积
sigma=prf1(:,1);        %sigma-profile横坐标数据
prfLevel=2;             %profile拆分成的维度，这里取2表示profile拆分成非氢键和氢键部分（N O F以及与其相连的H）
numSig=prfLevel*length(sigma);   %新profile横坐标数据点数
%
%
%计算纯组分片段的活度系数
%     if i==1
    convPure=ones(numComp,numSig);
    seggammaPureOld=ones(numComp,numSig);
    %+构造纯组分的sigma-profile矩阵
    purePrf=zeros(numSig,numComp);  %单位依然为面积单位，A^2
    purePrf(:,1)=[prf1(:,2);prf1(:,3)];
    purePrf(:,2)=[prf2(:,2);prf2(:,3)];
    purePrf(:,3)=[prf3(:,2);prf3(:,3)];
    seggammaPure=ones(numComp,numSig);
    for j=1:numComp        
%         seggammaPure=ones(numComp,numSig);
        istep=0;  %记录迭代次数
        while (sum(convPure(j,:))>=tol)
            istep=istep+1;
            convPure(j,:)=zeros(numSig,1);    %有必要重置为0吗？
            seggammaPureOld(j,:)=seggammaPure(j,:);
            for m=1:numSig
                summation=0.0;
                    for n=1:numSig
                        if (purePrf(n,j)~=0)
                            summation=summation+(purePrf(n,j)/Acosmo(j))*seggammaPureOld(j,n)*exp(-deltaW(m,n,sigma,parameterFlag)/rgas/T);  %要检查单位
                        end                            
                    end
                    seggammaPure(j,m)=exp(-log(summation));
                    seggammaPure(j,m)=(seggammaPure(j,m)+seggammaPureOld(j,m))/2.0;
            end
             convPure(j,:)=(seggammaPure(j,:)-seggammaPureOld(j,:)).^2;  %by Tu 2016-04-28
         end
            %convPure(j,:)=abs(seggammaPure(j,:)-seggammaPureOld(j,:))./seggammaPureOld(j,:);
           %convPure(j,:)=(seggammaPure(j,:)-seggammaPureOld(j,:)).^2;  %by Tu 2016-04-28
     end
    
%     end
    %计算纯组分片段活度系数结束
%
gammaR=ones(numPoint,numComp);
lngammaR=zeros(numPoint,numComp);
sumgamma=zeros(numPoint,numComp);
for i=1:numPoint
    %计算混合profile
    %+分配混合profile大小
     TempMixprf=zeros(numSig,1);
    % mix_area=area1*X(1)+area2*X(2);%这里的面积：[总面积，非氢键部分面积，氢键部分面积]
    mixprf=prf1(:,2:3)*x(i,1)+prf2(:,2:3)*x(i,2)+prf3(:,2:3)*x(i,3);
    %构造临时混合prf
    TempMixprf=[mixprf(:,1);mixprf(:,2)];
    %计算混合片段活度系数
    seggamma=ones(1,numSig);
    converge=ones(1,numSig);
    while sum(converge)>=tol
        seggammaold=seggamma;        
            for m=1:numSig
                summation=0.0;   %要考虑放在这是否合适                
                    for n=1:numSig
                        if (mixprf(n)~=0)
                            summation=summation+mixprf(n)/sumA(i)*seggammaold(n)*exp(-deltaW(m,n,sigma,parameterFlag)/rgas/T);
                        end
                    end                
                seggamma(m)=exp(-log(summation));
                seggamma(m)=(seggamma(m)+seggammaold(m))/2.0;
            end        
        converge=(seggamma-seggammaold).^2;
    end
    %计算混合片段活度系数结束
%
    %计算活度系数剩余项
%     gamma=ones(1,numComp);
%     lngamma=ones(1,numComp);
%     sumgamma=ones(1,numComp);
    for k=1:numComp        
            for j=1:numSig
                sumgamma(i,k)=sumgamma(i,k)+((purePrf(j,k)/Aeff)*(log(seggamma(j)/seggammaPure(k,j)))); %注意检查这个公式是否正确
            end        
       gammaR(i,k)=exp(sumgamma(i,k));  %         gamma(i,k)=exp(sumgamma(i,k)+lngamma(i,k));检查这个公式的正确性
       lngammaR(i,k)=log(gammaR(i,k));
    end
end
    %活度系数剩余项计算结束
%     
% ++开始计算组合项活度系数
gammaC=ones(numPoint,numComp);
lngammaC=zeros(numPoint,numComp);
for i=1:numPoint
[gammaC(i,:),lngammaC(i,:)]=gammacomb(Vcosmo,Acosmo,x(i,:),combFlag);
end

% 计算活度系数
lngamma=lngammaR+lngammaC;
gamma=exp(lngamma);

end