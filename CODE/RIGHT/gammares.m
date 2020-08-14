function [gamma,lngamma]=gammares(Vcosmo,Acosmo,X,prf1,prf2,T)
%第1个参数Vcosmo为组分的体积向量，Vcosmo(1)对应组分1的cosmo体积，Vcosmo(2)对应组分2的cosmo体积
%第2个参数Acosmo为组分的面积向量，Acosmo(1)对应组分1的cosmo面积，Acosmo(2)对应组分2的cosmo面积
%第3个参数X为组分的组成矩阵，每一行代表一个组成点；
%第4个参数prf1对应的是组分1的sigma-profile，为2维矩阵，column1对应横坐标，column2为非氢键的profile，column3为氢键部分的profile，单位为面积单位，A^2
%第5个参数prf2对应的是组分2的sigma-profile，为2维矩阵，column1对应横坐标，column2为非氢键的profile，column3为氢键部分的profile，单位为面积单位，A^2
%第6个参数mixprf对应的是组分1与组分2混合的sigma-profile，为2维矩阵，column1对应横坐标，column2为非氢键的profile，column3为氢键部分的profile，单位为面积单位，A^2
%第7个参数为体系的温度，单位为K
global Aeff;
%constant for calculation
fpol=0.6917;            %极化因子,from the paper of COSMO-SAC(2004)，COSMO-SAC(2002)的文献取0.64
%Aeff=7.25;              %有效标准片段面积，与片段平均化中数值保持一致,from the paper of COSMO-SAC(2007)
eps0=0.0002395;         %真空介电常数(e^2 mol)/((kcal)/(A) 
alpha=fpol*0.3*(Aeff^1.5)/2.0/eps0; %非氢片段相互作用能系数
chb=3484.42;            %氢键片段相互作用系数，kcal/(mol A^4)/e^2, from the paper of COSMO-SAC(2007)
% alpha=10072;             %来自杨犁博士论文
% chb=2797;               %来自杨犁博士论文
rgas=0.001987;          %气体常数，kcal/mol/K
tol=0.000001;            %收敛标准
numComp=length(Vcosmo); %组分数
numPoint=length(X(:,1));%数据点数
sumA=X*(Acosmo');       %混合后总面积,列向量，对应每个点的混合面积
sigma=prf1(:,1);        %sigma-profile横坐标数据
numSig=length(sigma);   %横坐标数据点数
prfLevel=2;             %profile拆分成的维度，这里取2表示profile拆分成非氢键和氢键部分（N O F以及与其相连的H）
%
%
%创建矩阵大小
deltaW=zeros(numSig,numSig,prfLevel,prfLevel);
%计算deltaW矩阵
for t=1:prfLevel
    for m=1:numSig
        for s=1:prfLevel
            for n=1:numSig
                if (s==2&t==2&(sigma(m)*sigma(n)<0.0))  %同时为氢键profile的情况
                    deltaW(m,n,t,s)=alpha*(sigma(m)+sigma(n))^2-chb*(sigma(m)-sigma(n))^2;
                else
                    deltaW(m,n,t,s)=alpha*(sigma(m)+sigma(n))^2;
                end
            end
        end
    end
end
%计算deltaW矩阵结束
%
    %计算纯组分片段的活度系数
%     if i==1
    convPure=ones(numSig,numComp,prfLevel);
    seggammaPureOld=ones(numSig,numComp,prfLevel);
    %+构造纯组分的sigma-profile矩阵
    purePrf=zeros(numSig,numComp,prfLevel);  %单位依然为面积单位，A^2
    purePrf(:,1,:)=prf1(:,2:3);
    purePrf(:,2,:)=prf2(:,2:3);
    for j=1:numComp
        seggammaPure=ones(numSig,numComp,prfLevel);
        iter=0;
        while (max(max(max(convPure(:,j,:))))>=tol)
            iter=iter+1;
            convPure(:,j,:)=zeros(numSig,prfLevel);    %有必要重置为0吗？
            seggammaPureOld(:,j,:)=seggammaPure(:,j,:);
            for m=1:numSig
                for t=1:prfLevel
                    summation=0.0;
                    for s=1:prfLevel
                        for n=1:numSig
                            summation=summation+(purePrf(n,j,s)/Acosmo(j))*seggammaPureOld(n,j,s)*exp(-deltaW(m,n,t,s)/rgas/T);  %要检查单位
                        end
                    end
                    seggammaPure(m,j,t)=exp(-log(summation));
                    seggammaPure(m,j,t)=(seggammaPure(m,j,t)+seggammaPureOld(m,j,t))/2.0;
                end
            end
            convPure(:,j,:)=abs(seggammaPure(:,j,:)-seggammaPureOld(:,j,:))./seggammaPureOld(:,j,:);
        end
    end
%     end
    %计算纯组分片段活度系数结束
%
gamma=ones(numPoint,numComp);
lngamma=zeros(numPoint,numComp);
sumgamma=zeros(numPoint,numComp);
for i=1:numPoint
    %计算混合profile
    %+分配混合profile大小
     mixprf=zeros(numSig,prfLevel);
    % mix_area=area1*X(1)+area2*X(2);%这里的面积：[总面积，非氢键部分面积，氢键部分面积]
    mixprf=prf1(:,2:3)*X(i,1)+prf2(:,2:3)*X(i,2);
    %计算混合片段活度系数
    seggamma=ones(numSig,prfLevel);
    converge=ones(numSig,prfLevel);
    while max(max(converge))>=tol
        seggammaold=seggamma;
        for t=1:prfLevel
            for m=1:numSig
                summation=0.0;   %要考虑放在这是否合适
                for s=1:prfLevel
                    for n=1:numSig
                        summation=summation+mixprf(n,s)/sumA(i)*seggammaold(n,s)*exp(-deltaW(m,n,t,s)/rgas/T);
                    end
                end
                seggamma(m,t)=exp(-log(summation));
                seggamma(m,t)=(seggamma(m,t)+seggammaold(m,t))/2.0;
            end
        end
        converge=abs((seggamma-seggammaold)./seggammaold);
    end
    %计算混合片段活度系数结束
%
    %计算活度系数剩余项
%     gamma=ones(1,numComp);
%     lngamma=ones(1,numComp);
%     sumgamma=ones(1,numComp);
    for k=1:numComp
        for b=1:prfLevel
            for j=1:numSig
                sumgamma(i,k)=sumgamma(i,k)+((purePrf(j,k,b)/Aeff)*(log(seggamma(j,b)/seggammaPure(j,k,b)))); %注意检查这个公式是否正确
            end
        end
       gamma(i,k)=exp(sumgamma(i,k));  %         gamma(i,k)=exp(sumgamma(i,k)+lngamma(i,k));检查这个公式的正确性
       lngamma(i,k)=log(gamma(i,k));
    end
end
end
    %活度系数剩余项计算结束
    
        
        