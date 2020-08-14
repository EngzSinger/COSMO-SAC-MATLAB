function [gamma,lngamma]=gammares(Vcosmo,Acosmo,X,prf1,prf2,T)
%��1������VcosmoΪ��ֵ����������Vcosmo(1)��Ӧ���1��cosmo�����Vcosmo(2)��Ӧ���2��cosmo���
%��2������AcosmoΪ��ֵ����������Acosmo(1)��Ӧ���1��cosmo�����Acosmo(2)��Ӧ���2��cosmo���
%��3������XΪ��ֵ���ɾ���ÿһ�д���һ����ɵ㣻
%��4������prf1��Ӧ�������1��sigma-profile��Ϊ2ά����column1��Ӧ�����꣬column2Ϊ�������profile��column3Ϊ������ֵ�profile����λΪ�����λ��A^2
%��5������prf2��Ӧ�������2��sigma-profile��Ϊ2ά����column1��Ӧ�����꣬column2Ϊ�������profile��column3Ϊ������ֵ�profile����λΪ�����λ��A^2
%��6������mixprf��Ӧ�������1�����2��ϵ�sigma-profile��Ϊ2ά����column1��Ӧ�����꣬column2Ϊ�������profile��column3Ϊ������ֵ�profile����λΪ�����λ��A^2
%��7������Ϊ��ϵ���¶ȣ���λΪK
global Aeff;
%constant for calculation
fpol=0.6917;            %��������,from the paper of COSMO-SAC(2004)��COSMO-SAC(2002)������ȡ0.64
%Aeff=7.25;              %��Ч��׼Ƭ���������Ƭ��ƽ��������ֵ����һ��,from the paper of COSMO-SAC(2007)
eps0=0.0002395;         %��ս�糣��(e^2 mol)/((kcal)/(A) 
alpha=fpol*0.3*(Aeff^1.5)/2.0/eps0; %����Ƭ���໥������ϵ��
chb=3484.42;            %���Ƭ���໥����ϵ����kcal/(mol A^4)/e^2, from the paper of COSMO-SAC(2007)
% alpha=10072;             %�������粩ʿ����
% chb=2797;               %�������粩ʿ����
rgas=0.001987;          %���峣����kcal/mol/K
tol=0.000001;            %������׼
numComp=length(Vcosmo); %�����
numPoint=length(X(:,1));%���ݵ���
sumA=X*(Acosmo');       %��Ϻ������,����������Ӧÿ����Ļ�����
sigma=prf1(:,1);        %sigma-profile����������
numSig=length(sigma);   %���������ݵ���
prfLevel=2;             %profile��ֳɵ�ά�ȣ�����ȡ2��ʾprofile��ֳɷ������������֣�N O F�Լ�����������H��
%
%
%���������С
deltaW=zeros(numSig,numSig,prfLevel,prfLevel);
%����deltaW����
for t=1:prfLevel
    for m=1:numSig
        for s=1:prfLevel
            for n=1:numSig
                if (s==2&t==2&(sigma(m)*sigma(n)<0.0))  %ͬʱΪ���profile�����
                    deltaW(m,n,t,s)=alpha*(sigma(m)+sigma(n))^2-chb*(sigma(m)-sigma(n))^2;
                else
                    deltaW(m,n,t,s)=alpha*(sigma(m)+sigma(n))^2;
                end
            end
        end
    end
end
%����deltaW�������
%
    %���㴿���Ƭ�εĻ��ϵ��
%     if i==1
    convPure=ones(numSig,numComp,prfLevel);
    seggammaPureOld=ones(numSig,numComp,prfLevel);
    %+���촿��ֵ�sigma-profile����
    purePrf=zeros(numSig,numComp,prfLevel);  %��λ��ȻΪ�����λ��A^2
    purePrf(:,1,:)=prf1(:,2:3);
    purePrf(:,2,:)=prf2(:,2:3);
    for j=1:numComp
        seggammaPure=ones(numSig,numComp,prfLevel);
        iter=0;
        while (max(max(max(convPure(:,j,:))))>=tol)
            iter=iter+1;
            convPure(:,j,:)=zeros(numSig,prfLevel);    %�б�Ҫ����Ϊ0��
            seggammaPureOld(:,j,:)=seggammaPure(:,j,:);
            for m=1:numSig
                for t=1:prfLevel
                    summation=0.0;
                    for s=1:prfLevel
                        for n=1:numSig
                            summation=summation+(purePrf(n,j,s)/Acosmo(j))*seggammaPureOld(n,j,s)*exp(-deltaW(m,n,t,s)/rgas/T);  %Ҫ��鵥λ
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
    %���㴿���Ƭ�λ��ϵ������
%
gamma=ones(numPoint,numComp);
lngamma=zeros(numPoint,numComp);
sumgamma=zeros(numPoint,numComp);
for i=1:numPoint
    %������profile
    %+������profile��С
     mixprf=zeros(numSig,prfLevel);
    % mix_area=area1*X(1)+area2*X(2);%����������[�������������������������������]
    mixprf=prf1(:,2:3)*X(i,1)+prf2(:,2:3)*X(i,2);
    %������Ƭ�λ��ϵ��
    seggamma=ones(numSig,prfLevel);
    converge=ones(numSig,prfLevel);
    while max(max(converge))>=tol
        seggammaold=seggamma;
        for t=1:prfLevel
            for m=1:numSig
                summation=0.0;   %Ҫ���Ƿ������Ƿ����
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
    %������Ƭ�λ��ϵ������
%
    %������ϵ��ʣ����
%     gamma=ones(1,numComp);
%     lngamma=ones(1,numComp);
%     sumgamma=ones(1,numComp);
    for k=1:numComp
        for b=1:prfLevel
            for j=1:numSig
                sumgamma(i,k)=sumgamma(i,k)+((purePrf(j,k,b)/Aeff)*(log(seggamma(j,b)/seggammaPure(j,k,b)))); %ע���������ʽ�Ƿ���ȷ
            end
        end
       gamma(i,k)=exp(sumgamma(i,k));  %         gamma(i,k)=exp(sumgamma(i,k)+lngamma(i,k));��������ʽ����ȷ��
       lngamma(i,k)=log(gamma(i,k));
    end
end
end
    %���ϵ��ʣ����������
    
        
        