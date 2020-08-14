function [gamma]=seggamma(prf,A,T,parameterFlag)
% Calculates the activity coeffecient of segments
% gamma is the segments activiy coeffecient vector 
% prf is sigma-profile: [sigma non-hbPrf hbPrf]
% A is cosmo area of compund
% T is temperature,K
% parameterFlag is a flag to choose the interaction parameters for deltaW()

%% const
rgas=0.001987;          %���峣����kcal/mol/K
tol=0.000000001;            %������׼
%% sigma-profile informations
[numSig,prfLevel]=size(prf);
sigma=prf(:,1); % screening charge density
prfLevel=prfLevel-1;
numSig=prfLevel*numSig;
%% rebuild sigma-profile
% newPrf=zeros(numSig,1);
newPrf=[prf(:,2);prf(:,3)];
%% Now calculating the activity coeffecient of segments
gamma=ones(1,numSig);
converge=ones(1,numSig);
while (sum(converge)>=tol)
    gammaold=gamma;
    for m=1:numSig
        summation=0.0;
        for n=1:numSig
            if (newPrf(n)~=0)  % ��Ҫ�c0����ȱ��^���Ƿ�ԓ���]�Q��newPrf>eps
                summation=summation+newPrf(n)/A*gammaold(n)*exp(-deltaW(m,n,sigma,parameterFlag)/rgas/T);
            end
        end
        gamma(m)=exp(-log(summation));
        gamma(m)=(gamma(m)+gammaold(m))/2.0;  % �Ƿ��@���Z��ĳɣ�gamma=(gamma-gammaold)/2.0, ���������һ��ѭ����
    end
    converge=(gamma-gammaold).^2;
end

end
