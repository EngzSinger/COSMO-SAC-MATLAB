function [gamma]=seggamma(prf,A,T,parameterFlag)
% Calculates the activity coeffecient of segments
% gamma is the segments activiy coeffecient vector 
% prf is sigma-profile: [sigma non-hbPrf hbPrf]
% A is cosmo area of compund
% T is temperature,K
% parameterFlag is a flag to choose the interaction parameters for deltaW()

%% const
rgas=0.001987;          %气体常数，kcal/mol/K
tol=0.000000001;            %收敛标准
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
            if (newPrf(n)~=0)  % 不要與0用相等比較，是否應該考慮換成newPrf>eps
                summation=summation+newPrf(n)/A*gammaold(n)*exp(-deltaW(m,n,sigma,parameterFlag)/rgas/T);
            end
        end
        gamma(m)=exp(-log(summation));
        gamma(m)=(gamma(m)+gammaold(m))/2.0;  % 是否將這個語句改成：gamma=(gamma-gammaold)/2.0, 將其放在下一层循环外
    end
    converge=(gamma-gammaold).^2;
end

end
