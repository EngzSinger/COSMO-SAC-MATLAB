function [gamma,lngamma] = GammaRes2007(solute,cation,anion,X,T,para,varargin)
% input arg should be first checked
% solute and il are must while some optional input
%tic;
%para = Para();
%toc
%tic;
Aeff = str2num(para('Aeff_2007'));
rgas = str2num(para('rgas'));
%toc
%tic;
p = inputParser;

% Required arg must be provided
p.addRequired('solute');
p.addRequired('cation');
p.addRequired('anion');
p.addRequired('X');
p.addRequired('T');
p.addRequired('para');


% parameter arg should be provide otherwise defult value would be utilized
defultSolutePath = '/home/zxwu/WORK_SPACE/COSMO_SAC/DATABASE/COSMO/SOLUTE/';
defultAnionPath = '/home/zxwu/WORK_SPACE/COSMO_SAC/DATABASE/COSMO/ANION/';
defultCationPath = '/home/zxwu/WORK_SPACE/COSMO_SAC/DATABASE/COSMO/CATION/';
defultILFlag = 1;
defultAlpha = str2num(para('Alpha'));
defultChb = str2num(para('Chb'));

p.addParameter('SolutePath',defultSolutePath);
p.addParameter('AnionPath',defultAnionPath);
p.addParameter('CationPath',defultCationPath);
p.addParameter('IL_Flag',defultILFlag);
p.addParameter('Alpha',defultAlpha);
p.addParameter('Chb',defultChb);

p.parse(solute,cation,anion,X,T,para,varargin{:});

cationName = p.Results.cation;
anionName = p.Results.anion;
soluteName = p.Results.solute;
cationPath = p.Results.CationPath;
anionPath = p.Results.AnionPath;
solutePath = p.Results.SolutePath;
alpha = p.Results.Alpha;
chb = p.Results.Chb;
%toc

% method of ionic liquid treatment
% IL_Flag = 1 ;ionic liquid as netral compound
% IL_Flag = 2 ;ionic liquid as totally departure
% IL_Flag = 3 ; ionic liquid as partly compound and partly ions
IL_Flag = p.Results.IL_Flag;

[Asolute,Vsolute,Psolute,nsolute] = Cosmo2Profile2007(solute,solutePath,para);
[Acation,Vcation,Pcation,ncation] = Cosmo2Profile2007(cation,cationPath,para);
[Aanion,Vanion,Panion,nanion] = Cosmo2Profile2007(anion,anionPath,para);
%il treatment prfC + prfA
AIL = Acation+Aanion;
VIL = Vcation+Vanion;
PIL = Pcation;
PIL(:,2:3) = PIL(:,2:3)+Panion(:,2:3);
%PIL
n = nsolute + ncation + nanion;
npure = [nsolute,ncation+nanion];

Acosmo = [Asolute(1),AIL(1)];
Vcosmo = [Vsolute(1),VIL(1)];

tol = 0.000000001;
numComp = length(Vcosmo);
%numPoint = length(X(:,1));
sumA = X*(Acosmo');
sigma = Psolute(:,1);
numsig = 2*length(sigma);

W = ones(numsig,numsig);
for i = 1:numsig
	for j = 1:numsig
		% check if the symmetry value is calculated,and save nearly half calculation
		if(i<=j)
			W(i,j) = deltaW(i,j,alpha,chb);
		else
			W(i,j) = W(j,i);
		end
	end
end

ConvPure = ones(numComp,numsig);
seggammaPure = ones(numComp,numsig);
PurePrf = zeros(numsig,numComp);
PurePrf(:,1) = [Psolute(:,2);Psolute(:,3)];
PurePrf(:,2) = [PIL(:,2);PIL(:,3)];
PurePrf
Acosmo
for i = 1:numComp
	istep = 0;
	while (sum(ConvPure(i,:)) >= tol)
		%istep = istep + 1
		ConvPure(i,:) = zeros(numsig,1);
		seggammaOld = seggammaPure(i,:);
		for j = 1:numsig
			sumation = 0.0;
			for m = 1:numsig
				if (PurePrf(m,i)~=0)
				sumation = sumation+(PurePrf(m,i)/Acosmo(i))*seggammaOld(m)*exp(-W(j,m)/rgas/T);
				end
			end
			seggammaPure(i,j) = (exp(-log(sumation))+seggammaOld(j))/2.0;
		end
		ConvPure(i,:) = (seggammaPure(i,:)-seggammaOld).^2;
	end
end

seggammaMix = ones(numComp,numsig);
ConvMix = ones(numComp,numsig);
Asum = sum(Acosmo);
gamma=ones(1,numComp);
lngamma=zeros(1,numComp);
sumgamma=zeros(1,numComp);


TempMixPrf = Psolute(:,2:3)*X(1)+PIL(:,2:3)*X(2);
MixPrf = [TempMixPrf(:,1);TempMixPrf(:,2)];
%MixPrf
%MixPrf = [TempMixPrf(:,2);TempMixPrf(:,3)];
for i = 1:numComp
	istep = 0;
       	while (sum(ConvMix(i,:)) >= tol)
		%istep = istep + 1
               	ConvMix(i,:) = zeros(numsig,1);
               	seggammaOld = seggammaMix(i,:);
               	for j = 1:numsig
                       	summation = 0.0;
                       	for m = 1:numsig
                               	if (MixPrf(m)~=0)
                               	sumation = sumation+(MixPrf(m)/sumA)*seggammaOld(m)*exp(-W(j,m)/rgas/T);
                               	end
                       	end
                       	seggammaMix(i,j) = (exp(-log(sumation))+seggammaOld(j))/2.0;
               	end
               	ConvMix(i,:) = (seggammaMix(i,:)-seggammaOld).^2;
       	end
end
for i = 1:numComp
	for j = 1:numsig
		%seggammaMix(i,j)
		%seggammaPure(i,j)
		ab = (seggammaMix(i,j)/seggammaPure(i,j))
		MixPrf(j)
		MixPrf(j)/Aeff
		sumgamma(1,i) = sumgamma(1,i)+((MixPrf(j)/Aeff*(log(seggammaMix(i,j)/seggammaPure(i,j)))));
		%sumgamma(i) = sumgamma(i)+((PurePrf(j,i)/Aeff*(log(seggammaMix(i,j)/seggammaPure(i,j)))));
		%sumgamma
		%sumgamma(l,i) = sumgamma(l,i)+(n*(MixPrf(j)/Asum*(log(seggammaMix(i,j)))))-(npure(i)*(PurePrf(j,i)/Acosmo(i)*(log(seggammaPure(i,j)))));
	end
	gamma(1,i) = exp(sumgamma(1,i));
	lngamma(1,i) = log(gamma(1,i));
end
end
