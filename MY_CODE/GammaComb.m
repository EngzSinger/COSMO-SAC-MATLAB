function [gamma,lngamma] = GammaComb(solute,cation,anion,x,para,varargin)

%para = Para();
q = str2num(para('anorm'));
r = str2num(para('vnorm'));
z = str2num(para('z'));
defultCationPath = para('CationPath');
defultAnionPath = para('AnionPath');
defultSolutePath = para('SolutePath');

p = inputParser;
p.addRequired('solute');
p.addRequired('cation');
p.addRequired('anion');
p.addRequired('x');
p.addOptional('CationPath',defultCationPath);
p.addOptional('AnionPath',defultAnionPath);
p.addOptional('SolutePath',defultSolutePath);

p.parse(solute,cation,anion,x,varargin{:});

CationPath = p.Results.CationPath;
AnionPath = p.Results.AnionPath;
SolutePath = p.Results.SolutePath;

%q = 79.53;
% area normalization constant A^2
%r = 66.69;
% volume normalization constant A^3
%z = 10.0;
% cooridination number 

[Asolute,Vsolute,Psolute] = Cosmo2Profile2007(solute,SolutePath,para);
[Acation,Vcation,Pcation] = Cosmo2Profile2007(cation,CationPath,para);
[Aanion,Vanion,Panion] = Cosmo2Profile2007(anion,AnionPath,para);

%il treatment prfC + prfA
AIL = Acation+Aanion;
VIL = Vcation+Vanion;
%PIL = Pcation+Panion;

Acosmo = [Asolute(1),AIL(1)];
Vcosmo = [Vsolute,VIL];
%Acosmo
%Vcosmo

%numPoint = length(x(:,1));
compNum = max(size(Vcosmo));
gamma = ones(compNum);
lngamma = zeros(compNum);

qi = Acosmo/q;
%qi
ri = Vcosmo/r;
%ri
phi = zeros(compNum);
theta = zeros(compNum);
l = (z/2.0)*(ri-qi)-(ri-1.0);
%l
	theta = x.*qi/sum(x.*qi);
	phi = x.*ri/sum(x.*ri);
%	theta
%	phi
	lngamma = log(phi./x)+z/2.0*qi.*log(theta./phi)+l-phi./x*sum(x.*l);
gamma=exp(lngamma);

end
