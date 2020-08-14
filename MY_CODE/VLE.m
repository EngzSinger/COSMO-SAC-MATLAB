function P = VLE(solute,cation,anion,x,T,para,varargin)
%calculate vapor liquid equlibrium data 

defultAlpha = para('Alpha');
defultChb = para('Chb');

p = inputParser;

p.addRequired('solute');
p.addRequired('cation');
p.addRequired('anion');
p.addRequired('x');
p.addRequired('T');
p.addRequired('para');

p.addParameter('Alpha',defultAlpha);
p.addParameter('Chb',defultChb);

p.parse(solute,cation,anion,x,T,para,varargin{:});

alpha = p.Results.Alpha;
chb = p.Results.Chb;

%saturation vapor pressure calculation
Psat = @(T) 10^(5.20409-(1581.341/(T-33.50)))*100;



[gammaC,lngammaC] = GammaComb(solute,cation,anion,x,para);
[gammaR,lngammaR] = GammaRes2007(solute,cation,anion,x,T,para,'Alpha',alpha,'Chb',chb);
size(lngammaR)
size(lngammaC)
lngammaC
lngammaR
lngamma = lngammaC + lngammaR;
gamma = exp(lngamma);
activity = x(:,1).*gamma(:,1);
activity(activity>1) = 1;
P = activity.*Psat(T);
