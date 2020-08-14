function W = deltaW(i,j,Alpha,Chb)
%tic;
%para = Para();
%defultAlpha = str2num(para('Alpha'));
%defultChb = str2num(para('Chb'));
%toc
%tic;
%p = inputParser;

%p.addRequired('i');
%p.addRequired('j');
%p.addParameter('Alpha',defultAlpha);
%p.addParameter('Chb',defultChb);

%p.parse(i,j,varargin{:});
%toc
%Alpha = p.Results.Alpha;
%Chb = p.Results.Chb;
sigRange = 0.035;
interval = 0.001;
pointNum = 71;
sigma = -sigRange:interval:sigRange;
sigma = sigma';

countSigma = length(sigma);
Ehb = 0.0;

sigmaNew = [sigma;sigma];
Emf = Alpha*(sigmaNew(i)+sigmaNew(j))^2;
if (i>countSigma & j>countSigma & sigmaNew(i)*sigmaNew(j)<0.0)
	Ehb = Chb*(sigmaNew(i)-sigmaNew(j))^2;
end
W = Emf-Ehb;
end
