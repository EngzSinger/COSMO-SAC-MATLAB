function [P,TLimit,Parameters]=antoine(SoluteName,T,varargin)
% This function calculates the Saturation Vapor Pressure of compounds 
% at a specific Temperature using Antoine Equation:
% log10(P)=A-B/(T+C), where A,B,C are parameters of Antoine Equation
% Parameters is a input vector consists of A,B,C as the format:[A B C]
% TLimit is a vector as the format:[LowerTemperature UpperTemperature],
% which defines the suitable temperature range for Antoine Equation
% T is the temperature that the user want to calculate
% varargin is a variable argument consists of a experimental saturation
% vapor pressure be used to calculate VLE when the temperature T out of the
% range TLimit
% P is the calculating saturation vapor pressure be used to caculate VLE
global SoluteIndex;
nVarargs=length(varargin);
SoluteFlag=-1;
P=zeros(length(T),1);
for i=1:length(SoluteIndex{1,1})
    if strcmp(SoluteIndex{1,2}{i},SoluteName)
       SoluteFlag=i;
        break;
    end
end
if (SoluteFlag<0)
    disp(['ANTOINE EQUATION MESSAGES:CAN''T NOT FIND',SoluteName,'IN SouteIndex FILE!']);
end
LowTemp=SoluteIndex{1,3}(SoluteFlag);
UpTemp=SoluteIndex{1,4}(SoluteFlag);
A=SoluteIndex{1,5}(SoluteFlag);
B=SoluteIndex{1,6}(SoluteFlag);
C=SoluteIndex{1,7}(SoluteFlag);
TLimit=[LowTemp,UpTemp];
Parameters=[A B C];

for j=1:length(T)
    if (0<(T(j)-LowTemp))&&(0<(UpTemp-T(j)))
        log10P=A-B/(T(j)+C);
        P(j)=100000*10^log10P;
    else if (nVarargs==1)
        P(j)=varargin{1};
        disp(['The temperature of ',num2str(j),'th data out of the temperarure range!']);
    else
        error('The temperature T is out of the temperarure range and the variable input argument is more then one!');
        end
    end
end

end
