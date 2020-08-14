function [vap,para,vapFlag]=vapAntoineWagner(vapPath,SoluteName,T,varargin)
% CompName: Compund Name to calculate, e.g. Methanol
% T: Temperature to calculate, K, e.g. 303.15 or [298.15 303.15]
%% OPEN THE .vap FILE TO READ ANTOINE COEFFICIENTS AND WAGNER COEFFICIENTS
% SoluteName='Methanol';
vapFlag=0;
para=[];
fid=fopen([vapPath,SoluteName,'.vap']);
if fid<0
    error(['Fail to open ',SoluteName,'.vap. Please check whether this file exist in ',vapPath,' !']);
end
% vap=0;
while ~feof(fid)
    tl=fgetl(fid);
    if strfind(tl,'Antoine')
        AntoineNum=fscanf(fid,'%f',1);   
        tl=fgetl(fid);
        tl=fgetl(fid);
        if AntoineNum>0
            AntoinePara=textscan(fid,'%f %f %f %f %f',AntoineNum);
        end
    end
    if strfind(tl,'Wagner')
        WagnerNum=fscanf(fid,'%f',1);
        tl=fgetl(fid);
        tl=fgetl(fid);
        if WagnerNum>0
            WagnerPara=textscan(fid,'%f %f %f %f %f %f %f %f',WagnerNum);
        end
    end    
end
fclose(fid);
%% CLACULATING THE VAPOR PRESSURE
    Tmin=min(T);
    Tmax=max(T);
% ANTOINE EQUATION
for i=1:AntoineNum
    LowTemp=AntoinePara{1,1}(i);
    UpTemp=AntoinePara{1,2}(i);
    A=AntoinePara{1,3}(i);
    B=AntoinePara{1,4}(i);
    C=AntoinePara{1,5}(i);
    if (0<(Tmin-LowTemp) && 0<(UpTemp-Tmax))
        log10P=A-B./(T+C);
        vap=100*10.^log10P;  % Unit is KPa
        vapFlag=1;
        para=[LowTemp UpTemp A B C];
        break;
    end
end

% WAGNER EQUATION
if vapFlag<1
    for i=1:WagnerNum
       LowTemp=WagnerPara{1,1}(i);
       UpTemp=WagnerPara{1,2}(i);
       A=WagnerPara{1,3}(i);
       B=WagnerPara{1,4}(i);
       C=WagnerPara{1,5}(i);
       D=WagnerPara{1,6}(i);
       E=WagnerPara{1,7}(i);
       F=WagnerPara{1,8}(i);       
       if (0<(Tmin-LowTemp) && 0<(UpTemp-Tmax))
            t=1-T./B;
            lnP=log(A)+1./(1-t).*(C*t+D*t.^1.5+E*t.^3+F*t.^6);
            vap=exp(lnP);  % Unit is KPa
            vapFlag=2;
            para=[LowTemp UpTemp A B C D E F];
            break;
        end
    end
end

% EXPERIMENTAL PRESSURE
if ((vapFlag<1) &&(1==length(varargin)) && (varargin{1}>0))
    vap=varargin{1};
else if (vapFlag<1)
%         varargin{1}
        error('The temperature T is out of the temperarure range and there is no experimental pressure input!')
    end
end

end











