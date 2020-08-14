function  [area,volume,prf]=ReadProfile(CompName,Path)
% This function read the profile data from file CompNameprf.txt, e.g.
% cation001prf.txt
% You can find the CompNameprf.txt in the direction:   E:\MATLAB\PROFILES\07
% arae is a vector consists of three parts cosmo area of CompName:[AreaSum AreaNHB AreaHB]
% volume is a scalar to storge the value of cosmo volume
% 
% The format of cation001prf.txt as folows:
% Name:  (2-ethoxyethyl)-ethyl-dimethylammonium_cation
%    AreaSum(A^2)    AreaNHB(A^2)     AreaHB(A^2)      Voume(A^3)
%      205.475210      202.994211        2.480999      218.507295
% The Information of profile is listed as follows:
%    sigma(e/A^2)     PrfNHB(A^2)      PrfHB(A^2)
%          -0.035      0.00000000      0.00000000
%          -0.034      0.00000000      0.00000000
%          -0.033      0.00000000      0.00000000
%          -0.032      0.00000000      0.00000000
%               :                   :                       :
%               :                   :                       :
pointNum=71;
area=zeros(3,1);
prf=zeros(pointNum,3);
[Index]=findIndex(Path,CompName);
fid=fopen([Path,Index,'prf.txt']);
if fid<0
    error(['PLEASE CHECK WHETHER EXIST ',Index,'prf.txt (',CompName,'prf.txt) FILE IN ',Path]);
end
while ~feof(fid)
    tline=fgetl(fid);
    if strfind(tline,'AreaSum')
        [C]=textscan(fid,'%f %f %f %f',1);
        area(1)=C{1,1};
        area(2)=C{1,2};
        area(3)=C{1,3};
        volume=C{1,4};
    end
    if strfind(tline,'sigma')
        [C]=textscan(fid,'%f %f %f',pointNum);
        prf(:,1)=C{1,1};
        prf(:,2)=C{1,2};
        prf(:,3)=C{1,3};
    end
end
fclose(fid);
end