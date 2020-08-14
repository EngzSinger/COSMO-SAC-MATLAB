function [area,volume,prf]=PurePrf(CompName,Path,cosmoFlag)
%This function get profile with two ways:
% if cosmoFlag=1(or 2 or 3), generate the profile from cosmo file
% if cosmoFlag>4, get profile from profile database
%areaΪ����������һ����ΪCOSMO�����ڶ�����Ϊ�����Ƭ�ε�����3����Ϊ����Ƭ�ε����A^2��
%prfΪSegNum��3�еľ��󣬵�1��Ϊsigma������꣩����2��3�зֱ�Ϊ���������ֵ�profile��A^2��
% global  ProfilePath;
if (cosmoFlag<4)
     %[area,volume,prf]=Cosmo2Profile(CompName,Path,cosmoFlag);
     [area,volume,prf]=Cosmo2Profile(CompName,Path,cosmoFlag);
     %��prf����������txt�ļ�
    %+����prf������·��
%     pathPrf=[ProfilePath,CompName,'prf.txt'];
%     fid=fopen(pathPrf,'w');
%     fprintf(fid,'%8.3f %15.8f %15.8f\n',prf');
%     fclose(fid);
    %prf���������
    
    %��prf����������excel�ļ�
    %+����prf������·��
    % pathPrf=['/home/whtu/COSMO-MATLAB/PROFILES/07/',CompName,'.xls'];
    % xlswrite(pathPrf,prf);
    %prf���������
    % plot(sigma,prf(:,2),'r');
    % hold on;
    % plot(sigma,prf(:,3),'b');
else
     [area,volume,prf]=ReadProfile(CompName,Path);
end

end








