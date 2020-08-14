function [atom,segments,Vcosmo,Acosmo,atomNum,segNum]=readVTcosmo(CompName)
% This function can read the CompName's atom information, segments
% information, COSMO volume (A^3),COSMO area(A^2) from VT's .cosmo files
% �˳�����Ҫ�Ǵ�������VT��cosmo�ļ�
% ����Vcosm��Acosmo�ĵ�λ�ֱ�ΪA^3��A^2
% ����Ҫ������������
% CompName='ButyricAcid'; %������Ʋ�Ҫ��������
%ͳ�Ƴ������ʱ�俪ʼ
% tic;
%����cosmo�ļ�·��
[Index]=findIndex(Path,CompName);
% pathCOSMO=['E:\MATLAB\COSMO\VT\',CompName,'.cosmo'];
pathCOSMO=[Path,Index,'.cosmo'];
fid=fopen(pathCOSMO);
if fid<0
    error(['PLEASE CHECK WHETHER EXIST ',Index,'.cosmo (',CompName,'.cosmo) FILE IN ',Path]);
end
%��ȡ�����е�ԭ�Ӹ����Ƭ�θ���
while ~feof(fid)
   tline=fgetl(fid);
   %+��ȡ�����е�ԭ�Ӹ���
   if strfind(tline,'!DATE')
       countAtom=0;
       while ~feof(fid)
           tline=fgetl(fid);
           if strfind(tline,'end')
                break;
           end
            countAtom=countAtom+1;
       end
   end
   %+��ȡ������ԭ�Ӹ������
   
   %+��ȡ������Ƭ����Ŀ
   if strfind(tline,'segments:')
       [str1,str2,str3,str4,segNum]=strread(tline,'%s%s%s%s%f'); %total number of segments:    958
       
%        break
   end
   %+��ȡ������Ƭ����Ŀ����
end 
atomNum=countAtom;
fseek(fid,0,-1);  %����ָ���ļ�ͷ
atom=zeros(atomNum,4);
segments=zeros(segNum,9);
while ~feof(fid)
    tline=fgetl(fid); 
%     Total surface area of cavity (A**2)     =    88.40645
% 
%   Total volume of cavity (A**3)           =    70.19948
    %+��ȡcosmo���
    if strfind(tline,'(A**2)')
        [str1,str2,str3,str4,str5,str6,str7,Acosmo]=strread(tline,'%s%s%s%s%s%s%s%f');
    end
    %+��ȡcosmo������
    
    %+��ȡcosmo���
    if strfind(tline,'(A**3)')
        [str1,str2,str3,str4,str5,str6,Vcosmo]=strread(tline,'%s%s%s%s%s%s%f');
    end
    %+��ȡcosmo������
    
    %+��ȡԭ�������Ϣ
    if strfind(tline,'!DATE')
            [C]=textscan(fid,'%s %f %f %f %s %d %s %s %f',atomNum);
            atom(:,2)=C{1,2};
            atom(:,3)=C{1,3};
            atom(:,4)=C{1,4};
            str=C{1,8};
            for j=1:atomNum
                str4=char(str(j));
                atom(j,1)=1*strcmp(str4,'H')+2*strcmp(str4,'B')+3*strcmp(str4,'C')+4*strcmp(str4,'N')+5*strcmp(str4,'O')+6*strcmp(str4,'F')+7*strcmp(str4,'Al')+8*strcmp(str4,'P')+9*strcmp(str4,'S')+10*strcmp(str4,'Cl')+11*strcmp(str4,'Fe')+12*strcmp(str4,'Zn')+13*strcmp(str4,'Ga')+14*strcmp(str4,'In');
            end
    end
    %+��ȡԭ�������Ϣ����
    
    %+��ȡƬ����Ϣ
    if strfind(tline,'charge/area')
       for i=1:2  %���������2���Ը��cosmo�ļ��ľ������ݽ��и��
          tline=fgetl(fid); 
       end
       [C]=textscan(fid,'%d %d %f %f %f %f %f %f %f',segNum);
       for i=1:9
           segments(:,i)=C{1,i};
       end
       break
    end
   %+��ȡƬ����Ϣ����
end
fclose(fid);
%��txt�ļ����ԭ����Ϣatom��Ƭ����Ϣsegments
%+����ԭ����Ϣtxt�ļ�·��
pathAtom=['/home/zxwu/WORK_SPACE/MATLAB_COSMO/ATOMS/07/',CompName,'Atom.txt'];
%+����Ƭ����Ϣtxt�ļ�·��
pathSegments=['/home/zxwu/WORK_SPACE/MATLAB_COSMO/AVERAGE/07/',CompName,'Seg.txt'];
%��ʽ�����ԭ�������Ϣ
fid=fopen(pathAtom,'w');
fprintf(fid,'%3d %15.10f %15.10f %15.10f\n',atom');
fclose(fid);
%��ʽ�����ԭ��Ƭ����Ϣ
fid=fopen(pathSegments,'w');
fprintf(fid,'%5d  %5d  %20.9f  %20.9f  %20.9f  %20.9f  %20.9f  %20.9f  %20.9f\n',segments');
fclose(fid);
%ͳ�Ƴ��������ʱ�����
% timeCost=toc
