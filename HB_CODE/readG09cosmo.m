function [atom,segments,Vcosmo,Acosmo,atomNum,segNum]=readG09cosmo(CompName,Path)
% This function can read the CompName's atom information, segments
% information, COSMO volume (A^3),COSMO area(A^2) from G09's .cosmo files
% �˳�����Ҫ�Ǵ�������G09���ɵ�cosmo�ļ�
% ����Vcosmo��Acosmo�ĵ�λ�ֱ�ΪA^3��A^2
% [FileName,PathName]=uigetfile('*.cosmo','Select the txt files');
%����Ҫ�������������
%CompName='butanol'; %�������Ʋ�Ҫ��������
%����cosmo�ļ�·��
[Index]=findIndex(Path,CompName);
 pathCOSMO=[Path,Index,'.cosmo'];
fid=fopen(pathCOSMO);
if fid<0
    error(['PLEASE CHECK WHETHER EXIST ',Index,'.cosmo (',CompName,'.cosmo) FILE IN ',Path]);
end
%   pathCOSMO=[Path,CompName,'.cosmo'];
% if fid<0
%     error(['PLEASE CHECK WHETHER EXIST ',CompName,'.cosmo FILE IN ',Path]);
% end
u2A=0.529177249;   %��ԭ�ӵ�λת��ΪA

%��ȡ�����е�ԭ�Ӹ���
while ~feof(fid)
   tline=fgetl(fid);
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
   if strfind(tline,'length')
       for i=1:4  %���������4���Ը���cosmo�ļ��ľ������ݽ��и���
          tline=fgetl(fid); 
       end
       countSegment=0;
       while ~feof(fid)
           tline=fgetl(fid);
           countSegment=countSegment+1;
       end
       break
   end
end 
atomNum=countAtom;
segNum=countSegment;
%��ȡ����ԭ�Ӹ�������
fseek(fid,0,-1);  %����ָ���ļ�ͷ
atom=zeros(atomNum,4);
segments=zeros(segNum,9);
while ~feof(fid)
    tline=fgetl(fid);
    %+��ȡcosmo��������
    if strfind(tline,'$cosmo_data')
        for i=1:4
            tline=fgetl(fid);
            %++��ȡcosmo���������ԭ�ӵ�λת��ΪA
            if strfind(tline,'area')
                [str1,str2,Acosmo]=strread(tline,'%s%s%f');
                Acosmo=Acosmo*u2A*u2A;
            end
            %++��ȡcosmo�������
            
            %++��ȡcosmo���������ԭ�ӵ�λת��ΪA
            if strfind(tline,'volume=')
                [str1,Vcosmo]=strread(tline,'%s%f');
                Vcosmo=Vcosmo*u2A*u2A*u2A;
            end
            %++��ȡcosmo�������            
        end                  
    end
    
    %+��ȡԭ��������Ϣ
    if strfind(tline,'!DATE')
            [C]=textscan(fid,'%s %f %f %f %s %d %s %s %f',atomNum);
            atom(:,2)=C{1,2};  %��֪��G09���ɵ�ԭ������ĵ�λ��ɶ��
            atom(:,3)=C{1,3};  %�𣺵�λΪA
            atom(:,4)=C{1,4};
            str=C{1,8};
            for j=1:atomNum
                str4=char(str(j));
                atom(j,1)=1*strcmp(str4,'H')+2*strcmp(str4,'B')+3*strcmp(str4,'C')+4*strcmp(str4,'N')+5*strcmp(str4,'O')+6*strcmp(str4,'F')+7*strcmp(str4,'Al')+8*strcmp(str4,'P')+9*strcmp(str4,'S')+10*strcmp(str4,'Cl')+11*strcmp(str4,'Fe')+12*strcmp(str4,'Zn')+13*strcmp(str4,'Ga')+14*strcmp(str4,'In');
            end
    end
    %+��ȡԭ��������Ϣ����
    
    %+��ȡƬ����Ϣ
    if strfind(tline,'length')
       for i=1:4  %���������4���Ը���cosmo�ļ��ľ������ݽ��и���
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
pathAtom=['/home/zxwu/WORK_SPACE/COSMO_SAC/ATOM/07/',CompName,'Atom.txt'];
%pathAtom=['/home/whtu/COSMO-MATLAB/ATOMS/07/',CompName,'Atom.txt'];
%+����Ƭ����Ϣtxt�ļ�·��
pathSegments=['/home/zxwu/WORK_SPACE/COSMO_SAC/AVERAGE/07/',CompName,'Seg.txt'];
%pathSegments=['/home/whtu/COSMO-MATLAB/AVERAGE/07/',CompName,'Seg.txt'];
%��ʽ�����ԭ��������Ϣ
fid=fopen(pathAtom,'w');
fprintf(fid,'%3d %15.10f %15.10f %15.10f\n',atom');
fclose(fid);
%��ʽ�����ԭ��Ƭ����Ϣ
fid=fopen(pathSegments,'w');
fprintf(fid,'%5d  %5d  %20.9f  %20.9f  %20.9f  %20.9f  %20.9f  %20.9f  %20.9f\n',segments');
fclose(fid);
