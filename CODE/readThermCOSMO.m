function [atom,segments,Vcosmo,Acosmo,atomNum,segNum]=readThermCOSMO(CompName,Path)
% �˳�����Ҫ�Ǵ�������COSMOtherm��cosmo�ļ�
% ����Vcosmo�ĵ�λΪA^3
% ����Acosmo�ĵ�λΪA^2
% hang=[4 17];
% lie=[3 8];
% [FileName,PathName]=uigetfile('*.cosmo','Select the txt files');
%����Ҫ������������
%CompName='butanol'; %������Ʋ�Ҫ��������
[Index]=findIndex(Path,CompName);
% pathCOSMO=['E:\MATLAB\COSMO\THERM\',CompName,'.cosmo'];
 pathCOSMO=[Path,Index,'.cosmo'];
fid=fopen(pathCOSMO);
if fid<0
    error(['PLEASE CHECK WHETHER EXIST ',Index,'.cosmo (',CompName,'.cosmo) FILE IN ',Path]);
end
u2A=0.529177249;   %��ԭ�ӵ�λת��ΪA
% temp=textscan(fid,'%s %s %s %s %s %s %s %s');
% fclose(fid);
% for i=1:(hang(2)-hang(1)+1)
%     for j=1:(lie(2)-lie(1)+1)
%         b=temp{j+lie(1)-1}{i+hang(1)-1};
%         A(i,j)=str2num(b);
%     end
% end
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
       for i=1:4  %���������4���Ը��cosmo�ļ��ľ������ݽ��и��
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
%��ȡ����ԭ�Ӹ������
fseek(fid,0,-1);  %����ָ���ļ�ͷ
atom=zeros(atomNum,4);
segments=zeros(segNum,9);
while ~feof(fid)
    tline=fgetl(fid);
    if strfind(tline,'$cosmo_data')
        [C]=textscan(fid,'%s %f',8);  %���������8�Ǹ��cosmo�ļ��ṹ����ģ���ͬ��cosmo�ļ���Ҫ����Ӧ�ĵ���                                      
        c=C{1,2};
        Acosmo=c(7)*u2A*u2A;  %ԭ�ӵ�λ���������ֵ7�Ǹ�������ݿ�cosmo�ļ���ʽ���ģ�          
        Vcosmo=c(8)*u2A*u2A*u2A;  %ԭ�ӵ�λ���������ֵ7�Ǹ�������ݿ�cosmo�ļ���ʽ����
                  
    end
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
    if strfind(tline,'length')
       for i=1:4  %���������4���Ը��cosmo�ļ��ľ������ݽ��и��
          tline=fgetl(fid); 
       end
       [C]=textscan(fid,'%d %d %f %f %f %f %f %f %f',segNum);
       for i=1:9
           segments(:,i)=C{1,i};
       end
       break
   end
end
fclose(fid);
%��txt�ļ����ԭ����Ϣatom��Ƭ����Ϣsegments
%+����ԭ����Ϣtxt�ļ�·��
%pathAtom=['/home/whtu/COSMO-MATLAB/ATOMS/07/',CompName,'Atom.txt'];
pathAtom=['/home/zxwu/WORK_SPACE/MATLAB_COSMO/ATOMS/07/',CompName,'Atom.txt'];
%+����Ƭ����Ϣtxt�ļ�·��
%pathSegments=['/home/whtu/COSMO-MATLAB/AVERAGE/07/',CompName,'Seg.txt'];
pathSegments=['/home/zxwu/WORK_SPACE/MATLAB_COSMO/AVERAGE/07/',CompName,'Seg.txt'];
%��ʽ�����ԭ�������Ϣ
fid=fopen(pathAtom,'w');
fprintf(fid,'%3d %15.10f %15.10f %15.10f\n',atom');
fclose(fid);
%��ʽ�����ԭ��Ƭ����Ϣ
fid=fopen(pathSegments,'w');
fprintf(fid,'%5d  %5d  %20.9f  %20.9f  %20.9f  %20.9f  %20.9f  %20.9f  %20.9f\n',segments');
fclose(fid);
