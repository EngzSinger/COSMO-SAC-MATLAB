function [atom,segments,Vcosmo,Acosmo,atomNum,segNum]=readG09cosmo(CompName,Path)
% This function can read the CompName's atom information, segments
% information, COSMO volume (A^3),COSMO area(A^2) from G09's .cosmo files
% 此程序主要是处理来自G09生成的cosmo文件
% 返回Vcosmo、Acosmo的单位分别为A^3、A^2
% [FileName,PathName]=uigetfile('*.cosmo','Select the txt files');
%输入要处理的物质名称
%CompName='butanol'; %物质名称不要出现数字
%构造cosmo文件路径
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
u2A=0.529177249;   %将原子单位转化为A

%获取分子中的原子个数
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
       for i=1:4  %这里的数字4可以根据cosmo文件的具体内容进行更改
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
%获取分子原子个数结束
fseek(fid,0,-1);  %重新指向文件头
atom=zeros(atomNum,4);
segments=zeros(segNum,9);
while ~feof(fid)
    tline=fgetl(fid);
    %+读取cosmo面积和体积
    if strfind(tline,'$cosmo_data')
        for i=1:4
            tline=fgetl(fid);
            %++读取cosmo面积，并将原子单位转换为A
            if strfind(tline,'area')
                [str1,str2,Acosmo]=strread(tline,'%s%s%f');
                Acosmo=Acosmo*u2A*u2A;
            end
            %++读取cosmo面积结束
            
            %++读取cosmo体积，并将原子单位转换为A
            if strfind(tline,'volume=')
                [str1,Vcosmo]=strread(tline,'%s%f');
                Vcosmo=Vcosmo*u2A*u2A*u2A;
            end
            %++读取cosmo体积结束            
        end                  
    end
    
    %+读取原子坐标信息
    if strfind(tline,'!DATE')
            [C]=textscan(fid,'%s %f %f %f %s %d %s %s %f',atomNum);
            atom(:,2)=C{1,2};  %不知道G09生成的原子坐标的单位是啥？
            atom(:,3)=C{1,3};  %答：单位为A
            atom(:,4)=C{1,4};
            str=C{1,8};
            for j=1:atomNum
                str4=char(str(j));
                atom(j,1)=1*strcmp(str4,'H')+2*strcmp(str4,'B')+3*strcmp(str4,'C')+4*strcmp(str4,'N')+5*strcmp(str4,'O')+6*strcmp(str4,'F')+7*strcmp(str4,'Al')+8*strcmp(str4,'P')+9*strcmp(str4,'S')+10*strcmp(str4,'Cl')+11*strcmp(str4,'Fe')+12*strcmp(str4,'Zn')+13*strcmp(str4,'Ga')+14*strcmp(str4,'In');
            end
    end
    %+读取原子坐标信息结束
    
    %+读取片段信息
    if strfind(tline,'length')
       for i=1:4  %这里的数字4可以根据cosmo文件的具体内容进行更改
          tline=fgetl(fid); 
       end
       [C]=textscan(fid,'%d %d %f %f %f %f %f %f %f',segNum);
       for i=1:9
           segments(:,i)=C{1,i};
       end
       break
    end
   %+读取片段信息结束
end
fclose(fid);
%以txt文件输出原子信息atom和片段信息segments
%+构造原子信息txt文件路径
pathAtom=['/home/zxwu/WORK_SPACE/COSMO_SAC/ATOM/07/',CompName,'Atom.txt'];
%pathAtom=['/home/whtu/COSMO-MATLAB/ATOMS/07/',CompName,'Atom.txt'];
%+构造片段信息txt文件路径
pathSegments=['/home/zxwu/WORK_SPACE/COSMO_SAC/AVERAGE/07/',CompName,'Seg.txt'];
%pathSegments=['/home/whtu/COSMO-MATLAB/AVERAGE/07/',CompName,'Seg.txt'];
%格式化输出原子坐标信息
fid=fopen(pathAtom,'w');
fprintf(fid,'%3d %15.10f %15.10f %15.10f\n',atom');
fclose(fid);
%格式化输出原子片段信息
fid=fopen(pathSegments,'w');
fprintf(fid,'%5d  %5d  %20.9f  %20.9f  %20.9f  %20.9f  %20.9f  %20.9f  %20.9f\n',segments');
fclose(fid);
