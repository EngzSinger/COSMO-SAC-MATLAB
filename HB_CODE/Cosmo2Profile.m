function  [area,volume,prf]=Cosmo2Profile(CompName,Path,cosmoFlag,HB_Flag)
%This function first average the segments with a effective surface,Aeff
%And generate the sigma-profile, divided into two parts,nonhydrogent bond
%part and hydrogent bond part(belongs to N O F atoms and H connected with these atoms)
%1-H;2-B;3-C;4-N;5-O;6-F;7-Al;8-P;9-S;10-Cl;11-Fe;12-Zn;13-Ga;14-In
%constants
global Aeff;
Aeff=7.25;
sigRange=0.035;
interval=0.001;
pointNum=2*sigRange/interval+1;
sigma0=0.007;
% Aeff=7.5;  %有效片段面积，for COSMO-SAC(2002)
%Aeff=7.25;  %Effective Fragment area，for COSMO-SAC(2007)
% f_decay=1; %衰退因子,for COSMO-SAC(2002)
f_decay=3.57; %衰退因子，for COSMO-SAC(2007)
u2u=0.529177249; %Unit transforming factor
% totatom=3;       %总原子数
rCovalent=1.5;   
%load the informations of atoms
switch cosmoFlag
    case 1
        [atom,segments,volume,Acosmo,totatom,segNum]=readThermCOSMO(CompName,Path); %get the information of compound
    case 2
        [atom,segments,volume,Acosmo,totatom,segNum]=readVTcosmo(CompName,Path);
    otherwise
         [atom,segments,volume,Acosmo,totatom,segNum]=readG09cosmo(CompName,Path);
end
atom=[atom,zeros(totatom,1)];
% atomConnect=zeros(totatom,totatom);
%distinguish hydrogent bond atom

%hb_index=[find(atom(:,1)==4);find(atom(:,1)==5);find(atom(:,1)==6)];
%atom(hb_index,5)=-1;
%H_index=find(atom(:,1)==1);
%for i=1:length(H_index)
%    for j=1:length(hb_index)        
%        if norm(atom(H_index(i),2:4)-atom(hb_index(j),2:4))<rCovalent %norm(A):求向量A的长度
%            atom(H_index(i),5)=1;  %确定与N O F原子相连的H原子
%        end
%    end
%end
switch HB_Flag
	%normal HB mode
	case 1
		hb_index=[find(atom(:,1)==4);find(atom(:,1)==5);find(atom(:,1)==6)];
		atom(hb_index,5)=-1;
		H_index=find(atom(:,1)==1);
		for i=1:length(H_index)
			for j=1:length(hb_index)
        			if norm(atom(H_index(i),2:4)-atom(hb_index(j),2:4))<rCovalent 
            				atom(H_index(i),5)=1;  
        			end
    			end
		end
	%remove F atom
	case 2
		
                hb_index=[find(atom(:,1)==4);find(atom(:,1)==5)];
                atom(hb_index,5)=-1;
                H_index=find(atom(:,1)==1);
                for i=1:length(H_index)
                        for j=1:length(hb_index)
                                if norm(atom(H_index(i),2:4)-atom(hb_index(j),2:4))<rCovalent
                                        atom(H_index(i),5)=1;
                                end
                        end
                end
	%add some H atom 
	case 3
		hb_index=[find(atom(:,1)==4);find(atom(:,1)==5);find(atom(:,1)==6)];
		N_index=[find(atom(:,1)==4);find(atom(:,1)==10)];
		C_index=[find(atom(:,1)==3)];
		%Cl_index=[find(atom(:,1)==10);
		atom(hb_index,5)=-1;
		H_index=find(atom(:,1)==1);
		for i=1:length(H_index)
			for l=1:length(hb_index)
                                if norm(atom(H_index(i),2:4)-atom(hb_index(l),2:4))<rCovalent
                                        atom(H_index(i),5)=1;
                                end
                        end
			for j=1:length(C_index)
				if norm(atom(H_index(i),2:4)-atom(C_index(j),2:4))<rCovalent
					count=0;
					for k=1:length(N_index)
						if norm(atom(N_index(k),2:4)-atom(C_index(j),2:4))<rCovalent
							count=count+1;
						end
					end
					if count>1
						atom(H_index(i),5)=1;
					end
				end
			end
		end
	%add H remove F
	case 4
		hb_index=[find(atom(:,1)==4);find(atom(:,1)==5)];
                N_index=[find(atom(:,1)==4);find(atom(:,1)==10)];
                C_index=[find(atom(:,1)==3)];
                %Cl_index=[find(atom(:,1)==10);
                atom(hb_index,5)=-1;
                H_index=find(atom(:,1)==1);
                for i=1:length(H_index)
			for l=1:length(hb_index)
                                if norm(atom(H_index(i),2:4)-atom(hb_index(l),2:4))<rCovalent
                                        atom(H_index(i),5)=1;
                                end
                        end
                        for j=1:length(C_index)
                                if norm(atom(H_index(i),2:4)-atom(C_index(j),2:4))<rCovalent
                                        count=0;
                                        for k=1:length(N_index)
                                                if norm(atom(N_index(k),2:4)-atom(C_index(j),2:4))<rCovalent
                                                        count=count+1;
                                                end
                                        end
                                        if count>1
                                                atom(H_index(i),5)=1;
                                        end
                                end
                        end
                end
end
format long;
% segments=load(path_segments);
%transform the unit of position of segments
segments(:,3:5)=segments(:,3:5)*u2u;
%average segments
segments_avg=segments;
r_eff2=Aeff/pi;
r_n2=segments(:,7)/pi;
% segNum=length(segments(:,1));  %Segments number
for i=1:segNum
    dmn=zeros(segNum,3);
    position=repmat(segments(i,3:5),segNum,1);  %segment_m的坐标
    dmn=(position-segments(:,3:5));
    dmn=dmn.^2;
    dmn2=sum(dmn')';
    num=sum(segments(:,8).*(r_n2*r_eff2./(r_n2+r_eff2)).*exp(-f_decay*(dmn2./(r_n2+r_eff2))));
    denom=sum((r_n2*r_eff2./(r_n2+r_eff2)).*exp(-f_decay*(dmn2./(r_n2+r_eff2))));
     segments_avg(i,8)=num/denom;
end
%generate profile
sigma=-sigRange:interval:sigRange;
sigma=sigma';
prf_nhb=zeros(pointNum,1);
prf_hb=zeros(pointNum,1);
for i=1:segNum
    temp=floor((segments_avg(i,8)-sigma(1))/interval);
    lval=segments_avg(i,7)*(sigma(temp+2)-segments_avg(i,8))/interval;
    hval=segments_avg(i,7)*(segments_avg(i,8)-sigma(temp+1))/interval;
    if (atom(segments_avg(i,2),5)*segments_avg(i,8)<0)
        prf_hb(temp+1)=prf_hb(temp+1)+lval;
        prf_hb(temp+2)=prf_hb(temp+2)+hval;
    else
         prf_nhb(temp+1)=prf_nhb(temp+1)+lval;
        prf_nhb(temp+2)=prf_nhb(temp+2)+hval;
    end
end
%not yet multiply hydrogent bond probability
prf=[sigma,prf_nhb,prf_hb];
%multiply hydrogent bond probability
Phb=1-exp(-sigma.^2/2/sigma0/sigma0);
% [sigma,Phb]
prf(:,2)=prf_nhb+prf_hb.*(1-Phb);
prf(:,3)=prf_hb.*Phb;

area_nhb=sum(prf(:,2));
area_hb=sum(prf(:,3));
area=[area_nhb+area_hb,area_nhb,area_hb];

end
