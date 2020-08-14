function [area,volume,prf,n] = Cosmo2Profie(CompName,Path,para)
% This function can transform segment information from cosmo file to profile

sigRange = 0.035;
interval = 0.001;
pointNum = 71;
sigma0 = 0.007;

% Read global parameter from a file 
%para = Para();

% note that f_decay value from file is string type
% data type conversion must be completed
f_decay = str2num(para('f_decay'));
u2u = str2num(para('u2A'));
rCovalent = str2num(para('rCovalent'));
Aeff = str2num(para('Aeff_2007'));


[atom,segments,volume,Acosmo,atomNum,segNum] = ReadG09(CompName,Path,para);
atom = [atom,zeros(atomNum,1)];

H_Index = [find(atom(:,1)==1)];
A_Index = [find(atom(:,1)==4);find(atom(:,1)==5);find(atom(:,1)==6)];

atom(A_Index,5) = -1;
for i = 1:length(H_Index)
	for j = 1:length(A_Index)
		% check if there is a bond formation between two atom
		if norm(atom(H_Index(i),2:4)-atom(A_Index(j),2:4)) < rCovalent
			atom(H_Index(i),5)=1;
		end
	end
end

format long;
segments(:,3:5) = segments(:,3:5)*u2u;
segment_avg = segments;
r_n2 = segments(:,7)/pi;

% reff value needs to be discussed
% Aeff should be the average area of segments on the molecular surface,not the real Aeff that refers to effective contact surface area
%f_decay is factor to rescale dmn2 in order to allow units of dmn2 agree with that of r_n2
%f_decay=3.57 is the case for previous research, while here f_decay sholud be 1/3.57

%Aeff = Acosmo/segNum;
%f_decay = 1/f_decay;
r_eff2 = Aeff/pi;



for i = 1:segNum
	dmn  = zeros(segNum,3);
	position = repmat(segments(i,3:5),segNum,1);
	dmn = position - segments(:,3:5);
	dmn = dmn.^2;
	dmn2 = sum(dmn')';
	num = sum(segments(:,8).*(r_n2*r_eff2./(r_n2+r_eff2)).*exp(-f_decay*(dmn2./(r_n2+r_eff2))));
	denum = sum((r_n2*r_eff2./(r_n2+r_eff2)).*exp(-f_decay*(dmn2./(r_n2+r_eff2))));
	%num/denum
	segment_avg(i,8) = num/denum;
end

sigma = -sigRange:interval:sigRange;
sigma = sigma';
prf_hb = zeros(pointNum,1);
prf_nhb = zeros(pointNum,1);
for i = 1:segNum
	local = floor((segment_avg(i,8)+sigRange)/interval);
	% a reverse order must be considered
	right_val = segment_avg(i,7)*(segment_avg(i,8)-sigma(local+1))/interval;
	left_val = segment_avg(i,7)*(sigma(local+2)-segment_avg(i,8))/interval;
	if (atom(segment_avg(i,2),5)*segment_avg(i,8)<0)
		prf_hb(local+1) = prf_hb(local+1)+left_val;
		prf_hb(local+2) = prf_hb(local+2)+right_val;
	else
		prf_nhb(local+1) = prf_nhb(local+1)+left_val;
                prf_nhb(local+2) = prf_nhb(local+2)+right_val;
	end
end
prf = [sigma,prf_nhb,prf_hb];
%prf
Phb = 1-exp(-sigma.^2/2/sigma0/sigma0);
prf(:,2) = prf_nhb+prf_hb.*(1-Phb);
prf(:,3) = prf_hb.*Phb;

area_nhb = sum(prf(:,2));
area_hb = sum(prf(:,3));
area = [area_nhb+area_hb,area_nhb,area_hb];
n = segNum;

end
