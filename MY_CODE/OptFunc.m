function aard = OptFunc(ac,para,listcell)

alpha = ac(1)
chb = ac(2)
Nsys = length(listcell{1,1});
Nsys

ard = zeros(Nsys,1);

Psat = @(T) 10^(5.20409-(1581.341/(T-33.50)))*100;

Mypar = parpool('local',24);
parfor i = 1:Nsys
	T = listcell{1,4}(i);
	solute = listcell{1,1}{i};
	cation = listcell{1,2}{i};
	anion = listcell{1,3}{i};
	Pexp = listcell{1,6}(i);
	x0 = listcell{1,5}(i);
	x = [x0,1-x0];
	Pcal = VLE(solute,cation,anion,x,T,para,'Alpha',alpha,'Chb',chb);
	disp(Pexp)
	disp(Pcal)	
	ard(i) = abs(Pcal-Pexp)/Pexp;
end
ard
disp(ard)
aard = sum(ard)/Nsys
delete(Mypar);
end
