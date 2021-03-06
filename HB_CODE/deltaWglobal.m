function deltaW=deltaWglobal(m,k,sigma,parameterFlag,C,T)
%���������㽻���໥������
%����m��k�ֱ�Ϊsigma�����
%����tempPrfΪ��֣����ϣ���sigma-profile���������ֵ�ƴ�ӣ�����Ϊԭ�����ȵ�2������һ��Ϊsigma���ڶ���Ϊprofile����λΪ���λA**2
%parameterFlagΪѡ����㽻���໥���õı�־
%
global Aeff;
%�����໥�����ѡ��parameterFlag��=1��ȡ��COSMO-SAC(2007)�����£�=2��ȡ��Wang Shu��ʿ���ĵ�C++����=3��ȡ������Ĳ�ʿ����)
switch parameterFlag
    case 1
        fpol=0.6917;                        %��������,from the paper of COSMO-SAC(2004)��COSMO-SAC(2002)������ȡ0.64
        %Aeff=7.25;                          %��Ч��׼Ƭ�������Ƭ��ƽ������ֵ����һ��,from the paper of COSMO-SAC(2007)
        eps0=0.0002395;                     %��ս�糣��(e^2 mol)/((kcal)/(A) 
        alpha=fpol*0.3*(Aeff^1.5)/2.0/eps0; %����Ƭ���໥������ϵ����������̩����
        chb=3484.42;                        %���Ƭ���໥����ϵ��kcal/(mol A^4)/e^2, from the paper of COSMO-SAC(2007)�� 
    case 2
        alpha=8502.526564103386;             %����Wang Shu C++����
        chb=3484.423824950168;               %��������Wang Shu C++����
    case 3
        %alpha=18785.014;             %�������粩ʿ����
        %alpha=10072;             %�������粩ʿ����
	alpha=5891.6;
	chb=-131.7;
        %alpha=9392.51;
	%chb=3631.760;               %�������粩ʿ����
        %chb=2797;               %�������粩ʿ����
    case 4
		% Optimize two parameter simultaneously
%         disp('You should choice a suit of parameter to calculate exchange interaction!');
        alpha=C(1);
        chb=C(2);
	case 5
		% Optimize alpha only. The VLE systems should not Hydrogen Bond system
		alpha=C;		% C is a scalar
		chb=0;
	case 6
		% Optimize chb only. The VLE systems should contain Hydrogen Bond interation
%		alpha=1.422620466857909e+03;	% This value is optimized alpha based on case 5
	%	alpha=9392.51;
		alpha=1422.6;
		chb=C;			% C is a scalar
    case 7
        alpha=C(1)+C(2)/T/T;
        chb=C(3);
end
%disp(alpha);
%disp(chb);
countSigma=length(sigma);  %sigma��ݵ���
Ehb=0.0;                  %����ܳ�ʼ��

%�����µ�sigma
sigmaNew=[sigma;sigma];
%���㲻ƥ����
Emf=alpha*(sigmaNew(m)+sigmaNew(k))^2;
%���������
if (m>countSigma & k>countSigma & sigmaNew(m)*sigmaNew(k)<0.0)
Ehb=chb*(sigmaNew(m)-sigmaNew(k))^2;
end
deltaW=Emf-Ehb;
