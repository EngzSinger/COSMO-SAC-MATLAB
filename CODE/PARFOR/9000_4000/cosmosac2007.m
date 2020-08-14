% 1. The format of calling PurePrf() function: [area,volume,prf]=PurePrf(CompName,Path,cosmoFlag)
% (1) Input variable CompName is the name of substance, e.g. Water
% (2) Output variable area is a cosmo area vector of a substance: [SumArea,non-HB arae,HB area]
% (3) Output variable prf is a sigma-profile 2-D matrix,[sigma,non-HB profile(A^2),HB profile(A^2)]
%2. Calling format of MixPrf() function: [area1,area2,mix_area,prf1,prf2,mix_prf]=MixPrf(comp1,comp2,X)
% (1) PurePrf() function is called in MixPrf() function
% (2) Input variables comp1 and comp2 corresponding to the substances' names, e.g. Water, Methanol
% (3) Input variable X is a molar fraction vector for substances, e.g X=[0.42,0.58]
% 3. [gamma,lngamma]=gammares2007(Vcosmo,Acosmo,X,prf1,prf2,T,parameterFlag)
% (1)Vcosmo is a cosmo volume vector
% (2) Acosmo is a cosmo area vector
% (3) Input variable X is a molar fraction vector for substances, e.g. X=[0.42,0.58]
% (4) prf1 and prf2 are the profiles of substance1 and substance2,respectively. The unit is A^2
% (5) T is system temperature, Unit is K, e.g. 298.15K
% (6) parameterFlag is the flag to determine which interation parameters are
%       used in the DeltaW() function
% (7) Output variables gamma and lngamma are activity coefficients and
%        their natural logarithm
% 4. [gamma,lngamma]=gammacomb(Vcosmo,Acosmo,x,flag)
% (1)Vcosmo is a cosmo volume vector
% (2) Acosmo is a cosmo area vector
% (3) Input variable X is a molar fraction vector for substances, e.g. X=[0.42,0.58]
% (4)flag is aan integer to determine the method to calculate combinatial activity coefficients
%    flag=0:STAVERMAN-GUGGENHEIM equation
%    flag=1:UNIFAC equation
%

clear;clc;

%% constants for choose the model and parameters
global combFlag parameterFlag cosmoFlag Aeff;

%+ Set the value of Aeff.The vaue should be changed if it was optimized
Aeff=7.25;

%+ The flag for chosing equation to calculate combinatorial activity coefficients
combFlag=1;
%+Choice for system types, ILflag��=1: ILs containing; =0: Non ILs )
ILflag=1;
%+ Choice of parameters for calculating interaction of segments in deltaW() function, parameterFlag(=1: from COSMO-SAC(2007) paper; =2: from Wang Shu C++ program; =3: from Yang Li ph.D dissertation)
parameterFlag=4;

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
        alpha=10072;             %�������粩ʿ����
        chb=2797;               %�������粩ʿ����
    otherwise
%         disp('You should choice a suit of parameter to calculate exchange interaction!');
%          ���Ż���Ĳ�������
          alpha=9392.51;
          chb=3631.76;
end
        
%+Choice for cosmo file generation; cosmoFlag(=1: from COSMOThermX database; =2: from VT-2005; =3: from Gaussian09 D.01 generation; =4: Read profile files)
cosmoFlag=3;
%+select the calculation of thermodynamic stye,1-vle;2-gamma and plot the
%graphic
thermo=1;
%% Index for the program
global SystemVle;
%% Direction for Files
global ExpDataPath CalcDataPath;
ExpDataPath='/home/whtu/COSMO-MATLAB/EXPDATA/VLE/';
CalcDataPath='/home/whtu/COSMO-MATLAB/RESULTS/VLE/qc/txt/20190202/';
VLEpngPath='/home/whtu/COSMO-MATLAB/RESULTS/VLE/qc/png/20190202/';
CationProfilePath='/home/whtu/COSMO-MATLAB/PROFILES/07/qc/compound/cation/6-311gd/';
AnionProfilePath='/home/whtu/COSMO-MATLAB/PROFILES/07/qc/compound/anion/6-311pgd/';
CationCosmoPath='/home/whtu/COSMO-MATLAB/COSMO/G09/Cation/6-311Gd/';
AnionCosmoPath='/home/whtu/COSMO-MATLAB/COSMO/G09/Anion/6-311pGd/';
SoluteCosmoPath='/home/whtu/COSMO-MATLAB/COSMO/G09/Solute/6-311Gd/';
vapPath='/home/whtu/COSMO-MATLAB/VAP/';
VleListPath='/home/whtu/COSMO-MATLAB/CODE/inFile/';
VleListFileName='SystemVleList1.txt';
%% PARAMETERS FOR FIGURES
markersize=2;
fontsize=12;
linewidth=1.5;
% begin to record the time
tic;
switch thermo
 
    case 1
     %% Calculate the VLE data
    %Read the VLE systerms
    SystemVle=ReadSystemVleList([VleListPath,VleListFileName]);
    % the format of SystemVleList.txt
    % VleDataIndex Solute Cation Anion Temperature SaturationVaporPressure
    Nsys=length(SystemVle{1,1});   % GET THE NUMBER OF SYSTERMS
    for i=1:Nsys  % A LOOP FOR CALCULATING EACH SYSTEM
        tic;
        disp([SystemVle{1,1}{i},' IS CALCULATING..............']);
         T=SystemVle{1,5}(i);  % SYSTEM TEMPERATURE, K
         comp1=SystemVle{1,2}{i}; %   Solute Name e.g. Metanol
         %+COMPONENT 1 PROFILE GENERATION
         [area1,volume1,prf1]=PurePrf(comp1,SoluteCosmoPath,cosmoFlag);
              
    switch ILflag
        case 1            
            cation=SystemVle{1,3}{i};     % Cation Abbreviation e.g. C4MIm
            anion=SystemVle{1,4}{i};      % Aation Abbreviation e.g. Tf2N
            comp2=['[',cation,'][',anion,']'];       %  Ionic Liquids  Abbreviation e.g. [C4MIm][Tf2N]        
            %+ILs PROFILES GENERATION
            if cosmoFlag<4
                CationPath=CationCosmoPath;
                AnionPath=AnionCosmoPath;
            else
                CationPath=CationProfilePath;
                AnionPath=AnionProfilePath;
            end
            [Acation,Vcation,prf_cation]=PurePrf(cation,CationPath,cosmoFlag);
            [Aanion,Vanion,prf_anion]=PurePrf(anion,AnionPath,cosmoFlag);
            % GENERATING THE ILs' PROFILES BY TREATING ILs AS IONIC PAIRS 
            area2=Acation+Aanion;
            volume2=Vcation+Vanion;
            prf2=prf_cation+prf_anion;
            prf2(:,1)=prf2(:,1)/2.0;
        otherwise
            comp2=SystemVle{1,3}{i};            
            %+OOMPONENT 2 PROFILE GENERATION
            [area2,volume2,prf2]=PurePrf(comp2,SoluteCosmoPath,cosmoFlag);
    end
    
    %+Read experimental data
    ExpPath=[ExpDataPath,SystemVle{1,1}{i},'.txt'];
    ExpData=load(ExpPath);
    if (T<0)
        T=ExpData(:,1);
        x1=ExpData(:,2);
        Pexp=ExpData(:,3);
    else
        x1=ExpData(:,1);
        Pexp=ExpData(:,2);
    end      
    X=[x1,1-x1];
	
%     Calculate the solute saturation vapor pressure
%           [Psat,TLimit,Parameters]=antoine( comp1,T,SystemVle{1,6}(i));         % Calculate the Saturation Vapor Pressure of Solute    
          [Psat,Parameters,vapFlag]=vapAntoineWagner(vapPath,comp1,T,SystemVle{1,6}(i));
    %cosmo area and volume matrices construction
    Acosmo=[area1(1),area2(1)];
    Vcosmo=[volume1,volume2];
    % Call gammares2007() to calculate the residual activity coefficient
     [gammaR,lngammaR]=gammares2007(Vcosmo,Acosmo,X,prf1,prf2,T,parameterFlag);
    %Call gammacomb() to calculate the combinatorial acticvity coefficient
     [gammaC,lngammaC]=gammacomb(Vcosmo,Acosmo,X,combFlag);
   
    lngamma=lngammaR+lngammaC;
    gamma=exp(lngamma);       
    activity=x1.*gamma(:,1);
    activity(activity>1)=1;
    Pcal=activity.*Psat;  
    [PosDataNum,P_RAAD]=RAAD(X,Pexp,Pcal);
    time=toc;
   %% Output the results
   optKeyword_cation='# freq opt b3lyp/6-311G(d)';
   optKeyword_anion='# freq opt b3lyp/6-311+G(d)';
   optKeyword_solute='# freq opt b3lyp/6-311G(d)';
   cosmoKeyword_cation='#p opt b3lyp/6-311g(d) int=ultrafine';
   cosmoKeyword_anion='#p opt b3lyp/6-311+g(d) int=ultrafine';
   cosmoKeyword_solute='#p opt b3lyp/6-311g(d) int=ultrafine';
    ResultPath=[CalcDataPath,SystemVle{1,1}{i},'.txt'];
    fid=fopen(ResultPath,'w');
    
    fprintf(fid,'%s %s\n%s %f %s\n','Date: ',date,'Time expend: ',time,'Second.');
    fprintf(fid,'%s\n','The geometry optimization keywrds for each compunds:');
    fprintf(fid,'%10s %40s %40s\n','Compund','optKeyword','cosmoKeyword');
    fprintf(fid,'%10s %40s %40s\n','Cation',optKeyword_cation,cosmoKeyword_cation);
    fprintf(fid,'%10s %40s %40s\n','Anion',optKeyword_anion,cosmoKeyword_anion);
    fprintf(fid,'%10s %40s %40s\n\n','Solute',optKeyword_solute,cosmoKeyword_solute);
    if length(T)>1
         fprintf(fid,'%s%s%s%s%s%10.3f%s%10.3f%s\n\n','VLE results of the systerm:',comp1,'(1)+',comp2,'(2) at ',min(T),'-',max(T),' K.');
    else       
        fprintf(fid,'%s%s%s%s%s%10.3f%s\n\n','VLE results of the systerm:',comp1,'(1)+',comp2,'(2) at ',T,' K.');
    end
    % output how to get the vapor pressure
    switch vapFlag
        case 0
            fprintf(fid,'%s\n','The Saturation Vapor Pressure from experiments!');
        case 1
            fprintf(fid,'%s\n','The Saturation Vapor Pressure calculated with  ANTOINE EQUATION [ log10P=A-B/(T+C) ], The parameters are: ');
            fprintf(fid,'%s %8.5f %10.3f %8.3f %s %7.2f %s %7.2f %s\n','[ A B C ]= [',Parameters(3:5),' ];  The suitable temperature ranges is ',Parameters(1),' K to ',Parameters(2),' K');
        case 2
             fprintf(fid,'%s\n','The Saturation Vapor Pressure calculated with  WAGNER EQUATION, The parameters are: !');
             fprintf(fid,'%s %8.2f %8.2f %10.6f %10.6f %10.6f %10.6f %s %7.2f %s %7.2f %s\n','[ A B C D E F]= [',Parameters(3:8),' ];  The suitable temperature ranges is ',Parameters(1),' K to ',Parameters(2),' K');
    end
%     fprintf(fid,'%s %s %s\n','The ANTOINE EQUATION [ log10P=A-B/(T+C) ] parameters of ',comp1,' as follows:');
%     fprintf(fid,'%s %8.5f %10.3f %8.3f %s %7.2f %s %7.2f %s\n','[ A B C D E F]= [',Parameters,' ];  The suitable temperature ranges is ',TLimit(1),' K to ',TLimit(2),' K');
    fprintf(fid,'%s%d%s%f%f%s%d%s%s\n','Program contral: parameterFlag=',parameterFlag,', [alpha chb]=',alpha,chb,' ; cosmoFlag=',cosmoFlag,'    Calculated on: ',date);
    fprintf(fid,'%s %6.3f%s%d%s\n',' Relative Absolute Average Deviation of Pressure(%) is :',P_RAAD,' of ',PosDataNum,' positive data.');
    if length(T)>1
        fprintf(fid,'%7s%10s%10s%15s%15s%15s%15s%15s%15s\n','T(K)',' x1',' x2','Gamma1','Gamma2','lnGamma1','lnGamma2','Pcal(kPa)','Pexp(kPa)');
        fprintf(fid,'%s\n','============================================================================================================');
        fprintf(fid,'%7.2f%10.4f%10.4f%15.5f%15.5f%15.5f%15.5f%15.6f%15.6f\n',[T,X,gamma,lngamma,Pcal,Pexp]');
    else        
        fprintf(fid,'%7s%7s%15s%15s%15s%15s%15s%15s\n',' x1',' x2','Gamma1','Gamma2','lnGamma1','lnGamma2','Pcal(kPa)','Pexp(kPa)');
        fprintf(fid,'%s\n','================================================================================================');
        fprintf(fid,'%10.4f%10.4f%15.5f%15.5f%15.5f%15.5f%15.6f%15.6f\n',[X,gamma,lngamma,Pcal,Pexp]');
    end
    fclose(fid);
    %%  PLOT THE FIGURES OF VLE RESULTS
    plot(X(:,1),Pexp,'ko','MarkerSize',markersize);hold on;
    plot(X(:,1),Pcal,'r-','LineWidth',linewidth,'MarkerSize',markersize);
    set(gcf,'color',[1 1 1]);
    legend({'Exp.','COSMO-SAC'},'Location','Best');
    xlabel('x_1','FontSize',fontsize);
    ylabel('P (kPa)','FontSize',fontsize);
    if length(T)>1
        title(['VLE Results of ',comp1,'(1)+',comp2,'(2) at ',num2str(min(T)),'-',num2str(max(T)),'K'],'FontSize',fontsize);
    else
        title(['VLE Results of ',comp1,'(1)+',comp2,'(2) at ',num2str(T),'K'],'FontSize',fontsize);
    end
    saveas(gcf,[VLEpngPath,SystemVle{1,1}{i},'.png'],'png');hold off;
    close(gcf);
    end

    case 2
        systemGamma=importdata('systemGammaList.txt');
        Nsys=length(systemGamma.data(:,1));
        for i=1:Nsys
            tic;
            T=systemGamma.data(i,1);  %SYSTEM TEMPERATURE, K
            comp1=systemGamma.textdata{i,1};comp2=systemGamma.textdata{i,2};
            %+COMPONENT 1 PROFILE GENERATION
            [area1,volume1,prf1]=PurePrf(comp1,SolteCosmoPath,cosmoFlag);
            %+COMPONENT 2 PROFILE GENERATION
            [area2,volume2,prf2]=PurePrf(comp2,SolteCosmoPath,cosmoFlag);
            % %+INPUT MOLE FRACTION
            x1=(linspace(0.00001,0.99999,30))';
            X=[x1,1-x1];

            %COSMO AREA AND VOUME MATRICES CONSTRUCTIN
            Acosmo=[area1(1),area2(1)];
            Vcosmo=[volume1,volume2];
            %CALL gammares2007()
            [gammaR,lngammaR]=gammares2007(Vcosmo,Acosmo,X,prf1,prf2,T,parameterFlag);
            %CALL gammacomb()
            [gammaC,lngammaC]=gammacomb(Vcosmo,Acosmo,X,flag);
            % gammaR;
            % lngammaR;
            % XgammaR=X.*gammaR
            % gammaC;
            % lngammaC;
            % x1'
            lngamma=lngammaR+lngammaC;
            gamma=exp(lngamma);
            gammax=X.*gamma;
            G=X.*log(gammax);
            Gtot=G(:,1)+G(:,2);            
            %OUTPUT THE INFORMATION OF ACTIVITY COEFFICIENTS WITH TXT FORMAT
            %+PATH OF ACTIVITY COEFFICIENT RESULTS
            timeCost1=toc;
            pathGamma=['E:\MATLAB\GAMMA\',comp1,'_',comp2,'Act.txt'];
            fid=fopen(pathGamma,'w');
            fprintf(fid,'%s%s%s%s%s%10.3f%s\n\n','Results for activity coefficient predictions of the systerm:',comp1,'(1)+',comp2,'(2) at ',T,' K.');            
            fprintf(fid,'%s%d%s%d%s%s\n','Program contral: parameterFlag=',parameterFlag,' ; cosmoFlag=',cosmoFlag,'    Calculated on: ',date);            
            fprintf(fid,'%s%10.2f%s\n','Expend ',timeCost1,' seconds');
            fprintf(fid,'%15s%15s%15s%15s%15s%15s%15s%15s%15s\n','x1','x2','Gamma1','Gamma2','lnGamma1','lnGamma2','ln(Gamma1x1)','ln(Gamma2x2)','Gtot');
            fprintf(fid,'%15.6f%15.6f%15.8f%15.8f%15.8f%15.8f%15.8f%15.8f%15.8f\n',[X,gamma,lngamma,G,Gtot]');
            fclose(fid);
        end
    case 3
        xinf=0.00000001;
        x=[0.207 0.793 xinf;0.9999297 0.0000703 xinf];
%         systemOW=importdata('systemOWlist.txt');
        comp3='benzene';
        [gamma]=gammaOW(comp3,x);
        logKow=log10(gamma(2,3)/gamma(1,3))+log10(0.140214)
       
end
% END OF RECORDING TIME OF PROGRAM
timeCost=toc
