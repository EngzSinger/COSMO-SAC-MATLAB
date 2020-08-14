%function [x,fval,exitflag,output] = CosmoParaOptim(x0)
%% This is an auto generated MATLAB file from Optimization Tool.

%% Start with the default options
options = optimset;
%% Modify options setting
options = optimset(options,'Display', 'iter','TolFun',1e-3,'TolX',1e-3);
%x0=[241.87,769.94];
x0=[9000,4000];
%% Recod time
tic;
disp('Optimizing! PLEASE WAIT>>>>>>>>>>>>>>>');
[x,fval,exitflag,output] = ...
fminsearch(@ObjFunc,x0,options)
TimeCost=toc
