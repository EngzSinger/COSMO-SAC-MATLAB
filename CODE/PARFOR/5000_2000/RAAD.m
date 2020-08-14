function [PositiveDataNum,RAAD]=RAAD(X,Pexp,Pcal)
% This function calcuates the Relative Absolute Average Deviation
DataNum=length(X(:,1));
PositiveDataNum=sum(X(:,1)>0);
DeltaDataNum=DataNum-PositiveDataNum;
Err=abs(Pcal((DeltaDataNum+1):DataNum)-Pexp((DeltaDataNum+1) ...
    :DataNum))./Pexp((DeltaDataNum+1):DataNum).*100.0/PositiveDataNum;
RAAD=sum(Err);
end