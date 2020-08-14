function [str2]=clearNum2UL(str1,modFlag)
% str1 is a string yu want to treat, e.g. str1='ab13D'
% modFlag what behaviour do you want to do:
%       modFlag=1: UPPER the characters; 
%       modFlag=2: LOWER the characters; 
% EXAMPLE:
% str1=str1='ab13D';
% [str2]=clearNum2UL(str1,1)
% str2=
%           ABD
% [str2]=clearNum2UL(str1,2)
% str2=
%           abd
k=length(str1);
str2='';
    for m=1:k
        c=str1(m);
        f=str2num(c);
        if isempty(f)
            str2=strcat(str2,c);
        end
    end
    switch modFlag
        case 1
            str2=upper(str2);
        case 2
            str2=lower(str2);
    end
end