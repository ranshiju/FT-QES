function [para]=Check_para(para)
if(rem(para.L,2)==1)
    warning('para.L should be even! Change it from %g to %g. \n',para.L,para.L+1);
    para.L=para.L+1;
end
end