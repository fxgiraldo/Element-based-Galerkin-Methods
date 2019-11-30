function [alpha,beta] = construct_WSRK_coefficients(stages)

alpha=zeros(stages,stages);
beta=zeros(stages,1);

for i=1:stages-1
    d=stages-i+1;
    alpha(i+1,i)=1.0/d;
end
beta(stages)=1;


end

