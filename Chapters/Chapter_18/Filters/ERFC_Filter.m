%---------------------------------------------------------------------%
%This code computes the ERFC Filter Function
%Written by F.X. Giraldo on 1/2016
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [sigma] = ERF_Filter(N,alpha,F,k,s,kmax,erfc_log)

%Constants
sigma=zeros(kmax,1);
eps=1e-8;

% 
%Compute Weight
for i=1:kmax
    if (k(i) <= s)
        sigma(i)=1;
    elseif (k(i) > s)
        x=(k(i)-s)/(N-s);
        omega=abs(x) - 0.5;
        temp=sqrt( -log(1-4*omega^2)/(4*omega^2) );
        temp=2*sqrt(F)*omega*temp;
        if (erfc_log > 0)
            sigma(i)=0.5*erfc(temp);
        else
         sigma(i)=0.5*erfc(2*sqrt(F)*omega);
        end
    end
end



      
