% SSPRK linear m stage m-1 order methods 
%
% method described at 115p
% @article{gottlieb:2005high,
% 	Author = {Gottlieb, S.},
% 	Journal = {J. Sci. Comput.},
% 	Number = {1},
% 	Pages = {105-128},
% 	Title = {On high order strong stability preserving Runge--Kutta and multi step time discretizations},
% 	Volume = {25},
% 	Year = {2005}}

function alpha = compute_ti_aux(order)

m = order + 1;
mm = m+1;
%mm = m;
alpha = zeros(mm);

alpha(3,1) = 0;
alpha(3,2) = 1;

for kk = 3:m
    alpha(kk+1,kk) = 2/kk * alpha(kk,kk-1);
    for k = 2:kk-1
        alpha(kk+1,k) = 2/(k-1)*alpha(kk,k-1);
    end
    sum_temp = 0;
    for ii=2:kk
        sum_temp = sum_temp + alpha(kk+1,ii);
    end
    alpha(kk+1,1) = 1-sum_temp;
end

end