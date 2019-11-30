%---------------------------------------------------------------------%
%This code computes the Legendre_Gauss_Lobatto Filter Matrix
%Written by F.X. Giraldo on 4/2000
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function f = filter_init_modal(P,xmu,iplot_filter)

%Constants
p=P-1;
ph=floor( (p+1)/2 );
alpha=17;
order=18;
order=6;

%Initialize
leg=zeros(P,P);
f=zeros(P,P);

%Compute Weight
% for i=1:P
%    weight(i)=exp(  - alpha*( (i-1)/p )^order );
% end
%ibeg=round(2/3*P);
ibeg=1;
iend=P;
idelta=iend-ibeg;
%idelta=0.05/idelta;
for i=1:ibeg
    weight(i)=1;
end
for i=ibeg+1:iend
    weight(i)=1.0 - ((i-ibeg)/(iend-ibeg))^20;
end
weight

%Construct 1D Filter Matrix
for i=1:P
    f(i,i)=xmu*weight(i) + (1-xmu);
end %i	 

%Plot Solution
if (iplot_filter == 1)
    figure;
    plot_handle=plot(0:p,weight,'r-');
    set(plot_handle,'LineWidth',2);
    axis([0 p 0 1]);
    xlabel('Modes','FontSize',18);
    ylabel('Filter Weight','FontSize',18);
    set(gca, 'FontSize', 18);
end
      
