%---------------------------------------------------------------------%
%This code computes the Interpolation using LGL points.
%Written by F.X. Giraldo on 10/2003
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
clear all; 
close all;

%ipoints=1 are Lobatto Points
%ipoints=2 are Legendre Points
%ipoints=3 are Equi-spaced Points
%ipoints=4 are Chebyshev Points

for ipoints=1:4 %Loop to do Lobatto, Legendre, and ES points, and Chebyshev
inop=0;
for nop=1:64
inop=inop + 1;
nopp(inop)=nop;    %Interpolation Order
noq=nop;

%Store Constants
ngl=nop + 1;
nq=noq + 1;
Ns=51;
%Ns=101;
c=pi/2;

%Plot Solutions
iplot_sol=0;
iplot_sol_deriv=0;
iplot_interp=1;
iplot_interp_deriv=1;
iplot_integration=1;

%Compute LGL Points
if (ipoints == 1)
    [xgl,wgl]=legendre_gauss_lobatto(ngl);
elseif (ipoints == 2)
    [xgl,wgl]=legendre_gauss(ngl);
elseif (ipoints == 3)
    xgl=linspace(-1,1,ngl);
elseif (ipoints == 4)
    [xgl,wgl]=chebyshev_basis(ngl);
end

%--------------------------------------------------%
%Interpolation
%--------------------------------------------------%
%Compute Sample Space
xs=zeros(1,Ns);
xs=linspace(-1,1,Ns);
[psi,dpsi] = lagrange_basis_v2(ngl,Ns,xgl,xs);

%Compute Expansion Coefficients
for i=1:ngl
  x=xgl(i);
  q_coeff(i)=cos(c*x);
end %i

%Compute Nth Order Interpolant
for i=1:Ns
  qsum=0;
  for j=1:ngl
    qsum=qsum + psi(j,i)*q_coeff(j);
  end %j
  qn1(i)=qsum;
end %i
%qn1=psi'*q_coeff';

%Compute Exact Solution
for i=1:Ns
  x=xs(i);
  qe1(i)=cos(c*x);
end %i

%Compute L1, L2, and L8 Norm
l1_top=0; l1_bot=0;
l2_top=0; l2_bot=0;
l8_top=-1000; l8_bot=-1000;
for i=1:Ns
   l1_top=l1_top + abs(qn1(i)-qe1(i));
   l1_bot=l1_bot + abs(qe1(i));
   l2_top=l2_top + (qn1(i)-qe1(i))^2;
   l2_bot=l2_bot + qe1(i)^2;
   l8_top=max(l8_top, abs(qn1(i)-qe1(i)));
   l8_bot=max(l8_bot, abs(qe1(i)));
end
l1_norm_interpolation(inop)=( l1_top/l1_bot );
l2_norm_interpolation(inop)=sqrt( l2_top/l2_bot );
l8_norm_interpolation(inop)=( l8_top/l8_bot );

%--------------------------------------------------%
%Derivative
%--------------------------------------------------%

%Compute Nth Order Derivative
for i=1:Ns
  qsum=0;
  for j=1:ngl
    qsum=qsum + dpsi(j,i)*q_coeff(j);
  end %j
  qn1_x(i)=qsum;
end %i
%qn1_x=dpsi'*q_coeff';

%Compute Exact Solution
for i=1:Ns
  x=xs(i);
  qe1_x(i)=-c*sin(c*x);
end %i

%Compute L1, L2, and L8 Norm
l1_top=0; l1_bot=0;
l2_top=0; l2_bot=0;
l8_top=-1000; l8_bot=-1000;
for i=1:Ns
   l1_top=l1_top + abs(qn1_x(i)-qe1_x(i));
   l1_bot=l1_bot + abs(qe1_x(i));
   l2_top=l2_top + (qn1_x(i)-qe1_x(i))^2;
   l2_bot=l2_bot + qe1_x(i)^2;
   l8_top=max(l8_top, abs(qn1_x(i)-qe1_x(i)));
   l8_bot=max(l8_bot, abs(qe1_x(i)));
end
l1_norm_derivative(inop)=( l1_top/l1_bot );
l2_norm_derivative(inop)=sqrt( l2_top/l2_bot );
l8_norm_derivative(inop)=( l8_top/l8_bot );

%--------------------------------------------------%
%Integration
%--------------------------------------------------%

%Get Quadrature roots
weight=ones(nq);
if ipoints == 1
    [xnq,wnq] = legendre_gauss_lobatto(nq);
elseif ipoints == 2
    [xnq,wnq] = legendre_gauss(nq);
elseif ipoints == 3
    [xnq,wnq] = legendre_gauss(nq);
elseif ipoints == 4
%     [xnq,wnq] = chebyshev_gauss(nq);
%     for i=1:nq
%         weight(i)=1.0/( sqrt(1.0 - xnq(i)*xnq(i)) );
%     end
    [xnq,wnq] = legendre_gauss(nq);
end 

%Compute Basis Functions
[psi,dpsi] = lagrange_basis_v2(ngl,nq,xgl,xnq);

%Compute Integral
qn=0;
qn_int=0;
for i=1:nq
  wq=wnq(i)*weight(i);
  qsum=0;
  for j=1:ngl
    qsum=qsum + psi(j,i)*q_coeff(j);
  end %j
  %qn(i)=qsum;
  qn_int=qn_int + wq*qsum;
end %i

%Compute Exact Integral
qe_int=4/pi;
%qe_int=(sin(c) - sin(-c))/c;

l1_norm_integration(inop)=abs( ( qn_int - qe_int )/qe_int );
l2_norm_integration(inop)=sqrt( ( qn_int - qe_int )^2/qe_int^2 );

end %nop

%Plot Interpolation
if iplot_sol == 1
    h=figure;
    figure(h);
    plot_handle=plot(xs,qn1,'r-');
    set(plot_handle,'LineWidth',2);
    hold on
    plot_handle=plot(xs,qe1,'b--');
    set(plot_handle,'LineWidth',2);
    title_text=[' N = ' num2str(nop) ', Interpolation L2 Norm = ' ...
               num2str(l2_norm_interpolation(inop))];

    title([title_text],'FontSize',14);
    xlabel('x','FontSize',18);
    ylabel('q','FontSize',18);
    legend('Numerical','Exact');
    set(gca, 'FontSize', 18);
end

%Plot Derivative of Solution
if iplot_sol_deriv == 1
    h=figure;
    figure(h);
    plot_handle=plot(xs,qn1_x,'r-');
    set(plot_handle,'LineWidth',2);
    hold on
    plot_handle=plot(xs,qe1_x,'b--');
    set(plot_handle,'LineWidth',2);
    title_text=[' N = ' num2str(nop) ', Derivative L2 Norm = ' ...
               num2str(l2_norm_derivative(inop))];

    title([title_text],'FontSize',14);
    xlabel('x','FontSize',18);
    ylabel('q_x','FontSize',18);
    legend('Numerical','Exact');
    set(gca, 'FontSize', 18);
end

%Plot Interpolant Error
    if iplot_interp == 1
    h=figure;
    figure(h);
    plot_handle=semilogy(nopp,l1_norm_interpolation,'r-');
    set(plot_handle,'LineWidth',2);
    hold on;
    plot_handle=semilogy(nopp,l2_norm_interpolation,'b--');
    set(plot_handle,'LineWidth',2);
    hold on;
    plot_handle=semilogy(nopp,l8_norm_interpolation,'k-.');
    set(plot_handle,'LineWidth',2);
    if ipoints == 1
        title_text=[' Interpolation Error for Legendre-Gauss-Lobatto ' ];
    elseif ipoints == 2
        title_text=[' Interpolation Error for Legendre-Gauss '];
    elseif ipoints == 3
        title_text=[' Interpolation Error for Equally-Spaced '];
    elseif ipoints == 4
        title_text=[' Interpolation Error for Chebyshev '];
    end
    title([title_text],'FontSize',14);
    xlabel('N','FontSize',18);
    ylabel('Error Norm','FontSize',18);
    legend('L1','L2','L\infty');
    set(gca, 'FontSize', 18);
end

%Plot Derivative Error
if iplot_interp_deriv == 1
    h=figure;
    figure(h);
    plot_handle=semilogy(nopp,l1_norm_derivative,'r-');
    set(plot_handle,'LineWidth',2);
    hold on;
    plot_handle=semilogy(nopp,l2_norm_derivative,'b--');
    set(plot_handle,'LineWidth',2);
    hold on;
    plot_handle=semilogy(nopp,l8_norm_derivative,'k-.');
    if ipoints == 1
        title_text=[' Derivative Error for Legendre-Gauss-Lobatto ' ];
    elseif ipoints == 2
        title_text=[' Derivative Error for Legendre-Gauss '];
    elseif ipoints == 3
        title_text=[' Derivative Error for Equally-Spaced '];
    elseif ipoints == 4
        title_text=[' Derivative Error for Chebyshev '];
    end
    title([title_text],'FontSize',14);
    xlabel('N','FontSize',18);
    ylabel('Error Norm','FontSize',18);
    legend('L1','L2','L\infty');
    set(gca, 'FontSize', 18);
end

%Plot Integration Error
if iplot_integration == 1
    h=figure;
    figure(h);
    plot_handle=semilogy(nopp,l1_norm_integration,'r-');
    set(plot_handle,'LineWidth',2);
    hold on;
    plot_handle=semilogy(nopp,l2_norm_integration,'b--');
    set(plot_handle,'LineWidth',2);
    hold on;
    % plot_handle=semilogy(nopp,l8_norm_derivative,'k-.');
    % set(plot_handle,'LineWidth',2);
    if ipoints == 1
        title_text=[' Integration Error for Legendre-Gauss-Lobatto ' ];
    elseif ipoints == 2
        title_text=[' Integration Error for Legendre-Gauss '];
    elseif ipoints == 3
        title_text=[' Integration Error for Equally-Spaced '];
    elseif ipoints == 4
        title_text=[' Integration Error for Chebyshev '];
    end
    title([title_text],'FontSize',14);
    xlabel('N','FontSize',18);
    ylabel('Error Norm','FontSize',18);
    legend('L1','L2');
    set(gca, 'FontSize', 18);
end

end %ipoints

%file_ps=['se_n' num2str(nelem) 'p' num2str(nop)];
%eval(['print ' file_ps ' -depsc']);


