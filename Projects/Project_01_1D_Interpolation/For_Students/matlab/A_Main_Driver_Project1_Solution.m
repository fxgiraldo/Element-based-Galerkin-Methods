%---------------------------------------------------------------------%
%This code computes the Interpolation, Differentiation, amd Integration 
%of a known function using Lobatto, Legendre, Equally-Spaced, and 
%Chebyshev points.
%
%This is the solution template to Project 1 in Element-based Galerkin Methods
%
%Written by F.X. Giraldo on 7/2021
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
nop_max=64;
nopp=zeros(nop_max,1);

l1_norm_interpolation=zeros(nop_max,1);
l2_norm_interpolation=zeros(nop_max,1);
l8_norm_interpolation=zeros(nop_max,1);
l1_norm_derivative=zeros(nop_max,1);
l2_norm_derivative=zeros(nop_max,1);
l8_norm_derivative=zeros(nop_max,1);
l1_norm_integration=zeros(nop_max,1);
l2_norm_integration=zeros(nop_max,1);

for nop=1:nop_max
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
xs=linspace(-1,1,Ns);

%-------------------------------------------------------------------------%
%-----------------Students include this File------------------------------%
[psi,dpsi] = lagrange_basis_v2(ngl,Ns,xgl,xs);
%-----------------Students include this File------------------------------%
%-------------------------------------------------------------------------%


%Compute Expansion Coefficients
q_coeff=zeros(ngl,1);
for i=1:ngl
  x=xgl(i);
  q_coeff(i)=cos(c*x);
end %i

%Compute Nth Order Interpolant
qn1=psi'*q_coeff;

%Compute Exact Solution
qe1=zeros(Ns,1);
for i=1:Ns
  x=xs(i);
  qe1(i)=cos(c*x);
end %i

%Compute L1, L2, and L8 Norm
l1_norm_interpolation(inop)=( norm(qn1-qe1,1)/norm(qe1,1) );
l2_norm_interpolation(inop)=( norm(qn1-qe1,2)/norm(qe1,2) );
l8_norm_interpolation(inop)=( norm(qn1-qe1,inf)/norm(qe1,inf) );

%--------------------------------------------------%
%Derivative
%--------------------------------------------------%

%Compute Nth Order Derivative
qn1_x=dpsi'*q_coeff;

%Compute Exact Solution
qe1_x=zeros(Ns,1);
for i=1:Ns
  x=xs(i);
  qe1_x(i)=-c*sin(c*x);
end %i

%Compute L1, L2, and L8 Norm
l1_norm_derivative(inop)=( norm(qn1_x-qe1_x,1)/norm(qe1_x,1) );
l2_norm_derivative(inop)=( norm(qn1_x-qe1_x,2)/norm(qe1_x,2) );
l8_norm_derivative(inop)=( norm(qn1_x-qe1_x,inf)/norm(qe1_x,inf) );

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
%-------------------------------------------------------------------------%
%-----------------Students include the global integral in qn_init---------%
%-----------------Students include the global integral in qn_init---------%
%-------------------------------------------------------------------------%

%Compute Exact Integral
qe_int=(sin(c) - sin(-c))/c; %for c=pi/2 equals 4/pi

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


