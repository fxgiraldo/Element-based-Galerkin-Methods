%---------------------------------------------------------------------%
%This code computes the Interpolation, Differentiation, amd Integration 
%of a known function using Lobatto, Legendre, Equally-Spaced, and 
%Chebyshev points.
%
%This is the solution to Project 1 in Element-based Galerkin Methods
%The algorithms here are used in Secs. 3.4 (Algorithms 3.1 and 3.2) and 
%4.4 (Eq. 4.13) in F.X. Giraldo's Introduction to Element-based Galerkin 
%Methods using Tensor-Product Bases: Analysis, Algorithms, and Applications.
%
%Written by F.X. Giraldo on 7/2021
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
clear all; 
close all;
addpath('lgr') 
addpath('fxg')

%ipoints=1 are Lobatto Points
%ipoints=2 are Legendre Points
%ipoints=3 are Equi-spaced Points
%ipoints=4 are Chebyshev Points
%ipoints=5 are Laguerre-Gauss-Radau Points

%Plot Solutions
iplot_sol=1;
iplot_sol_deriv=1;
iplot_interp=1;
iplot_interp_deriv=1;
iplot_integration=1;

for ipoints=4:4 %Loop to do Lobatto, Legendre, ES points, Chebyshev, and Laguerre-Gauss-Radau
inop=0;
nop_max=64;
nopp=zeros(nop_max,1);
if ipoints == 5
    [xgl,wgl]=laguerre_gauss_radau_eigenvalue(nop_max+1,true);
    xss_max=max(xgl);
end

l1_norm_interpolation=zeros(nop_max,1);
l2_norm_interpolation=zeros(nop_max,1);
l8_norm_interpolation=zeros(nop_max,1);
l1_norm_derivative=zeros(nop_max,1);
l2_norm_derivative=zeros(nop_max,1);
l8_norm_derivative=zeros(nop_max,1);
l1_norm_integration=zeros(nop_max,1);
l2_norm_integration=zeros(nop_max,1);

for nop=1:nop_max
%for nop=3:3
inop=inop + 1;
nopp(inop)=nop;    %Interpolation Order
noq=nop;

%Store Constants
ngl=nop + 1;
nq=noq + 1;
%Ns=10;
Ns=100;
if (ipoints <= 4)
    c=pi/2;
elseif (ipoints == 5)
    d=1; %doesn't matter if it is 1,10,100
    c=0.125;
end


%Compute LGL Points
if (ipoints == 1)
    %[xgl,wgl]=legendre_gauss_lobatto(ngl);
    [xgl,wgl]=legendre_gauss_lobatto_eigenvalue(ngl);
elseif (ipoints == 2)
    [xgl,wgl]=legendre_gauss(ngl);
elseif (ipoints == 3)
    xgl=linspace(-1,1,ngl);
elseif (ipoints == 4)
    [xgl,wgl]=chebyshev_gauss(ngl);
elseif (ipoints == 5)
    [xgl,wgl]=laguerre_gauss_radau_eigenvalue(ngl,true);
end

%--------------------------------------------------%
%Interpolation
%--------------------------------------------------%
%Compute Sample Space
if ipoints <= 4
    xs_min=-1; xs_max=1; 
elseif ipoints == 5
    xs_min=0; xs_max=xss_max;
end
xs_min=min(xgl);
xs_max=max(xgl);
xs=linspace(xs_min,xs_max,Ns);
% xs=xgl;
% Ns=ngl;
if (ipoints <= 4)
    [psi,dpsi] = lagrange_basis_v2(ngl,Ns,xgl,xs);
elseif (ipoints == 5)
    [psi,dpsi] = lagrange_basis_laguerre_v2(ngl,Ns,xgl,xs);
    %[psi] = lagrange_basis_laguerre_v3(ngl,Ns,xgl,xs);
end
disp([' nop= ',num2str(nop),' sum(psi)= ',num2str(sum(psi(:))),' sum(dpsi)= ',num2str(sum(dpsi(:)))])

%Compute Expansion Coefficients
q_coeff=zeros(ngl,1);
if (ipoints <= 4)
    for i=1:ngl
      x=xgl(i);
      q_coeff(i)=cos(c*x);
    end %i
elseif (ipoints == 5)
    for i=1:ngl
      x=xgl(i);
      q_coeff(i)=d*exp(-c*x);%*sin(pi*x);
    end %i
end

%Compute Nth Order Interpolant
qn1=psi'*q_coeff;

%Compute Exact Solution
qe1=zeros(Ns,1);
for i=1:Ns
  x=xs(i);
  qe1(i)=cos(c*x);
end %i

if (ipoints <= 4)
    for i=1:Ns
        x=xs(i);
        qe1(i)=cos(c*x);
    end %i
elseif (ipoints == 5)
    for i=1:Ns
        x=xs(i);
        qe1(i)=d*exp(-c*x);%*sin(pi*x);
    end %i
end

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
if (ipoints <= 4)
    for i=1:Ns
        x=xs(i);
        qe1_x(i)=-c*sin(c*x);
    end %i
elseif (ipoints == 5)
    for i=1:Ns
        x=xs(i);
        qe1_x(i)=-c*d*exp(-c*x);
    end %i
end

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
elseif (ipoints == 5)
    [xnq,wnq]=laguerre_gauss_radau_eigenvalue(nq,true);
end 

%Compute Basis Functions
if (ipoints <= 4)
    [psi2,dpsi2] = lagrange_basis_v2(ngl,nq,xgl,xnq);
elseif (ipoints == 5)
    [psi2,dpsi2] = lagrange_basis_laguerre_v2(ngl,nq,xgl,xnq);
end

%Compute Integral
qsum=psi2'*q_coeff;
qn_int=0;
for i=1:nq
  qn_int=qn_int + wnq(i)*weight(i)*qsum(i);
end %i

%Compute Exact Integral
if (ipoints <= 4)
    qe_int=(sin(c) - sin(-c))/c; %for c=pi/2 equals 4/pi
elseif (ipoints == 5)
    qe_int=-d/c*(exp(-c*xs_max) - exp(-c*xs_min));
end

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
    elseif ipoints == 5
        title_text=[' Interpolation Error for Laguerre-Gauss-Radau '];
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
    elseif ipoints == 5
        title_text=[' Derivative Error for Laguerre-Gauss-Radau '];
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
    elseif ipoints == 5
        title_text=[' Integration Error for Laguerre-Gauss-Radau '];
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


