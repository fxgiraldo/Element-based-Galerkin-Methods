%---------------------------------------------------------------------%
%This code computes the Interpolation, Differentiation, amd Integration 
%of a known function using Lobatto (LGL), Legendre (LG), and Laguerre (LGR)
%points.
%
%This is the solution to Project 1 in Element-based Galerkin Methods
%The algorithms here are used in Secs. 3.4 (Algorithms 3.1 and 3.2) and 
%4.4 (Eq. 4.13) in F.X. Giraldo's Introduction to Element-based Galerkin 
%Methods using Tensor-Product Bases: Analysis, Algorithms, and Applications.
%
%Written by F.X. Giraldo on 12/2023
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
clear all; 
close all;
addpath('fxg')

%Define plotting colors
blue = [114 147 203]./255;
red = [211 94 96]./255;
black = [128 133 133]./255;
green = [132 186 91]./255;
brown = [171 104 87]./255;
purple = [144 103 167]./255;

%ipoints=1 are Lobatto Points
%ipoints=2 are Legendre Points
%ipoints=3 are Equi-spaced Points
%ipoints=4 are Chebyshev Points
%ipoints=5 are Laguerre-Gauss-Radau Points

%----------------------------INPUT PARAMETERS---------------------------%

%Plot Solutions
%f(x)
iplot_sol=0;
iplot_interp_convergence=1;
iprint_interp_error=0;
%dfdx(x);
iplot_sol_deriv=0;
iplot_interp_deriv_convergence=1;
iprint_interp_deriv_error=0;
%Integral f(x) dx
iplot_integration_convergence=1;
iprint_integration_error=1;
%Plot Basis Functions Info
iplot_basis_functions=0;
iprint_basis_functions=0;
iplot_Lebesgue_function=0;
iprint_Lebesgue_constant=0;
%Sample Points
Ns=101;
%Max Polynomial Order
nop_max=64;
%poly1=p1*nop + a1
p1=1; a1=0; %for Function Interpolation x^N
%poly2=p2*nop + a2
p2=1; a2=0; %for Function Derivative Approximation to evaluate x^{N} (e^{-x})
%p2=0; a2=0; %for Function Derivative Approximation to evaluate e^{-x}
%----------------------------INPUT PARAMETERS---------------------------%

for ipoints=5:5 %Loop to do (1) Lobatto, (2) Legendre, (3) Equi-spaced, (4) Chebyshev, and (5) Laguerre-Gauss-Radau
inop=0;

%For Integration
p3=2; a3=0;
if ipoints == 1 
    p3=2; a3=-1; %Lobatto = 2N-1
elseif ipoints == 2 
    p3=2; a3=1; %Legendre = 2N+1
elseif ipoints == 5
    p3=2; a3=0; %Laguerre = 2N
end

nopp=zeros(nop_max,1);
l1_norm_interpolation=zeros(nop_max,1);
l2_norm_interpolation=zeros(nop_max,1);
l8_norm_interpolation=zeros(nop_max,1);
l1_norm_derivative=zeros(nop_max,1);
l2_norm_derivative=zeros(nop_max,1);
l8_norm_derivative=zeros(nop_max,1);
l1_norm_integration=zeros(nop_max,1);
l2_norm_integration=zeros(nop_max,1);
Lebesgue_function=zeros(Ns,nop_max);
Lebesgue_constant=zeros(nop_max,1);

for nop=1:nop_max
inop=inop + 1;
nopp(inop)=nop;    %Interpolation Order
e1=p1*nop + a1; c1=0; %for function and derivative interpolation
e2=p2*nop + a2; c2=0; %for function and derivative interpolation
%e2=3;
e3=p3*nop + a3; c3=0.5; %for function integration (if c3=0, the integral is zero for LGL and LG)

%Store Constants
ngl=nop + 1;

%Compute LGL Points
if (ipoints == 3)
    xgl=linspace(-1,1,ngl);
    method_text=['Equi-Spaced'];
else
    [xgl,wgl,method_text] = Gauss_Quadrature_Eigenvalue(ipoints,ngl);
end

%--------------------------------------------------%
%Interpolation e=N for f(x)=x^N
%--------------------------------------------------%
%Compute Sample Space
xs_min=min(xgl); xs_max=max(xgl);
xs=linspace(xs_min,xs_max,Ns);
% xs=xgl;
% Ns=ngl;
if (ipoints <= 4)
    [psi,dpsi] = lagrange_basis_v2(ngl,Ns,xgl,xs);
elseif (ipoints == 5)
    [psi,dpsi] = lagrange_basis_laguerre_v2(ngl,Ns,xgl,xs);
end

%Compute Expansion Coefficients
q_coeff1=zeros(ngl,1);
if (ipoints <= 4)
    for i=1:ngl
      x=xgl(i);
      q_coeff1(i)=x^e1;
    end %i
elseif (ipoints == 5)
    for i=1:ngl
      x=xgl(i);
      q_coeff1(i)=x^e1*exp(-0.5*x);
    end %i
end

%Compute Nth Order Interpolant
qn=psi'*q_coeff1;

%Compute Exact Solution
qe=zeros(Ns,1);
if (ipoints <= 4)
    for i=1:Ns
        x=xs(i);
        qe(i)=x^e1;
    end %i
elseif (ipoints == 5)
    for i=1:Ns
        x=xs(i);
        qe(i)=x^e1*exp(-0.5*x);
    end %i
end

%Compute L1, L2, and L8 Norm
l1_norm_interpolation(inop)=( norm(qn-qe,1)/norm(qe,1) );
l2_norm_interpolation(inop)=( norm(qn-qe,2)/norm(qe,2) );
l8_norm_interpolation(inop)=( norm(qn-qe,inf)/norm(qe,inf) );
if (iprint_interp_error==1)
   disp([' nop = ',num2str(nop),' norm(qn-qe)/norm(qe) = ',num2str(norm(qn-qe,2)/norm(qe,2))])
end
%Rescale order because if it is exactly 0 does not show up.
l1_norm_interpolation(inop)=max(l1_norm_interpolation(inop),1e-16);
l2_norm_interpolation(inop)=max(l2_norm_interpolation(inop),1e-16);
l8_norm_interpolation(inop)=max(l8_norm_interpolation(inop),1e-16);

%--------------------------------------------------%
%Derivative for d/dx( x^N )
%--------------------------------------------------%

%Derivative Coefficients
q_coeff2=zeros(ngl,1);
if (ipoints <= 4)
    for i=1:ngl
      x=xgl(i);
      q_coeff2(i)=x^e2;
    end %i
elseif (ipoints == 5)
    for i=1:ngl
      x=xgl(i);
      q_coeff2(i)=x^e2*exp(-0.5*x);
    end %i
end

%Compute Nth Order Derivative
qn_x=dpsi'*q_coeff2;

%Compute Exact Solution
qe_x=zeros(Ns,1);
if (ipoints <= 4)
    for i=1:Ns
        x=xs(i);
        qe_x(i)=e2*x^(e2-1);
    end %i
elseif (ipoints == 5)
    for i=1:Ns
        x=xs(i);
        if (e2==0)
            qe_x(i)=-exp(-x);
        else
            qe_x(i)=exp(-0.5*x)*x^(e2-1)*(e2-0.5*x);
        end
    end %i
end

%Compute L1, L2, and L8 Norm
l1_norm_derivative(inop)=( norm(qn_x-qe_x,1)/norm(qe_x,1) );
l2_norm_derivative(inop)=( norm(qn_x-qe_x,2)/norm(qe_x,2) );
l8_norm_derivative(inop)=( norm(qn_x-qe_x,inf)/norm(qe_x,inf) );
if (iprint_interp_deriv_error==1)
    %disp([' nop = ',num2str(nop),' norm(qn_x-qe_x)= ',num2str(norm(qn_x-qe_x,2)),' norm(qe_x) = ',num2str(norm(qe_x,2)),' norm(qn_x-qe_x)/norm(qe_x) = ',num2str(norm(qn_x-qe_x,2)/norm(qe_x,2))])
    disp([' nop = ',num2str(nop),' norm(qn_x-qe_x)/norm(qe_x) = ',num2str(norm(qn_x-qe_x,2)/norm(qe_x,2))])
end

%--------------------------------------------------%
%Integration Integral x^N e^{-x} dx = Gamma(N+1) for N=2*nop-2
%--------------------------------------------------%

%Get Quadrature roots
noq=nop;
nq=noq+1;
weight=ones(nq);
if (ipoints == 3)
    [xnq,wnq] = legendre_gauss(nq);
else
    %[xnq,wnq,~] = Gauss_Quadrature_Eigenvalue(ipoints,nq);
    xnq=xgl; wnq=wgl;
end

%Compute Basis Functions
if (ipoints <= 4)
    [psi2,dpsi2] = lagrange_basis_v2(ngl,nq,xgl,xnq);
elseif (ipoints == 5)
    [psi2,dpsi2] = lagrange_basis_laguerre_v2(ngl,nq,xgl,xnq);
end

%Compute Expansion Coefficients
q_coeff3=zeros(ngl,1);
if (ipoints <= 4)
    for i=1:ngl
      x=xgl(i);
      q_coeff3(i)=x^e3 + c3;
    end %i
elseif (ipoints == 5)
    for i=1:ngl
      x=xgl(i);
      q_coeff3(i)=x^e3*exp(-x);
    end %i
end

%Compute Integral
qsum=psi2'*q_coeff3;
qn_int=0;
for i=1:nq
  qn_int=qn_int + wnq(i)*weight(i)*qsum(i);
end %i

%Compute Exact Integral
if (ipoints <= 4)
    xs_min=-1; xs_max=+1;
    qe_int=1.0/(e3+1)*(xs_max^(e3+1)-xs_min^(e3+1)) + c3*(xs_max-xs_min);
elseif (ipoints == 5)
    %qe_int=gamma(e3+1);
    qe_int=factorial(e3);
end

l1_norm_integration(inop)=norm(qn_int - qe_int,1)/norm(qe_int,1);
l2_norm_integration(inop)=norm(qn_int - qe_int,2)/norm(qe_int,2);

if (iprint_integration_error==1)
   %disp([' nop = ',num2str(nop),' norm(qe_int) ',num2str(norm(qe_int,2))])
   disp([' nop = ',num2str(nop),' norm(qn_int - qe_int)/norm(qe_int) ',num2str(l2_norm_integration(inop))])
end

%--------------------------------------------------%
%Plot Basis Functions and Lebesgue Functions
%--------------------------------------------------%
if (iprint_basis_functions==1)
    disp([' nop= ',num2str(nop),' sum(psi)/ngl= ',num2str(sum(psi(:)/ngl)),' sum(dpsi)/ngl= ',num2str(sum(dpsi(:)/ngl))])
end
for i=1:Ns
    Lebesgue_function(i,nop)=Lebesgue_function(i,nop) + sum(abs(psi(:,i)));
end
Lebesgue_constant(nop)=max(Lebesgue_function(:,nop));
if (iprint_Lebesgue_constant==1)
    disp([' nop = ',num2str(nop),' Lebesgue Constant ',num2str(Lebesgue_constant(nop))])
end

end %nop

%Plot Interpolation
if iplot_sol == 1
    h=figure;
    figure(h);
    plot_handle=plot(xs,qn,'r-');
    set(plot_handle,'LineWidth',2);
    hold on
    plot_handle=plot(xs,qe,'b--');
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
    plot_handle=plot(xs,qn_x,'r-');
    set(plot_handle,'LineWidth',2);
    hold on
    plot_handle=plot(xs,qe_x,'b--');
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
if iplot_interp_convergence == 1
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
    title_text=[' Interpolation Error for ' method_text];
    
    title([title_text],'FontSize',14);
    xlabel('N','FontSize',18);
    ylabel('Error Norm','FontSize',18);
    legend('L1','L2','L\infty');
    set(gca, 'FontSize', 18);
end %plot_interp_convergence

%Plot Derivative Error
if iplot_interp_deriv_convergence == 1
    h=figure;
    figure(h);
    plot_handle=semilogy(nopp,l1_norm_derivative,'r-');
    set(plot_handle,'LineWidth',2);
    hold on;
    plot_handle=semilogy(nopp,l2_norm_derivative,'b--');
    set(plot_handle,'LineWidth',2);
    hold on;
    plot_handle=semilogy(nopp,l8_norm_derivative,'k-.');
    title_text=[' Derivative Error for ' method_text];
    
    title([title_text],'FontSize',14);
    xlabel('N','FontSize',18);
    ylabel('Error Norm','FontSize',18);
    legend('L1','L2','L\infty');
    set(gca, 'FontSize', 18);
end %iplot_interp_deriv_convergence

%Plot Integration Error
if iplot_integration_convergence == 1
    h=figure;
    figure(h);
    plot_handle=semilogy(nopp,l1_norm_integration,'r-');
    set(plot_handle,'LineWidth',2);
    hold on;
    plot_handle=semilogy(nopp,l2_norm_integration,'b--');
    set(plot_handle,'LineWidth',2);
    hold on;
    title_text=[' Integration Error for ' method_text];
    
    title([title_text],'FontSize',14);
    xlabel('N','FontSize',18);
    ylabel('Error Norm','FontSize',18);
    legend('L1','L2');
    %legend('L1');
    set(gca, 'FontSize', 18);
end %iplot_integration_convergence

%Plot Basis Funxtions
if iplot_basis_functions == 1
    h=figure;
    figure(h);
    plot_handle0=plot(xs,psi(1,:),'r-','LineWidth',2);
    hold on;
    plot_handle0=plot(xs,psi(2,:),'b-','LineWidth',2);
    plot_handle0=plot(xs,psi(3,:),'g-','LineWidth',2);
    plot_handle0=plot(xs,psi(4,:),'k-','LineWidth',2);
    plot_handle0=plot(xs,psi(5,:),'m-','LineWidth',2);
    plot_handle0=plot(xs,psi(6,:),'c-','LineWidth',2);
    plot_handle0=plot(xs,psi(7,:),'y-','LineWidth',2);
    plot_handle0=plot(xs,psi(8,:),'Color',[0.6706 0.4078 0.3412],'LineWidth',2); %brown
    plot_handle0=plot(xs,psi(9,:),'Color',[0.5176 0.7294 0.3569],'LineWidth',2); %purple
    hold on;
    plot_handle2=plot(xs,Lebesgue_function(:,nop),'--k','LineWidth',2);
    xlabel('X','FontSize',18);
    ylabel('\psi(x)', 'FontSize', 18);
    set(gca, 'FontSize', 18);
    %set(gca, 'XTick', [-1:0.5:1.0]);
    title_text=[' Basis for ' method_text];
    
    Lebesgue_text=[' \Lambda_N = ',num2str(Lebesgue_constant(nop))];
    main_text=[title_text  ':' Lebesgue_text];
    title([main_text],'FontSize',14);
    set(gca, 'FontSize', 18);
end %iplot_basis_functions

end %for ipoints=1:5

%file_ps=['se_n' num2str(nelem) 'p' num2str(nop)];
%eval(['print ' file_ps ' -depsc']);


