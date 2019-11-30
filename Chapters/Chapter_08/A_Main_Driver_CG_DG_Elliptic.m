%---------------------------------------------------------------------%
%This code computes the 1D Diffusion Equation using the 
%CG and DG methods with 3rd Order RK.
%This version constructs the Global Matrices which are good for 
%comparing CG and DG.
%Written by F.X. Giraldo on July 3, 2012
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
clear all; 
close all;

tic

%Input Data
nelem=4; %Number of Elements == 22
nop=16;    %Interpolation Order == 5

iplot_solution=1; %Switch to Plot or Not.
iplot_matrices=0;
iplot_title=1; %=1 plot title in Figures
integration_points=1; %=1 for LGL and =2 for LG
integration_type=1; %=1 is inexact and =2 is exact
space_method_type='dg'; %=1 for CG and =2 for DG
elliptic_method_type=2; %=1 IBP with Strong BC for CG and 
                        %=2 LDG for DG
                        %=3 SIP for DG
alpha=0.5; beta=1-alpha; %LDG weights
%tau=1e3; %SIP Penalty Weight
tau=0e8;

icase=1; %case number: 1=Homogeneous Dirichlet with Sine;
                      %2=Homogeneous Dirichlet with Half-Cosine;
                      %3=Non-Homogeneous Dirichlet with Cosine

%Store Constants
ngl=nop + 1;
method_text = [space_method_type];
if space_method_type == 'cg'
    npoin=nop*nelem + 1;
    sipdg_flag=1; %=0 turns off Boundary Integrals
elseif space_method_type == 'dg'
    npoin=ngl*nelem;
    sipdg_flag=1;
end

%Compute Interpolation and Integration Points
[xgl,wgl]=legendre_gauss_lobatto(ngl);
if (integration_points == 1)
    integration_text=['LGL'];
    if (integration_type == 1)
        noq=nop;
    elseif (integration_type == 2)
        noq=nop+1;
    end
    nq=noq + 1;
    [xnq,wnq]=legendre_gauss_lobatto(nq);
elseif (integration_points == 2)
    integration_text=['LG'];
    noq=nop;
    nq=noq + 1;
    [xnq,wnq]=legendre_gauss(nq);
end
main_text=[method_text ': ' integration_text];

%Compute Lagrange Polynomial and derivatives
[psi,dpsi] = lagrange_basis3(ngl,nq,xgl,xnq);
   
%Create Grid
[coord,intma_cg]=create_grid(ngl,nelem,xgl);

%Form Global Matrix Pointer
intma=zeros(ngl,nelem);
if space_method_type == 'cg'
    intma=intma_cg;
elseif space_method_type == 'dg'
    ip=0;
    for e=1:nelem
        for i=1:ngl
            ip=ip+1;
            intma(i,e)=ip;
        end
    end
end

%Compute Exact Solution
[qe,qe_x,fe] = exact_solution(intma,coord,npoin,nelem,ngl,icase);

%Create Local/Element Mass and Differentiation Matrices
Mmatrix = create_Mmatrix(intma,coord,npoin,nelem,ngl,nq,wnq,psi);

%Compute LMATRIX using RHS approach
Lmatrix=zeros(npoin,npoin);
Dhatmatrix_Q=zeros(npoin,npoin);
Dhatmatrix_q=zeros(npoin,npoin);
RHSvector=Mmatrix*fe;
if elliptic_method_type == 1 %IBP with Strong BC
    for i=1:npoin
        q=zeros(npoin,1);
        q(i)=1;
        Lmatrix(:,i)=create_Lmatrix_IBP(intma,coord,npoin,nelem,ngl,nq,wnq,dpsi,q);
    end

    %Apply Dirichlet BC
    RHSvector(1)=qe(1);
    RHSvector(npoin)=qe(npoin);
    Lmatrix(1,:)=0;
    Lmatrix(1,1)=1;
    Lmatrix(npoin,:)=0;
    Lmatrix(npoin,npoin)=1;
    
elseif elliptic_method_type == 2 %LDG
    for i=1:npoin
        q=zeros(npoin,1);
        q(i)=1;
        Dhatmatrix_Q(:,i)=create_Dhatmatrix(intma,coord,npoin,nelem,ngl,nq,wnq,psi,dpsi,q,qe_x,beta,alpha);
    end
    Q=Dhatmatrix_Q\RHSvector;
    for i=1:npoin
        q=zeros(npoin,1);
        q(i)=1;
        Dhatmatrix_q(:,i)=create_Dhatmatrix(intma,coord,npoin,nelem,ngl,nq,wnq,psi,dpsi,q,qe,alpha,beta);
    end
    %q=Dhatmatrix_q\(Mmatrix*Q);
    Lmatrix=Dhatmatrix_Q*inv(Mmatrix)*Dhatmatrix_q;

elseif elliptic_method_type == 3 %SIPDG
    for i=1:npoin
        q=zeros(npoin,1);
        q(i)=1;
        Lmatrix(:,i)=create_Lmatrix_SIPDG(intma,coord,npoin,nelem,ngl,nq,wnq,psi,dpsi,q,qe,qe_x,tau,sipdg_flag);
    end
    
    %Apply Dirichlet BC
%     RHSvector(1)=qe(1);
%     RHSvector(npoin)=qe(npoin);
%     Lmatrix(1,:)=0;
%     Lmatrix(1,1)=1;
%     Lmatrix(npoin,:)=0;
%     Lmatrix(npoin,npoin)=1;
end

%Solve Equations
q=Lmatrix\RHSvector;  

%Compute Norm
top=0;
bot=0;
error=zeros(npoin,1);
for i=1:npoin
   top=top + (q(i)-qe(i))^2;
   error(i)=(q(i)-qe(i));
   bot=bot + qe(i)^2;
end %i
l2_norm=sqrt( top/bot )
npoin
nelem
ngl
nq

%Compute a gridpoint solution
x_sol=zeros(npoin,1);
for ie=1:nelem
   for i=1:ngl
      ip=intma(i,ie);
      x_sol(ip)=coord(i,ie);
   end 
end

%Plot Solution
if (iplot_solution == 1)
    h=figure;
    figure(h);
    plot_handle=plot(x_sol,q,'r-');
    set(plot_handle,'LineWidth',2);
    hold on
    plot_handle=plot(x_sol,qe,'b--');
    set(plot_handle,'LineWidth',2);

    xlabel('x','FontSize',18);
    ylabel('q(x,t)','FontSize',18);

    title_text=[main_text ': Ne = ' num2str(nelem) ', N = ' num2str(nop) ', Q = ' num2str(noq) ', L2 Norm = ' num2str(l2_norm)];
    title([title_text],'FontSize',18);
    set(gca, 'FontSize', 18);
    
end

if (iplot_matrices == 1)

    %Plot E-values
    figure
    E=eig(Lmatrix);
    m=length(E);
    for i=1:m
        norm_E(i)=sqrt( conj(E(i))*E(i) );
    end
    max_norm_E=max(norm_E);
    E=E/max_norm_E;
    plot_handle=plot(real(E),imag(E),'ro');
    title_text=[' Lmatrix E-values with max(Re) = ', num2str(max(real(E))) ];
    if (iplot_title == 1)
        title([title_text],'FontSize',18);
    end
    set(plot_handle,'LineWidth',2);
    xlabel('Re','FontSize',18);
    ylabel('Im','FontSize',18);
    set(gca, 'FontSize', 18);
    axis([ -1 +1 -1 +1]);
%     file_ps=['DG_n' num2str(nelem) 'p' num2str(nop) 'q' num2str(noq) '_Evalues'];
%     eval(['print ' file_ps ' -depsc']);
    
    figure
    spy(Mmatrix);
    title_text=[' Mmatrix ' ];
    if (iplot_title == 1)
        title([title_text],'FontSize',18);
    end
    xlabel('Columns','FontSize',18);
    ylabel('Rows','FontSize',18);
    set(gca, 'FontSize', 18);
%     file_ps=['DG_n' num2str(nelem) 'p' num2str(nop) 'q' num2str(noq) '_Mmatrix'];
%     eval(['print ' file_ps ' -depsc']);
    
    figure
    spy(Dhatmatrix_q);
    title_text=[' Dhatmatrix_q ' ];
    if (iplot_title == 1)
        title([title_text],'FontSize',18);
    end
    xlabel('Columns','FontSize',18);
    ylabel('Rows','FontSize',18);
    set(gca, 'FontSize', 18);
%     file_ps=['DG_n' num2str(nelem) 'p' num2str(nop) 'q' num2str(noq) '_Mmatrix'];
%     eval(['print ' file_ps ' -depsc']);

    figure
    spy(Dhatmatrix_Q);
    title_text=[' Dhatmatrix_Q ' ];
    if (iplot_title == 1)
        title([title_text],'FontSize',18);
    end
    xlabel('Columns','FontSize',18);
    ylabel('Rows','FontSize',18);
    set(gca, 'FontSize', 18);
%     file_ps=['DG_n' num2str(nelem) 'p' num2str(nop) 'q' num2str(noq) '_Mmatrix'];
%     eval(['print ' file_ps ' -depsc']);

    figure
    spy(Lmatrix);
    title_text=[' Lmatrix ' ];
    if (iplot_title == 1)
        title([title_text],'FontSize',18);
    end
    xlabel('Columns','FontSize',18);
    ylabel('Rows','FontSize',18);
    set(gca, 'FontSize', 18);
%     file_ps=['DG_n' num2str(nelem) 'p' num2str(nop) 'q' num2str(noq) '_Lmatrix'];
%     eval(['print ' file_ps ' -depsc']);
    
end
