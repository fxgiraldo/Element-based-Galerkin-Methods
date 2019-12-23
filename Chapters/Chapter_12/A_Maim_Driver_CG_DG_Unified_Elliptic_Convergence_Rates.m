%---------------------------------------------------------------------%
%This code solves the 2D Poisson Equation using Unified CG/DG methods
%with tensor product of 1D basis function with either
%Exact or Inexact Integration and Using an NPOIN based data-structure
%Written by F.X. Giraldo on 9/2014
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
clear all; 
close all;

tic

%Input Data
nel=4; %Number of Elements
nop=8;    %Interpolation Order

plot_matrices=0;
plot_solution=0;
integration_type=1; %=1 is inexact and =2 is exact
space_method='dg'; %=cgc for CG continuous; 
                    %=cgd for CG discontinuous;
                    %=dgd for DG
elliptic_method='SIP';    %=LDG or SIP
bc_type='strong'; %strong=hard set; weak=weak form => For Dirichlet Only
mu_constant=1e5;    %Penalty term for SIP
alpha=1; beta=1-alpha;
icase=1; %1=2D with homogeneous BCs in x and y; 
         %2=1D with homogeneous BCs along x=-1/+1 and non-homogeneous along y=-1/+1 .
         %3=2D with non-homogeneous BCs along x and y.

nopp=[1 2 4 8 16 32];
inop_begin=1;
inop_end=4;
for inop=inop_begin:inop_end
    nop=nopp(inop);
    ngl=nop + 1;

switch nop
    case (1)
        ibeg=8;
        iskip=4;
        iend=32;
    case (2)
        ibeg=4;
        iskip=2;
        iend=20;
    case (4)
        ibeg=2;
        iskip=2;
        iend=12;
    case (8)
        ibeg=1;
        iskip=1;
        iend=3;
    case (12)
        ibeg=1;
        iskip=1;
        iend=3;
    case (16)
        ibeg=1;
        iskip=1;
        iend=2;
    case (32)
        ibeg=1;
        iskip=1;
        iend=2;
end %switch

icount=0;
for nel=ibeg:iskip:iend 
    icount=icount + 1;
    t0=cputime;
    nelx=nel;
    nely=nel;
    nelem=nelx*nely; %Number of Elements
    ngl=nop + 1;
    npts=ngl*ngl;
    npoin_CG=(nop*nelx + 1)*(nop*nely + 1);
    npoin_DG=npts*nelem;
    nboun=2*nelx + 2*nely;
    nside=2*nelem + nelx + nely;

    %Compute LGL Points
    [xgl,wgl]=legendre_gauss_lobatto(ngl);

    if (integration_type == 1)
        noq=nop;
        integration_text = ['Inexact'];
    elseif (integration_type == 2)
        noq=nop+1;
        integration_text = ['Exact'];
    end
    nq=noq + 1;
    main_text=[space_method ':'  integration_text];

    %Compute Legendre Cardinal functions and derivatives
    [psi,dpsi,xnq,wnq] = lagrange_basis(ngl,nq,xgl);

    %Create CG-Storage Grid
    [coord_CG,intma_CG,bsido_CG] = create_grid(npoin_CG,nelem,nboun,...
                                    nelx,nely,ngl,xgl);

    %Create CGDG-Storage Grid
    [coord,intma,DG_to_CG,npoin] = create_CGDG_Storage(space_method,npoin_CG,npoin_DG,coord_CG, ...
                                          intma_CG,nelem,ngl);

    %Compute Metric Terms
    [ksi_x,ksi_y,eta_x,eta_y,jac] = metrics(coord,intma,psi,dpsi,wnq,nelem,ngl,nq);

    %Compute Side/Edge Information
    [iside,jeside] = create_side(intma_CG,bsido_CG,npoin_CG,nelem,nboun,nside,ngl);
    [psideh,imapl,imapr] = create_side_dg(iside,intma_CG,nside,nelem,ngl);
    [nx,ny,jac_side] = compute_normals(psideh,intma_CG,coord_CG,...
                       nside,ngl,nq,wnq,psi,dpsi);

    %Compute Exact Solution
    [qe,qe_x,qe_y,fe]=exact_solution(coord,npoin,icase);

    %Create RHS Vector and LMatrix  
    Rvector = create_Source_Vector(jac,intma,psi,npoin,nelem,ngl,nq,fe);
    Lmatrix=zeros(npoin,npoin);
    Dmatrix_x=zeros(npoin,npoin);
    Dmatrix_y=zeros(npoin,npoin);
    Fmatrix_Q_x=zeros(npoin,npoin);
    Fmatrix_Q_y=zeros(npoin,npoin);
    Fmatrix_q_x=zeros(npoin,npoin);
    Fmatrix_q_y=zeros(npoin,npoin);
    if (elliptic_method == 'LDG')
        disp(['IN LDG LOOP ']);

        Mmatrix = create_Mmatrix(jac,intma,psi,npoin,nelem,ngl,nq);
        Mmatrix_inv=inv(Mmatrix);
        for i=1:npoin
            q=zeros(npoin,1);
            q(i)=1;
            rhs = create_Dmatrix(intma,jac,ksi_x,ksi_y,eta_x,eta_y,psi,dpsi,...
                  npoin,nelem,ngl,nq,q);
            Dmatrix_x(:,i)=rhs(:,1);
            Dmatrix_y(:,i)=rhs(:,2);
            rhs = compute_flux_LDG(q,psideh,nx,ny,jac_side,psi,npoin,...
                   nside,ngl,nq,imapl,imapr,intma,qe,alpha,beta);          
            Fmatrix_q_x(:,i)=rhs(:,1);
            Fmatrix_q_y(:,i)=rhs(:,2);
            rhs = compute_flux_LDG(q,psideh,nx,ny,jac_side,psi,npoin,...
                   nside,ngl,nq,imapl,imapr,intma,qe_x,beta,alpha);          
            Fmatrix_Q_x(:,i)=rhs(:,1);
            rhs = compute_flux_LDG(q,psideh,nx,ny,jac_side,psi,npoin,...
                   nside,ngl,nq,imapl,imapr,intma,qe_y,beta,alpha); 
            Fmatrix_Q_y(:,i)=rhs(:,2);
        end
        Dhatmatrix_q_x=Fmatrix_q_x + Dmatrix_x;
        Dhatmatrix_q_y=Fmatrix_q_y + Dmatrix_y;
        Dhatmatrix_Q_x=Fmatrix_Q_x + Dmatrix_x;
        Dhatmatrix_Q_y=Fmatrix_Q_y + Dmatrix_y;
        Lmatrix=Dhatmatrix_Q_x*Mmatrix_inv*Dhatmatrix_q_x + Dhatmatrix_Q_y*Mmatrix_inv*Dhatmatrix_q_y;
    else %CGC, CGD, or DG with SIP
        for i=1:npoin
            q=zeros(npoin,1);
            q(i)=1;
            rhs = create_Lmatrix_Vector(intma,jac,ksi_x,ksi_y,eta_x,eta_y,psi,dpsi,...
                  npoin,nelem,ngl,nq,q);
            rhs = compute_flux_SIP(rhs,q,psideh,nx,ny,jac_side,jac,psi,dpsi,...
               nside,ngl,nq,imapl,imapr,intma,ksi_x,ksi_y,eta_x,eta_y,mu_constant,qe,qe_x,qe_y);
            Lmatrix(:,i)=rhs(:);
        end
    end %if

    %Impose Homogeneous Dirichlet Boundary Conditions
    if (strcmp(bc_type,'strong'))
    %     if (space_method == 'cgc' | space_method == 'cgd') %| space_method == 'dgd')
    %        [Lmatrix,Rvector] = apply_Dirichlet_BC_Vector(Lmatrix,Rvector,psideh,...
    %                    nside,ngl,imapl,intma,qe);
    %     end
           [Lmatrix,Rvector] = apply_Dirichlet_BC_Vector(Lmatrix,Rvector,psideh,...
                       nside,ngl,imapl,intma,qe);
    end

    %Apply DSS to the Matrix and Vector 
    if (strcmp(space_method,'cgd') > 0)
       [Lmatrix,Rvector] = apply_dss_matrices_Vector(Lmatrix,Rvector,DG_to_CG,npoin,npoin_CG);
    end

    %Solve System 
    q0=Lmatrix\Rvector; 

    if (strcmp(space_method,'cgd') > 0)
       q0=apply_dss_vector(q0,DG_to_CG,npoin,npoin_CG);
    end

    %Compute Norm
    top=0;
    bot=0;
    for i=1:npoin
        top=top + (q0(i)-qe(i))^2;
        bot=bot + qe(i)^2;
    end %e
    l2_norm=sqrt( top/bot + eps );
    l2_norm_total(icount,inop)=l2_norm;
    npoin_total(icount,inop)=npoin;
    t1=cputime;
    dt=t1-t0;
    wallclock_time(icount,inop)=t1-t0;
    disp([' nop  = ' num2str(nop),' nel = ' num2str(nel),' cpu = ' num2str(dt) ]);

end %nel
end %inop

%Plot Solution
if (plot_solution == 1)
    h=figure;
    figure(h);
    xmin=min(coord(:,1)); xmax=max(coord(:,1));
    ymin=min(coord(:,2)); ymax=max(coord(:,2));
    xe=coord(:,1);
    ye=coord(:,2);
    nxx=100; nyy=100;
    dx=(xmax-xmin)/nxx;
    dy=(ymax-ymin)/nyy;
    [xi,yi]=meshgrid(xmin:dx:xmax,ymin:dy:ymax);
    qi=griddata(xe,ye,q0,xi,yi,'cubic');
    % [cl,h]=contourf(xi,yi,qi);
    surf(xi,yi,qi);
    colorbar('SouthOutside');
    xlabel('X','FontSize',18);
    ylabel('Y','FontSize',18);
    axis image
    title_text=[space_method ', Ne = ' num2str(nelem) ', N = ' num2str(nop) ', Q = ' num2str(noq) ', L2 Norm = ' num2str(l2_norm)];
    title([title_text],'FontSize',18);
    set(gca, 'FontSize', 18);
    %file_ps=['se_n' num2str(nelem) 'p' num2str(nop)];
    %eval(['print ' file_ps ' -depsc']);

    disp(['space_method = ',space_method]);
    disp(['nop = ',num2str(nop),'  nelem = ',num2str(nelem) ]);
    q_max=max(q0);
    q_min=min(q0);
    disp(['L2_Norm = ',num2str(l2_norm),' q_max = ',num2str(q_max),'  q_min = ',num2str(q_min) ]);
    disp(['npoin = ',num2str(npoin),' npoin_CG = ',num2str(npoin_CG),'  npoin_DG = ',num2str(npoin_DG) ]);
end

%Plot E-values
if (plot_matrices == 1)
    figure
    E=eig(Lmatrix);
    m=length(E);
    for i=1:m
        norm_E(i)=sqrt( conj(E(i))*E(i) );
    end
    max_norm_E=max(norm_E);
    E=E/max_norm_E;
    plot_handle=plot(real(E),imag(E),'ro');
    title_text=[' E-values with max(Re) = ', num2str(max(real(E))) ];
    title([title_text],'FontSize',18);
    set(plot_handle,'LineWidth',2);
    xlabel('Re','FontSize',18);
    ylabel('Im','FontSize',18);
    set(gca, 'FontSize', 18);
    axis([ -1 +1 -1 +1]);
   
    figure
    spy(Lmatrix);
    title_text=[' Lmatrix ' ];
    title([title_text],'FontSize',18);
    xlabel('Columns','FontSize',18);
    ylabel('Rows','FontSize',18);
    set(gca, 'FontSize', 18);
    
%     figure
%     spy(Mmatrix);
%     title_text=[' Mmatrix ' ];
%     title([title_text],'FontSize',18);
%     xlabel('Columns','FontSize',18);
%     ylabel('Rows','FontSize',18);
%     set(gca, 'FontSize', 18);
%     
%     figure
%     spy(Dmatrix_x);
%     title_text=[' Dmatrix_x ' ];
%     title([title_text],'FontSize',18);
%     xlabel('Columns','FontSize',18);
%     ylabel('Rows','FontSize',18);
%     set(gca, 'FontSize', 18);
%       
%     figure
%     spy(Dmatrix_y);
%     title_text=[' Dmatrix_y ' ];
%     title([title_text],'FontSize',18);
%     xlabel('Columns','FontSize',18);
%     ylabel('Rows','FontSize',18);
%     set(gca, 'FontSize', 18);

end

h=figure;
figure(h);
for inop=inop_begin:inop_end
    switch inop
    case (1)
        plot_handle=semilogy(npoin_total(:,inop),l2_norm_total(:,inop),'r-');
    case (2)
        plot_handle=semilogy(npoin_total(:,inop),l2_norm_total(:,inop),'b-o');
    case (3)
        plot_handle=semilogy(npoin_total(:,inop),l2_norm_total(:,inop),'g-x');
    case (4)
        plot_handle=semilogy(npoin_total(:,inop),l2_norm_total(:,inop),'k-+');
    case (5)
        plot_handle=semilogy(npoin_total(:,inop),l2_norm_total(:,inop),'m-*');
    case (6)
        plot_handle=semilogy(npoin_total(:,inop),l2_norm_total(:,inop),'c-s');
end %switch
set(plot_handle,'LineWidth',2);
hold on
end
if (noq == nop)
    title_text=[space_method ': Inexact Integration'];
elseif (noq > nop)
    title_text=[space_method ': Exact Integration']; 
end
title([title_text],'FontSize',18);
xlabel('N_p','FontSize',18);
ylabel('Normalized L^2 Error','FontSize',18);
legend('N=1','N=2','N=4','N=8','N=16');
set(gca, 'FontSize', 18);

h=figure;
figure(h);
for inop=inop_begin:inop_end
    switch inop
    case (1)
        plot_handle=semilogy(wallclock_time(:,inop),l2_norm_total(:,inop),'r-');
    case (2)
        plot_handle=semilogy(wallclock_time(:,inop),l2_norm_total(:,inop),'b-o');
    case (3)
        plot_handle=semilogy(wallclock_time(:,inop),l2_norm_total(:,inop),'g-x');
    case (4)
        plot_handle=semilogy(wallclock_time(:,inop),l2_norm_total(:,inop),'k-+');
    case (5)
        plot_handle=semilogy(wallclock_time(:,inop),l2_norm_total(:,inop),'m-*');
    case (6)
        plot_handle=semilogy(wallclock_time(:,inop),l2_norm_total(:,inop),'c-s');
end %switch
set(plot_handle,'LineWidth',2);
hold on
end
if (noq == nop)
    title_text=[space_method ': Inexact Integration'];
elseif (noq > nop)
    title_text=[space_method ': Exact Integration']; 
end
title([title_text],'FontSize',18);
xlabel('Wallclock Time (S)','FontSize',18);
ylabel('Normalized L^2 Error','FontSize',18);
legend('N=1','N=2','N=4','N=8','N=16');
set(gca, 'FontSize', 18);

toc
