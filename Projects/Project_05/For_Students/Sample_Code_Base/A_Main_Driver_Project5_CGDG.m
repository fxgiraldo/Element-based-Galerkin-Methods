%---------------------------------------------------------------------%
%This code computes the 2D Advection Equation using the CG/DG methods
%with 3rd Order RK and tensor product of 1D basis function with 
%Inexact Integration and Using an NPOIN based data-structure
%Written by F.X. Giraldo on 5/2021
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

%--------------------------------%
kstages=4;  %4=RK4 DO NOT CHANGE
%--------------------------------%
dt=1; %time-step, Changes automatically to keep Courant_max fixed!
Courant_max=0.5;
time_final=0.25; %final time in revolutions
nplots=4; %Number of Frames in movie
plot_movie=0;
plot_solution=1;
store_movie=0;
space_method='cg'; %=cg for CG or =dg for DG
icase=1; %case number: 1 is a Gaussian in CW;
         %2 is a Gaussian along x;
         %3 is a Gaussian along y;
         %4 is a Gaussian along diagonal
         %5 is a Square in CW:
         %6 is a Square along x
xmu=0.05; %filtering strength: 1 is full strength and 0 is no filter
ifilter=0; %time-step frequency that the filter is applied.
flux_method='rusanov'; %rusanov or centered

%Store Constants
if (icase ==1)
    c=2*pi;
elseif (icase ==2)
    c=1;
elseif (icase ==3)
    c=1;
elseif (icase ==4)
    c=1;
elseif (icase ==5)
    c=2*pi;
elseif (icase ==6)
    c=1;
end
%ntime=time_final/dt;
dt=dt*c;
time_final=time_final*c;
nelx=nel;
nely=nel;
nelem=nelx*nely; %Number of Elements
ngl=nop + 1;
npts=ngl*ngl;
npoin_CG=(nop*nelx + 1)*(nop*nely + 1);
npoin_DG=npts*nelem;

nboun=2*nelx + 2*nely;
nface=2*nelem + nelx + nely;
ntime=time_final/dt;

%Compute LGL Points
[xgl,wgl]=legendre_gauss_lobatto(ngl);
noq=nop;
nq=noq + 1;
main_text=[space_method];

%Compute Legendre Cardinal functions and derivatives
[psi,dpsi,xnq,wnq] = lagrange_basis(ngl,nq,xgl);

%Compute Filter Matrix
f = filter_init(ngl,xgl,xmu);

%Create CG-Storage Grid
[coord_CG,intma_CG,bsido_CG,iperiodic_CG] = create_grid_2d(npoin_CG,nelem,nboun, ...
                                nelx,nely,ngl,xgl);
                         
%Create CGDG-Storage Grid
[coord,intma,iperiodic,DG_to_CG,npoin] = create_CGDG_Storage(space_method,npoin_CG,npoin_DG,coord_CG, ...
                                      intma_CG,iperiodic_CG,nelem,ngl);
                       
%Compute Metric Terms
[ksi_x,ksi_y,eta_x,eta_y,jac] = metrics(coord_CG,intma_CG,psi,dpsi,nelem,ngl,nq);

%Compute Face/Side/Edge Information
[iside,jeside] = create_side(intma_CG,bsido_CG,npoin_CG,nelem,nboun,nface,ngl);
[face,mapL,mapR] = create_face(iside,intma_CG,nface,ngl);
[nx,ny,jac_face] = compute_normals(face,intma_CG,coord_CG,...
                   nface,ngl,nq,psi,dpsi); 
%Impose Face Periodicity
face=create_face_periodicity(iside,face,coord_CG,nface,nboun);

%Compute Exact Solution
time=0;
[qe,ue,ve]=exact_solution(coord,npoin,time,icase);

%Plot Exact Solution
xmin=min(coord_CG(:,1)); xmax=max(coord_CG(:,1));
ymin=min(coord_CG(:,2)); ymax=max(coord_CG(:,2));
xe=coord_CG(:,1);
ye=coord_CG(:,2);
nxx=200; nyy=200;
dx=(xmax-xmin)/nxx;
dy=(ymax-ymin)/nyy;
[xi,yi]=meshgrid(xmin:dx:xmax,ymin:dy:ymax);

%Create Mass Matrix  
Mmatrix = create_Mmatrix2d(jac,wnq,intma,iperiodic,npoin,nelem,ngl);

%Compute Courant Number
[~,vel,ds,dt] = compute_Courant(ue,ve,intma,coord,nelem,ngl,dt,Courant_max);
ntime=round(time_final/dt);
dt=time_final/ntime;
Courant=vel*dt/ds;
%pause(3);

%Initialize State Vector
q1=qe;
q0=qe;
qp=qe;
iplot=round(ntime/nplots);
iframe=0;

%Compute RK Time-Integration Coefficients
[a0,a1,beta] = compute_ti_coefficients(kstages);
rhs_rk4=zeros(npoin,4);

%Time Integration
for itime=1:ntime
    itime;
    time=time + dt;
    timec=time/(c);
    timec;

    if (mod(itime,iplot) == 0 )
        disp(['itime time courant = ',num2str(itime),' ',num2str(timec),' ',num2str(Courant)]);
    end
    
    for ik=1:kstages
       
        %------------Students Add Your Routines Here---------------%
        %------------Students Add Your Routines Here---------------%
        %Construct RHS vector
        rhs = construct_RHS_vector(qp,ue,ve,ksi_x,ksi_y,eta_x,eta_y,jac,...
        wnq,dpsi,intma,iperiodic,Mmatrix,face,nx,ny,jac_face,...
        mapL,mapR,npoin,nelem,nface,ngl,space_method,flux_method);
        %------------Students Add Your Routines Here---------------%
        %------------Students Add Your Routines Here---------------%

        %Evolve forward in Time
        qp=a0(ik)*q0 + a1(ik)*q1 + dt*beta(ik)*rhs;
        
        %Filter Solution
        if (mod(itime,ifilter) == 0)
          qp = apply_filter2D(qp,f,intma,iperiodic,jac,wnq,Mmatrix,npoin,nelem,ngl);
        end
        
        %Update
        rhs_rk4(:,ik)=rhs(:);
        q1=qp;
    end %ik

    %Do RK4 Solution
    if (kstages == 4)
        qp(:)=q0(:) + dt/6.0*( rhs_rk4(:,1) + 2*rhs_rk4(:,2)...
                    + 2*rhs_rk4(:,3) + rhs_rk4(:,4) );
    end
   
    %Update Q
    q0=qp;

    %PLOT Solution for MOVIES
    if (mod(itime,iplot) == 0 || itime == ntime)
        iframe=iframe + 1;
        %Compute gridpoint solution
        q_sol=zeros(npoin_CG,1);
        lhowm=zeros(npoin_CG,1);
        for i=1:npoin
            ip_CG=iperiodic_CG(DG_to_CG(i));
            q_sol(ip_CG)=q_sol(ip_CG) + q0(i);
            lhowm(ip_CG)=lhowm(ip_CG)+1;
        end %i
        for i=1:npoin_CG
           q_sol(i)=q_sol(i)/lhowm(i);
        end
        qi_movie(:,:,iframe)=griddata(xe,ye,q_sol,xi,yi,'cubic');
        time_movie(iframe)=timec;
    end %iplot

end %itime

%Compute Exact Solution
[qe,ue,ve] = exact_solution(coord,npoin,time,icase);

%Compute Norm
l2_norm=norm(q0-qe)/norm(qe);

%Compute gridpoint solution
q_sol=zeros(npoin_CG,1);
qe_sol=zeros(npoin_CG,1);
lhowm=zeros(npoin_CG,1);
for i=1:npoin
    ip_CG=iperiodic_CG(DG_to_CG(i));
    q_sol(ip_CG)=q_sol(ip_CG) + q0(i);
    qe_sol(ip_CG)=qe_sol(ip_CG) + qe(i);
    lhowm(ip_CG)=lhowm(ip_CG)+1;
end %i
for i=1:npoin_CG
    q_sol(i)=q_sol(i)/lhowm(i);
    qe_sol(i)=q_sol(i)/lhowm(i);
end

%Plot Solution
if (plot_solution == 1)
    h=figure;
    figure(h);
    qi=griddata(xe,ye,q_sol,xi,yi,'cubic');
    %[cl,h]=contourf(xi,yi,qi);
    mesh(xi,yi,qi);
    colorbar('SouthOutside');
    xlabel('X','FontSize',18);
    ylabel('Y','FontSize',18);
    axis image
    title_text=[space_method ', Ne = ' num2str(nelem) ', N = ' num2str(nop) ', Q = ' num2str(noq) ', L2 Norm = ' num2str(l2_norm)];
    title([title_text],'FontSize',18);
    set(gca, 'FontSize', 18);
    %file_ps=['se_n' num2str(nelem) 'p' num2str(nop)];
    %eval(['print ' file_ps ' -depsc']);
end

if (plot_movie == 1)
    figure;
    for i=1:iframe
        mesh(xi,yi,qi_movie(:,:,i));
        colorbar('SouthOutside');
        axis([-1 +1 -1 +1 0 1]);
        title_text=[space_method ', Ne = ' num2str(nelem) ', N = ' num2str(nop) ', Q = ' num2str(noq) ', Time = ' num2str(time_movie(i))];
        title([title_text],'FontSize',18);
        set(gca, 'FontSize', 18);
        M_i=getframe(gcf);
        M(i)=M_i;
        pause(0.2);
    end
end
if (store_movie == 1)
    movie2avi(M,'dg_inexact_TensorProduct_movie.avi','fps',5);
end

disp(['nop = ',num2str(nop),'  nelem = ',num2str(nelem) ]);
dt=dt/(c);
disp(['dt = ',num2str(dt),'  Courant = ',num2str(Courant) ]);
q_max=max(max(max(q0)));
q_min=min(min(min(q0)));
disp(['L2_Norm = ',num2str(l2_norm),' q_max = ',num2str(q_max),'  q_min = ',num2str(q_min) ]);
disp(['npoin = ',num2str(npoin),' npoin_CG = ',num2str(npoin_CG),'  npoin_DG = ',num2str(npoin_DG) ]);

toc