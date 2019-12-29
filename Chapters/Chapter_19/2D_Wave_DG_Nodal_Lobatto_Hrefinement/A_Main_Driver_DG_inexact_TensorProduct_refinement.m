%---------------------------------------------------------------------%
%This code computes the 2D Advection Equation using the CG/DG methods
%with 2nd/3rd Order RK and tensor product of 1D basis function with 
%Inexact Integration (Inexact Integration and Tensor-Product basis
%functions) with non-conforming h-refinement AMR
%Solver written by F.X. Giraldo on 6/2008
%AMR algorithm written by M.A. Kopera
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
clear all; 
close all;

tic

%Input Data
nel=2; %Number of Elements

nop=4;    %Interpolation Order
noq=nop; %DO NOT CHANGE!
space_method='dg'; %only works for DG
kstages=3;  %2=RK2, 3=RK3
dt=1; %time-step, Changes automatically to keep Courant_max fixed!
Courant_max=0.5;
time_final=0.25; %final time in revolutions
nplots=50; %Number of Frames in movie
store_movie=0;

icase=1; %case number: 1 is a Gaussian in CW, 2 is Gaussian along X, 
         %3 is Gaussian along Y, 4 is Gaussian along Diagonal, 
         %5 is square wave in CW, 6 is square wave along X.
xmu=0.05; %filtering strength: 1 is full strength and 0 is no filter
ifilter=1; %time-step frequency that the filter is applied.
iplot_grid = 1; %plot grid (1) or not (0)
iplot_result = 1;
inorm = 1;

maxlev = 4;
iadapt = 1;
icoarse = 1;
ref_threshold = [0.001, 0.001, 0.001, 0.001, 0.001];
cor_threshold = [0.001, 0.001, 0.001, 0.001, 0.001];
ref_freq = 5;

% allocation_excess = 21;

iref = [];
icor = [];

ae = 0;
for k=0:maxlev
    ae=ae+2^(2*k);
end
allocation_excess = ae

iref = [1:nel*nel*(allocation_excess-2^(2*maxlev))];
% iref = [1:25];

% iref = [1:100];
%  iref = [7,8,12,13,14,18,19]
%   iref = [7:9,12:14,17:19];
% iref = [7:9,12:14,17:19,42:45]

%icor = [26,42,106,125];
%icor = [75,79,95,99];
%icor = [26:125];
%iref = [1,2,8,57,64,63];
%iref = [28:29,36:37];
% iref = [10:15,18:23,26:31,34:39,42:47,50:55];
%iref = [10,19,28,37,46,55];
%iref=[12,14,16,18,23,25,27,29,32,34,36,38,43,45,47,49,52,54,56,58,63,65,67,69,72,74,76,78,83,85,87,89];
%iref=[42,43,52,53];
%iref = [12:19,22:29,32:39,42:49,52:59,62:69,72:79,82:89];
%iref = [42,52];
%iref = [23:27,33:37,43:47,53:57,63:67,73:77,115,116,121,122];
%
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
ntime=time_final/dt;
dt=dt*c;
time_final=time_final*c;
nelx=nel;
nely=nel;
nelem=nelx*nely; %Number of Elements
ngl=nop + 1;
nq=noq + 1;
npoin=(nop*nelx + 1)*(nop*nely + 1);
nelem=nelx*nely;
nboun=2*nelx + 2*nely;
nside=2*nelem + nelx + nely;
ncor = size(icor,2);
%ntime=time_final/dt;

%Compute LGL Points
[xgl,wgl]=legendre_gauss_lobatto(ngl);

%Compute Legendre Cardinal functions and derivatives
[psi,dpsi,xnq,wnq] = lagrange_basis(ngl,nq,xgl);

%Compute Filter Matrix
f = filter_init(ngl,xgl,xmu);

%Create Grid
[coord,intma,bsido,iperiodic,bsidoa] = create_grid_2d_dg(npoin,nelem,nboun, ...
                                nelx,nely,ngl,xgl,allocation_excess);
 
%Compute Side/Edge Information
[iside,jeside,jesideh] = create_side(intma,bsido,npoin,nelem,nboun,nside,ngl,allocation_excess);
[psideh,imapl,imapr] = create_side_dg(iside,intma,nside,nelem,ngl,nboun,bsidoa,jeside);

               
%Impose Periodicity
psideh=create_periodicity(iside,psideh,coord,nside,nboun);

%Create Faces
[face,facepa,facep] = create_face(nside,iside,psideh,allocation_excess);
jeside=reorganize_jeside(jeside,nelem);
jesideh=reorganize_jeside(jesideh,nelem);
%---
% [P1g2d,P2g2d,P3g2d,P4g2d,P1s2d,P2s2d,P3s2d,P4s2d] = l2_projection_2D_init(ngl,ngl+1);
% [qa,ua,va]=exact_solution(coord,npoin*4,0,icase);
%---

%Refine Mesh
elmptr = zeros(1,nelem);
tc = zeros(4,allocation_excess*nelem); %tree children
tp = zeros(1,allocation_excess*nelem); %tree parents
tm = zeros(1,allocation_excess*nelem); %tree markers
tl = zeros(1,allocation_excess*nelem); %tree levels

disp('start refine');
tm(1:nelem) = 1;
refel = nelem+1;
refpt = npoin+1;
% [refel,refpt,coord,intma,face,facep,facepa,iface,jeside,jesideh,iperiodic,tc,tp,tm,tl] = refine_elements_CG_DG(iref,refel,refpt,...
%                                                                      coord,intma,ngl,xgl,jeside,jesideh,nside,face,facepa,facep,...
%                                                                      iperiodic,tc,tp,tm,tl);

[refel,refpt,coord,intma,face,facepa,iface,jeside,jesideh,iperiodic,tc,tp,tm,tl] = refine_elements_CG_DG_new(iref,refel,refpt,coord,...
                                                intma,ngl,xgl,jeside,jesideh,nside,face,facepa,iperiodic,tc,tp,tm,tl);
nface = iface;
nelemT=refel-1;
npoinT=refpt-1;

% disp('start coarsen');
% [tm,facepa,face,jeside] = coarsen_elements(tc,tp,tm,icor,ncor,facepa,jeside,face);
% [tm,facepa,face,jeside] = refine_elements(tc,tm,7,1,facepa,jeside,face);

disp('start arrange face');
face = arrange_face(face,nface);


disp('start sfc');
%Generate the Space Filling Curve - SFC
[sfc,nsfc] = space_filling_curve_tree(nelem,tc,tm);
[ffc,nffc] = face_filling_curve(nface,facepa);

disp('start normals');
% [nx,ny,jac_side] = compute_normals(face,intma,coord,...
%                    nface,ngl,nq,wnq,psi,dpsi,ffc,nffc);
[nx,ny,jac_side] = compute_normals_all(face,intma,coord,...
                   nface,ngl,nq,wnq,psi,dpsi);
%Plot Grid
% plot_grid(coord,intma,sfc,nsfc,ngl,nop,iplot_grid);
% 
%***********

%Compute Metric Terms
% [ksi_x,ksi_y,eta_x,eta_y,jac] = metrics(coord,intma,psi,dpsi,wnq,nelemT,ngl,nq,sfc,nsfc);
[ksi_x,ksi_y,eta_x,eta_y,jac] = metrics_all(coord,intma,psi,dpsi,wnq,nelemT,ngl,nq);

%Compute Exact Solution
time=0;

disp('start exact');
[qa,ua,va]=exact_solution(coord,npoinT,time,icase);

q0=zeros(nelem,ngl,ngl);
qe=zeros(nelem,ngl,ngl);
ue=zeros(nelem,ngl,ngl);
ve=zeros(nelem,ngl,ngl);
qerr=zeros(nelem,ngl,ngl);

% Reshape Arrays
for e=1:nelemT
    for j=1:ngl
    for i=1:ngl
        ip=intma(e,i,j);
        qe(e,i,j)=qa(ip);
        ue(e,i,j)=ua(ip);
        ve(e,i,j)=va(ip);
    end %i
    end %j
end %e

pflag = zeros(npoinT,1);

% for is=1:nelem
%     ie = sfc(is); 
%     for i=1:ngl
%         for j=1:ngl
%             pflag(intma(ie,i,j)) = 1;
%         end
%     end
%     for f=1:4
%        ifa = jeside(ie,f);
%        if(face(ifa,9)==1)
%            
%            
%        end
%     end
% end
% 
%Plot Exact Solution
xmin=min(coord(:,1)); xmax=max(coord(:,1));
ymin=min(coord(:,2)); ymax=max(coord(:,2));

nxx=200; nyy=200;
dx=(xmax-xmin)/nxx;
dy=(ymax-ymin)/nyy;
[xi,yi]=meshgrid(xmin:dx:xmax,ymin:dy:ymax);

if iplot_result==1
h5 = figure;
figure(h5);
hold on;
for is=1:nsfc
    ie=sfc(is);
    for i=1:ngl
       x1(i) = coord(intma(ie,i,1),1);
       for j=1:ngl
           qq(j,i) = qe(ie,i,j);
            y1(j) = coord(intma(ie,1,j),2); 
       end
    end
    h5 = surf(x1,y1,qq);
    view(3);
end
title_text=['Exact Solution For: Ne = ' num2str(nelemT) ', N = ' num2str(nop) ', Q = ' num2str(noq)];
title([title_text],'FontSize',18);
xlabel('X','FontSize',18);
ylabel('Y','FontSize',18);
end

%Create Mass Matrix
[Mmatrix, Mmatrix_global] = create_Mmatrix2d_dg_TensorProduct_inexact_all(jac,intma,iperiodic,npoinT,nelemT,ngl,face,nface);  
% for ie=1:nelem
%     e=sfc(ie);
%     Mmatrix_inv(e,:,:)=1/Mmatrix(e,:,:);
% end
for e=1:nelemT
    Mmatrix_inv(e,:,:)=1/Mmatrix(e,:,:);
end

%Compute projection matrices
[P1g,P2g,P1s,P2s] = l2_projection_init(ngl,ngl+1);
[P1g2d,P2g2d,P3g2d,P4g2d,P1s2d,P2s2d,P3s2d,P4s2d] = l2_projection_2D_init(ngl,ngl+1);

%Compute Courant Number
[Courant,vel,ds,dt] = compute_Courant(ue,ve,intma,coord,ngl,dt,Courant_max,sfc,nsfc);
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

mass0 = total_mass(qp,jac,sfc,nsfc,ngl);
ii=1;
% hm = figure;
if iplot_result==1
scrsz = get(0,'ScreenSize');
if iplot_grid==1
    movie1 = figure('Position',[1,scrsz(4)/2,1280,480]);
else
    movie1 = figure('Position',[1,scrsz(4)/2,640,480]); 
end
end

% if iplot_result==1
%    %PLOT Solution for MOVIES
%        iframe= 1;
%        figure(movie1);
% %        qi_movie(:,:,:,iframe) = q0(:,:,:);
%         time_movie(iframe)=0;
%         if iplot_grid==1
%             subplot(1,2,1);
%         end
%         axis([-1 +1 -1 +1 0 1]);
%         title_text=[space_method ' : Ne = ' num2str(nelemT) ', N = ' num2str(nop) ', Q = ' num2str(noq) ', Time = ' num2str(time_movie(iframe))];
%         title([title_text],'FontSize',18);
%         set(gca, 'FontSize', 18);
%         %hold on;
%         for is=1:nsfc
%             ie=sfc(is);
%             for i=1:ngl
%                x1(i) = coord(intma(ie,i,1),1);
%                for j=1:ngl
% %                    qq(j,i) = qi_movie(ie,i,j,iframe);
%                     qq(j,i) = q0(ie,i,j);
%                     y1(j) = coord(intma(ie,1,j),2); 
%                end
%             end
%             surf(x1,y1,qq);
%             hold on;
%             view([1 -1 1]);
%             axis([-1 +1 -1 +1 0 1]);
%             
%         end
%         hold off;
%         if iplot_grid==1
%             hh2 = subplot(1,2,2);
%             cla(hh2);
%             x=zeros(5,1);
%             y=zeros(5,1);
%             hold on;
%             for is=1:nsfc
%                 el=sfc(is);
%                 for j=1:ngl-1
%                     for i=1:ngl-1
%                         i1=intma(el,i,j);
%                         i2=intma(el,i+1,j);
%                         i3=intma(el,i+1,j+1);
%                         i4=intma(el,i,j+1);
%                         x(1)=coord(i1,1); y(1)=coord(i1,2);
%                         x(2)=coord(i2,1); y(2)=coord(i2,2);
%                         x(3)=coord(i3,1); y(3)=coord(i3,2);
%                         x(4)=coord(i4,1); y(4)=coord(i4,2);
%                         x(5)=coord(i1,1); y(5)=coord(i1,2);
%                         plot_handle=plot(x,y,'-r');
%                         set(plot_handle,'LineWidth',1.5);
%                     end
%                 end
%                  i1=intma(el,1,1);
%                  i2=intma(el,ngl,1);
%                  i3=intma(el,ngl,ngl);
%                  i4=intma(el,1,ngl);
%                  x(1)=coord(i1,1); y(1)=coord(i1,2);
%                  x(2)=coord(i2,1); y(2)=coord(i2,2);
%                  x(3)=coord(i3,1); y(3)=coord(i3,2);
%                  x(4)=coord(i4,1); y(4)=coord(i4,2);
%                  x(5)=coord(i1,1); y(5)=coord(i1,2);
%                  plot_handle=plot(x,y,'-b');
%                  set(plot_handle,'LineWidth',2);
%             end
%             xlabel('X','FontSize',18);
%             ylabel('Y','FontSize',18);
%             set(gca,'FontSize',16);
%             axis image
%             hold off;
%         end
%         
%         M_i=getframe(gcf);
% %         MakeQTMovie('addfigure');
% %         MakeQTMovie('addaxes');
%         hold off;
%         M(iframe)=M_i;
%        
%    
%     end
% MakeQTMovie('start', 'test_movie.mov');
% MakeQTMovie('quality',1.0);
% MakeQTMovie('size',[640 480]);
%Time Integration
itime=0;
% for itime=1:ntime
while time<time_final
   itime=itime+1;
   time=time + dt;
   timec=time/(c);
   timec;
   
   
   %--- AMR ---
   
   [refel,refpt,coord,intma,face,facepa,iface,jeside,jesideh,...
          iperiodic,tc,tp,tm,tl,nx,ny,jac_side,ksi_x,ksi_y,...
          eta_x,eta_y,jac,Mmatrix,Mmatrix_inv,qp,q0,q1,ue,ve,ffc,nffc,sfc,...
          nsfc,Courant,vel,ds,dt,nelemT,npoinT,nface] = ...
                    AMR(ref_threshold,cor_threshold,iadapt,icoarse,itime,ref_freq,...
                        sfc,nsfc,ffc,nffc,q0,q1,qp,ue,ve,refel,refpt,xgl,...
                        jeside,jesideh,iface,face,facepa,iperiodic,tc,tp,tm,tl,...
                        nq,wnq,psi,dpsi,nx,ny,jac_side,intma,coord,Courant_max,...
                        maxlev,nelem,nface,ngl,ksi_x,ksi_y,eta_x,eta_y,jac,...
                        Mmatrix,Mmatrix_inv,nelemT,npoinT,Courant,vel,ds,dt,...
                        P1s2d,P2s2d,P3s2d,P4s2d,P1g2d,P2g2d,P3g2d,P4g2d);
   %--- end AMR ---
   

    for ik=1:kstages
        switch kstages
            case 1  %euler
                a0=1;
                a1=0;
                beta=1;
            case 2  %RK2
                switch ik
                    case 1
                        a0=1;
                        a1=0;
                        beta=1;
                    case (2)
                        a0=0.5;
                        a1=0.5;
                        beta=0.5;
                    end %ik
            case 3 %RK3
                switch ik
                    case 1
                        a0=1;
                        a1=0;
                        beta=1;
                    case (2)
                        a0=3.0/4.0;
                        a1=1.0/4.0;
                        beta=1.0/4.0;
                    case (3)
                        a0=1.0/3.0;
                        a1=2.0/3.0;
                        beta=2.0/3.0;
                    end %ik
      end %kstages
      dtt=dt*beta;
      
      %Construct RHS
      rhs = compute_rhs_dg_TensorProduct_inexact(qp,ue,ve,ksi_x,ksi_y,eta_x,eta_y,jac,...
	        dpsi,nelemT,ngl,sfc,nsfc);
        

        
      %Construct Communicator: DSS for CG or Fluxes for DG
      if (space_method == 'cg')
          rhs = apply_dss(rhs,intma,iperiodic,Mmatrix_global,ngl,npoinT,sfc,nsfc,face,ffc,nffc);
      elseif (space_method == 'dg')
        rhs = compute_flux_dg_TensorProduct_inexact_1(rhs,qp,ue,ve,nx,ny,jac_side,...
            ngl,imapl,imapr,ffc,nffc,face,P1g,P2g,P1s,P2s);
%         rhs = compute_flux_dg_TensorProduct_inexact_2(rhs,qp,ue,ve,psideh,nx,ny,jac_side,psi,...
%                nside,ngl,nq,imapl,imapr,ffc,nffc,face,wgl);
         for is=1:nsfc
             ie=sfc(is);
             for j=1:ngl
                 for i=1:ngl
                    rhs(ie,i,j)=rhs(ie,i,j)/Mmatrix(ie,i,j);
                 end
             end
         end
        
      end %if
      
      %Evolve forward in Time
      for is=1:nsfc
            e=sfc(is);
          qp(e,:,:)=a0*q0(e,:,:) + a1*q1(e,:,:) + dtt*rhs(e,:,:);
      end

      %Filter Solution
       if (mod(itime,ifilter) == 0)
          rhs = apply_filter2D_dg(qp,f,nelemT,sfc,nsfc,ngl);
          if (space_method == 'cg')
               rhs=rhs.*jac;
               rhs = apply_dss(rhs,intma,iperiodic,Mmatrix_global,ngl,npoinT,sfc,nsfc,face,ffc,nffc);
               %rhs = apply_dss(rhs,intma,iperiodic,Mmatrix_global,ngl,npoin,nelem);
          end
          qp=rhs;
       end

      %Update
      q1=qp;
      
   end %ik

   
   %Update Q
   q0=qp;
   
   if inorm==1
   %Compute Exact Solution
    [qa,ua,va] = exact_solution(coord,npoinT,time,icase);
    
    % Reshape Arrays
    for is=1:nsfc
        e=sfc(is);
        for j=1:ngl
        for i=1:ngl
            ip=intma(e,i,j);
            qe(e,i,j)=qa(ip);
            ue(e,i,j)=ua(ip);
            ve(e,i,j)=va(ip);
            qerr(e,i,j) = abs(qe(e,i,j)-q0(e,i,j));
        end %i
        end %j
    end %e

    %Compute Norms
    top=0;
    bot=0;
    for is=1:nsfc
        e=sfc(is);
        for j=1:ngl
        for i=1:ngl
            top=top + (q0(e,i,j)-qe(e,i,j))^2;
            bot=bot + qe(e,i,j)^2;
        end %i
        end %j
    end %e
    l2_norm=sqrt( top/bot );
    l_max = max(max(max(qerr)));

   
    mass(itime) = total_mass(qp,jac,sfc,nsfc,ngl);
     disp(['itime time courant mass l2 lmax= ',num2str(itime),' ',num2str(timec),' ',num2str(Courant),' ',num2str(mass(itime)-mass0),' ',num2str(l2_norm),' ',num2str(l_max)]);
 
   end
    if iplot_result==1
   %PLOT Solution for MOVIES
   if (mod(itime,iplot) == 0)
       iframe=iframe + 1;
       figure(movie1);
%        qi_movie(:,:,:,iframe) = q0(:,:,:);
        time_movie(iframe)=timec;
        if iplot_grid==1
            subplot(1,2,1);
        end
        axis([-1 +1 -1 +1 0 1]);
        title_text=[space_method ' : Ne = ' num2str(nelemT) ', N = ' num2str(nop) ', Q = ' num2str(noq) ', Time = ' num2str(time_movie(iframe))];
        title([title_text],'FontSize',18);
        set(gca, 'FontSize', 18);
        %hold on;
        for is=1:nsfc
            ie=sfc(is);
            for i=1:ngl
               x1(i) = coord(intma(ie,i,1),1);
               for j=1:ngl
%                    qq(j,i) = qi_movie(ie,i,j,iframe);
                    qq(j,i) = q0(ie,i,j);
                    y1(j) = coord(intma(ie,1,j),2); 
               end
            end
            surf(x1,y1,qq);
            hold on;
            view([1 -1 1]);
            axis([-1 +1 -1 +1 0 1]);
            xlabel('X','FontSize',18);
            ylabel('Y','FontSize',18);
        end
        hold off;
        if iplot_grid==1
            hh2 = subplot(1,2,2);
            cla(hh2);
            x=zeros(5,1);
            y=zeros(5,1);
            hold on;
            for is=1:nsfc
                el=sfc(is);
                for j=1:ngl-1
                    for i=1:ngl-1
                        i1=intma(el,i,j);
                        i2=intma(el,i+1,j);
                        i3=intma(el,i+1,j+1);
                        i4=intma(el,i,j+1);
                        x(1)=coord(i1,1); y(1)=coord(i1,2);
                        x(2)=coord(i2,1); y(2)=coord(i2,2);
                        x(3)=coord(i3,1); y(3)=coord(i3,2);
                        x(4)=coord(i4,1); y(4)=coord(i4,2);
                        x(5)=coord(i1,1); y(5)=coord(i1,2);
                        plot_handle=plot(x,y,':r');
                        set(plot_handle,'LineWidth',1.5);
                    end
                end
                 i1=intma(el,1,1);
                 i2=intma(el,ngl,1);
                 i3=intma(el,ngl,ngl);
                 i4=intma(el,1,ngl);
                 x(1)=coord(i1,1); y(1)=coord(i1,2);
                 x(2)=coord(i2,1); y(2)=coord(i2,2);
                 x(3)=coord(i3,1); y(3)=coord(i3,2);
                 x(4)=coord(i4,1); y(4)=coord(i4,2);
                 x(5)=coord(i1,1); y(5)=coord(i1,2);
                 plot_handle=plot(x,y,'-b');
                 set(plot_handle,'LineWidth',2);
            end
            xlabel('X','FontSize',18);
            ylabel('Y','FontSize',18);
            set(gca,'FontSize',16);
            axis image
            hold off;
        end
        
        M_i=getframe(gcf);
%         MakeQTMovie('addfigure');
%         MakeQTMovie('addaxes');
        hold off;
        M(iframe)=M_i;
       
   end 
    end
%here was AMR

   
end %itime

if inorm==1
plot_grid_mass(coord,intma,sfc,nsfc,ngl,nop,iplot_grid,time_final,itime,mass,mass0,nelem);
end
% dt = time_final/itime;
% t = [1:itime]*dt;
% figure;
% plot(t,mass-mass0);

%Compute Exact Solution
[qa,ua,va] = exact_solution(coord,npoinT,time,icase);

% Reshape Arrays
for is=1:nsfc
    e=sfc(is);
    for j=1:ngl
    for i=1:ngl
        ip=intma(e,i,j);
        qe(e,i,j)=qa(ip);
        ue(e,i,j)=ua(ip);
        ve(e,i,j)=va(ip);
    end %i
    end %j
end %e

%Compute Norm
top=0;
bot=0;
for is=1:nsfc
    e=sfc(is);
    for j=1:ngl
    for i=1:ngl
        top=top + (q0(e,i,j)-qe(e,i,j))^2;
        bot=bot + qe(e,i,j)^2;
    end %i
    end %j
end %e
l2_norm=sqrt( top/bot );

if iplot_result==1
h6 = figure;
figure(h6);
hold on;
for is=1:nsfc
    ie=sfc(is);
    for i=1:ngl
       x1(i) = coord(intma(ie,i,1),1);
       for j=1:ngl
           qq(j,i) = q0(ie,i,j);
            y1(j) = coord(intma(ie,1,j),2); 
       end
    end
    h5 = surf(x1,y1,qq);
    view(3);
end
title_text=['Solution For: Ne = ' num2str(nelemT) ', N = ' num2str(nop) ', Q = ' num2str(noq) ', time' num2str(time)];
title([title_text],'FontSize',18);
xlabel('X','FontSize',18);
ylabel('Y','FontSize',18);
end

if (store_movie == 1)
    movie2avi(M,'dg_inexact_TensorProduct_movie.avi','fps',10);
end

% MakeQTMovie('finish');

nop
nelem
dt=dt/(c);
dt
Courant
q_max=max(max(max(q0)))
q_min=min(min(min(q0)))
l2_norm

toc
