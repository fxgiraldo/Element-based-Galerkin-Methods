%---------------------------------------------------------------------%
%This code computes the 2D Advection Equation using the CG/DG methods
%with 2nd Order RK and tensor product of 1D basis function with 
%Inexact Integration (Inexact Integration and Tensor-Product basis
%functions)
%Written by F.X. Giraldo on 6/2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
clear all; 
close all;

tic

%Input Data
nel=5; %Number of Elements
nop=6;    %Interpolation Order
noq=nop; %DO NOT CHANGE!
space_method='dg'; %CG or DG
kstages=3;  %2=RK2, 3=RK3
dt=1; %time-step, Changes automatically to keep Courant_max fixed!
Courant_max=0.5;
time_final=1; %final time in revolutions
nplots=50; %Number of Frames in movie
store_movie=1;

icase=4; %case number: 1 is a Gaussian in CW, 2 is Gaussian along X, 
         %3 is Gaussian along Y, 4 is Gaussian along Diagonal, 
         %5 is square wave in CW, 6 is square wave along X.
xmu=0.0; %filtering strength: 1 is full strength and 0 is no filter
ifilter=1; %time-step frequency that the filter is applied.
iplot_grid = 1; %plot grid (1) or not (0)
maxlev = 1;
iadapt = 0;

iref = [];


%iref = [1:25];
  iref = [7:9,12:14,17:19];
% iref = [7:9,12:14,17:19,42:45];
icor = [];
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
                                nelx,nely,ngl,xgl);
 
%Compute Side/Edge Information
[iside,jeside,jesideh] = create_side(intma,bsido,npoin,nelem,nboun,nside,ngl);
[psideh,imapl,imapr] = create_side_dg(iside,intma,nside,nelem,ngl,nboun,bsidoa,jeside);

               
%Impose Periodicity
psideh=create_periodicity(iside,psideh,coord,nside,nboun);

%Create Faces
[face,facepa,facep] = create_face(nside,iside,psideh);

%---
% [P1g2d,P2g2d,P3g2d,P4g2d,P1s2d,P2s2d,P3s2d,P4s2d] = l2_projection_2D_init(ngl,ngl+1);
% [qa,ua,va]=exact_solution(coord,npoin*4,0,icase);
%---

%Refine Mesh
elmptr = zeros(1,nelem);
tc = zeros(4,4*nelem); %tree children
tp = zeros(1,4*nelem); %tree parents
tm = zeros(1,4*nelem); %tree markers
tl = zeros(1,4*nelem); %tree levels

tm(1:nelem) = 1;
refel = nelem+1;
refpt = npoin+1;
[refel,refpt,coord,intma,face,facep,facepa,iface,jeside,jesideh,iperiodic,tc,tp,tm,tl] = refine_elements_CG_DG(iref,refel,refpt,...
                                                                     coord,intma,ngl,xgl,jeside,jesideh,nside,face,facepa,facep,...
                                                                     iperiodic,tc,tp,tm,tl);
nface = iface;
nelemT=refel-1;
npoinT=refpt-1;

[tm,facepa,face,jeside] = coarsen_elements(tc,tp,tm,icor,ncor,facepa,jeside,face);
% [tm,facepa,face,jeside] = refine_elements(tc,tm,7,1,facepa,jeside,face);

face = arrange_face(face,nface);

%Generate the Space Filling Curve - SFC
[sfc,nsfc] = space_filling_curve_tree(nelem,tc,tm);
[ffc,nffc] = face_filling_curve(nface,facepa);

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
[Courant,vel,ds,dt] = compute_Courant(ua,va,intma,coord,ngl,dt,Courant_max,sfc,nsfc);
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

hm = figure;

% MakeQTMovie('start', 'test_movie.mov');
% MakeQTMovie('quality',1.0);
% MakeQTMovie('size',[640 480]);
%Time Integration
for itime=1:ntime
   itime;
   time=time + dt;
   timec=time/(c);
   timec;
   

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
  
   %PLOT Solution for MOVIES
   if (mod(itime,iplot) == 0)
       iframe=iframe + 1;
       
       qi_movie(:,:,:,iframe) = q0(:,:,:);
       time_movie(iframe)=timec;
       
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
                   qq(j,i) = qi_movie(ie,i,j,iframe);
                    y1(j) = coord(intma(ie,1,j),2); 
               end
            end
            surf(x1,y1,qq);
            hold on;
            view([1 -1 1]);
            axis([-1 +1 -1 +1 0 1]);
        end
        M_i=getframe(gcf);
%         MakeQTMovie('addfigure');
%         MakeQTMovie('addaxes');
        hold off;
        M(iframe)=M_i;
       
   end 

   
%    %mark for mesh adaptation
%    if(mod(itime,10) == 0 && iadapt == 1)
%        iref = 0;
%        ir = 0;
%        icor = 0;
%        ic = 0;      
%        for is=1:nelem
%            ie = sfc(is);
%            qoi = max(max(qp(ie,:,:)));
% 
%            if (qoi > 0.5 && tl(ie)<maxlev)
%                ir=ir+1;
%                iref(ir) = ie;
%            elseif (qoi < 0.1 && tl(ie)>0)
%                ic=ic+1;
%                icor(ic) = ie;
%            end
%        end
%        iref
%        icor
%        nref = ir;
%        ncor = ic;
%        if(nref>0)
%             [tm,facepa,face,jeside] = refine_elements(tc,tm,iref,nref,facepa,jeside,face);
%             [qp,ue,ve] = refine_project_data(qp,ue,ve,iref,nref,P1s2d,P2s2d,P3s2d,P4s2d,ngl,tc);
%        elseif(ncor>0)
%             [tm,facepa,face,jeside] = coarsen_elements(tc,tp,tm,icor,ncor,facepa,jeside,face);
%             [qp,ue,ve] = coarse_project_data(qp,ue,ve,icor,ncor,P1g2d,P2g2d,P3g2d,P4g2d,ngl,tc,tp);
%        end
%        q0=qp;
%        %Generate the Space Filling Curve - SFC
%        [sfc,nsfc] = space_filling_curve_tree(nelem,tc,tm);
%        [ffc,nffc] = face_filling_curve(nface,facepa);
% % %        Normals
% %        [nx,ny,jac_side] = compute_normals(face,intma,coord,...
% %                    nface,ngl,nq,wnq,psi,dpsi,ffc,nffc);
% % %        Compute Metric Terms
% %        [ksi_x,ksi_y,eta_x,eta_y,jac] = metrics(coord,intma,psi,dpsi,wnq,nelemT,ngl,nq,sfc,nsfc);
% % %        Mass Matrix
% %        [Mmatrix, Mmatrix_global] = create_Mmatrix2d_dg_TensorProduct_inexact(jac,intma,iperiodic,npoinT,nelemT,ngl,sfc,nsfc,face,ffc,nffc);  
% %         for ie=1:nelem
% %             e=sfc(ie);
% %             Mmatrix_inv(e,:,:)=1/Mmatrix(e,:,:);
% %         end
% 
%    end

   %mark for mesh adaptation
%    if(mod(itime,10) == 0 && iadapt == 1)
%        iref = 0;
%        ir = 0;
%        icor = 0;
%        ic = 0;      
%        for is=1:nelem
%            ie = sfc(is);
%            qoi = max(max(qp(ie,:,:)));
% 
%            if (qoi > 0.9 && tl(ie)<maxlev)
%                iref = ie;
%                ir=1;
%            elseif (qoi < 0.1 && tl(ie)>0)
%                icor = ie;
%                ic=1;
%            end
%        end
%        iref
%        icor
%        nref = ir;
%        ncor = ic;
%        if(nref>0)
%             [tm,facepa,face,jeside] = refine_elements(tc,tm,iref,nref,facepa,jeside,face);
%             [qp,ue,ve] = refine_project_data(qp,ue,ve,iref,nref,P1s2d,P2s2d,P3s2d,P4s2d,ngl,tc);
%        elseif(ncor>0)
%             [tm,facepa,face,jeside] = coarsen_elements(tc,tp,tm,icor,ncor,facepa,jeside,face);
%             [qp,ue,ve] = coarse_project_data(qp,ue,ve,icor,ncor,P1g2d,P2g2d,P3g2d,P4g2d,ngl,tc,tp);
%        end
%        q0=qp;
%        %Generate the Space Filling Curve - SFC
%        [sfc,nsfc] = space_filling_curve_tree(nelem,tc,tm);
%        [ffc,nffc] = face_filling_curve(nface,facepa);
% 
% 
%    end

% %  if(itime==10)
% % %       [tm,facepa,face,jeside] = coarsen_elements(tc,tp,tm,icor,ncor,facepa,jeside,face);
% % %       [qp,ue,ve] = coarse_project_data(qp,ue,ve,icor,ncor,P1g2d,P2g2d,P3g2d,P4g2d,ngl,tc,tp);
% % %       q0=qp;
% % %       %Generate the Space Filling Curve - SFC
% % %       [sfc,nsfc] = space_filling_curve_tree(nelem,tc,tm);
% % %       [ffc,nffc] = face_filling_curve(nface,facepa);
% % % 
% % %    elseif(itime==200)
% %       [tm,facepa,face,jeside] = refine_elements(tc,tm,7,1,facepa,jeside,face);
% %       [qp,ue,ve] = refine_project_data(qp,ue,ve,7,1,P1s2d,P2s2d,P3s2d,P4s2d,ngl,tc);
% %       q0=qp;
% %       face = arrange_face(face,nface);
% %       %Generate the Space Filling Curve - SFC
% %       [sfc,nsfc] = space_filling_curve_tree(nelem,tc,tm);
% %       [ffc,nffc] = face_filling_curve(nface,facepa);
% % 
% %    end
  
%    if(itime==100)
%       [tm,facepa,face,jeside] = coarsen_elements(tc,tp,tm,icor,ncor,facepa,jeside,face);
%       [qp,ue,ve] = coarse_project_data(qp,ue,ve,icor,ncor,P1g2d,P2g2d,P3g2d,P4g2d,ngl,tc,tp);
%       q0=qp;
%       %Generate the Space Filling Curve - SFC
%       [sfc,nsfc] = space_filling_curve_tree(nelem,tc,tm);
%       [ffc,nffc] = face_filling_curve(nface,facepa);
% 
%    elseif(itime==200)
%       [tm,facepa,face,jeside] = refine_elements(tc,tm,7,1,facepa,jeside,face);
%       [qp,ue,ve] = refine_project_data(qp,ue,ve,7,1,P1s2d,P2s2d,P3s2d,P4s2d,ngl,tc);
%       q0=qp;
%       %Generate the Space Filling Curve - SFC
%       [sfc,nsfc] = space_filling_curve_tree(nelem,tc,tm);
%       [ffc,nffc] = face_filling_curve(nface,facepa);
% 
%    end
   
end %itime

plot_grid_mass(coord,intma,sfc,nsfc,ngl,nop,iplot_grid,time_final,itime,mass,mass0,nelem);

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

%Compute gridpoint solution
% q_sol=zeros(npoinT,1);
% qe_sol=zeros(npoinT,1);
% lhowm=zeros(npoinT,1);
% for is=1:nsfc
%     ie=sfc(is);
%     
%     for j=1:ngl
%     for i=1:ngl
%       ip=intma(ie,i,j);
%       lhowm(ip)=lhowm(ip)+1;
%       q_sol(ip)=q_sol(ip) + q0(ie,i,j);
%       qe_sol(ip)=qe_sol(ip) + qe(ie,i,j);
%     end %i
%     end %j
% end
% for i=1:npoinT
%    q_sol(i)=q_sol(i)/lhowm(i);
%    qe_sol(i)=qe_sol(i)/lhowm(i);
% end

% %Plot Solution
% h=figure;
% figure(h);
% qi=griddata(xe,ye,q_sol,xi,yi,'cubic');
% [cl,h]=contourf(xi,yi,qi);
% %colorbar('SouthOutside');
% xlabel('X','FontSize',18);
% ylabel('Y','FontSize',18);
% axis image
% title_text=[space_method ' : Ne = ' num2str(nelem) ', N = ' num2str(nop) ', Q = ' num2str(noq) ', L2 Norm = ' num2str(l2_norm)];
% title([title_text],'FontSize',18);
% set(gca, 'FontSize', 18);
% %file_ps=['se_n' num2str(nelem) 'p' num2str(nop)];
% %eval(['print ' file_ps ' -depsc']);

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
title_text=['Solution For: Ne = ' num2str(nelemT) ', N = ' num2str(nop) ', Q = ' num2str(noq) ', time' time];
title([title_text],'FontSize',18);


% figure;
% for k=1:iframe
%     mesh(xi,yi,qi_movie(:,:,i));
%     axis([-1 +1 -1 +1 0 1]);
%     title_text=[space_method ' : Ne = ' num2str(nelemT) ', N = ' num2str(nop) ', Q = ' num2str(noq) ', Time = ' num2str(time_movie(k))];
%     title([title_text],'FontSize',18);
%     set(gca, 'FontSize', 18);
%     %hold on;
%     for is=1:nsfc
%         ie=sfc(is);
%         for i=1:ngl
%            x1(i) = coord(intma(ie,i,1),1);
%            for j=1:ngl
%                qq(j,i) = qi_movie(ie,i,j,k);
%                 y1(j) = coord(intma(ie,1,j),2); 
%            end
%         end
%         surf(x1,y1,qq);
%         hold on;
%         view([1 -1 1]);
%         axis([-1 +1 -1 +1 0 1]);
%     end
%     M_i=getframe(gcf);
%     hold off;
%     M(k)=M_i;
%     pause(0.2);
% end
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
