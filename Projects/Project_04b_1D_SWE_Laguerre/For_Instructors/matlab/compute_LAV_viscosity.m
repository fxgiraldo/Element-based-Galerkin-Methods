%---------------------------------------------------------------------%
%This function computes the LGL grid and elements.
%Written by F.X. Giraldo on 10/2003
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5502
%---------------------------------------------------------------------%
function [visc_elem,min_visc_elem,max_visc_elem] = compute_LAV_viscosity(rhs,qp,qb,intma,coord,npoin,nelem,ngl,gravity,visc,C1,C2)

%Initialize
visc_elem=zeros(npoin,2);
qh=zeros(npoin,1);
qu=zeros(npoin,1);

%Extract Variables
qh(:)=qp(:,1) + qb(:);
qu(:)=qp(:,2);

%Mean Values
h_mean=sum(qh)/npoin;
u_mean=sum(qu)/npoin;

% %Compute Global Infty-Norm
h_infty=norm(qh(:) - h_mean,inf);
u_infty=norm(qu(:) - u_mean,inf);

%Loop through Elements
for e=1:nelem
   
   %Store Coordinates
   for i=1:ngl(e)
      I=intma(i,e);
      h=qp(I,1);
      u=qp(I,2)/(h+qb(I));
      wave_speed(i)=abs(u) + sqrt(gravity*h);
      h_rhs(i)=rhs(I,1);
      u_rhs(i)=rhs(I,2);
   end
   
   %Jacobians
   dx_mean=(coord(ngl(e),e)-coord(1,e))/(ngl(e)-1);

   %Infinity-norm of Element
   h_norm_infty=norm(h_rhs,Inf);
   u_norm_infty=norm(u_rhs,Inf);
   max_wave_speed=max(wave_speed);
   
   %Compute Viscosity
   mu_max=C2*0.5*dx_mean*max_wave_speed;
   mu_1=C1*dx_mean^2*max(h_norm_infty/h_infty,u_norm_infty/u_infty);
   mu_n=max(0,min(mu_max,mu_1));
   mu_h=mu_n;
   mu_u=mu_n;
   
   %Store Viscosity
   for i=1:ngl(e)
       I=intma(i,e);
       visc_elem(I,1)=mu_h; 
       visc_elem(I,2)=mu_u; 
   end
end %e

%Constant Viscosity
if (C1<0 || C2<0)
    visc_elem(:,1)=visc;
    visc_elem(:,2)=visc;
end

min_visc_elem=min(visc_elem(:,1));
max_visc_elem=max(visc_elem(:,1));

%disp(['visc_min =  ',num2str(min_visc_elem),' visc_max = ', num2str(max_visc_elem)]);
   
