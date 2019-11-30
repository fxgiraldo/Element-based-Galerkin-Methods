%---------------------------------------------------------------------%
%This code computes the 1D Advection Equation using the 
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
nelem=3; %Number of Elements
nop_min=1;    %Interpolation Order
nop_max=8;
%Store Constants
d=3; %space-dimension
nelem_max=1000^(1/d);
%nelem_max=nelem;

for nelem=nelem_max:nelem_max
    for nop=nop_min:nop_max
        ngl=nop + 1;
        npoin_dg(nop)=(ngl*nelem)^d;
        npoin_cg(nop)=(nop*nelem + 1)^d;
        npoin_hdg(nop)=d*(nop*nelem+1)^(d-1)*(nelem+1);
        ndof_dg(nop)=(ngl*nelem)^(3*d);
        ndof_cg(nop)=(nop*nelem + 1)^(3*d);
        ndof_hdg(nop)=(d*(nop*nelem+1)^(d-1)*(nelem+1))^3 + nelem*ngl^(3*d);
        nop_vector(nop)=nop;
    end
end
npoin_max=max(npoin_dg);
ndof_min=min(ndof_hdg);
ndof_max=max(ndof_dg);
R=npoin_dg./npoin_cg

%Plot Complexity
h=figure;
figure(h)
plot(nop_vector,npoin_cg,'r-','Linewidth',2);
hold on;
plot(nop_vector,npoin_dg,'b-','Linewidth',2);
plot(nop_vector,npoin_hdg,'k-','Linewidth',2);
xlabel('N','FontSize',18);
ylabel('N_p','FontSize',18);
legend('CGc','DG','HDG');
axis([nop_min nop_max 0 npoin_max]); 
set(gca, 'FontSize', 18);

%Plot Operation Count
h=figure;
figure(h)
semilogy(nop_vector,ndof_cg,'r-','Linewidth',2);
hold on;
semilogy(nop_vector,ndof_dg,'b-','Linewidth',2);
semilogy(nop_vector,ndof_hdg,'k-','Linewidth',2);
xlabel('N','FontSize',18);
ylabel('Floating Point Operations','FontSize',18);
legend('CGc','DG','HDG');
axis([nop_min nop_max 0 ndof_max]); 
set(gca, 'FontSize', 18);