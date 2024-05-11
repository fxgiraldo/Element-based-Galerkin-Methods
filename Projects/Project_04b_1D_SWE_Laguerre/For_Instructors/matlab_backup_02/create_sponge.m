%---------------------------------------------------------------------%
%This function computes the mass and energy.
%Written by F.X. Giraldo on January 19, 2024.
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [Igamma,Igamma_volume,xmax_LGL,xt,xd]=create_sponge(coord,coord_cg,npoin,ngl_LGL,nelem_LGL,lsponge,sponge_method,sponge_shape,sponge_amp,sponge_exponent)

%Klemp-Lily sponge
xt = max(coord_cg);
xd = coord(ngl_LGL,nelem_LGL);
ds = max(coord_cg-xd,0);
if (sponge_shape==1)
    sponge_function = sponge_amp.*( sin(pi/2 .* ds./(xt -xd)) ).^sponge_exponent;
elseif (sponge_shape==2)
    sponge_function = sponge_amp.*( tanh(ds./(xt -xd)) ).^sponge_exponent;
end
Igamma=ones(npoin,1);
Igamma_volume=zeros(npoin,1);
if (lsponge==1)
    if (sponge_method==1) %Hard reset
        Igamma=1.0 - sponge_function;
    elseif (sponge_method==2) %gentle in CREATE_RHS_VOLUME
        Igamma_volume=1.0 - sponge_function;
    end
end
disp(['xmax =  ',num2str(xt),' xd = ', num2str(xd)]);
xmax_LGL=xd;
