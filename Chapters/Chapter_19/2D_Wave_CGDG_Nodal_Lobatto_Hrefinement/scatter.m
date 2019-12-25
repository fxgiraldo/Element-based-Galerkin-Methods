%---------------------------------------------------------------------%
%This code scatters the coefficients from parent face to children faces.
% u1,u2 - output vector of coefficients of children edges
%         the convention is that u1 corresponds to [-1,0] child
%         and u2 corresponds to [0,1] child
% u     - input vector of coefficients at the parent edge 
%
%Written by M.A. Kopera on 10/2011
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%

function [u1,u2]=scatter(uin,ngl)

    NP = ngl;
    
    xi = legendre_gauss_lobatto(NP); %get the gll points (possibly stored)
    xi1 = xi*0.5-0.5; % get children gll points
    xi2 = xi*0.5+0.5;   
    uina = zeros(1,ngl);
    for i=1:ngl
        uina(i) = uin(i);
    end
    L1 = lagrange_poly(xi1,xi); %compute projection matrices for children
    L2 = lagrange_poly(xi2,xi);

    u1 = mtimes(uina,L1); %project data onto children
    u2 = mtimes(uina,L2);
    
    
end