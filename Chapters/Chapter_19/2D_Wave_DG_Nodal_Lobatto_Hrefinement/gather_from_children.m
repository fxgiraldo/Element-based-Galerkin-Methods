%---------------------------------------------------------------------%
%This code gathers coefficients from two children faces to the parent face
% u1,u2 - vector of coefficients of children edges
%         the convention is that u1 corresponds to [-1,0] child
%         and u2 corresponds to [0,1] child
% u     - output vector of coefficients at the parent edge 
%
%Written by M.A. Kopera on 10/2011
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function u = gather_from_children(u1,u2)

    NP = size(u1,2);
    
    xi = legendre_gauss_lobatto(NP); %get the gll points (possibly stored)
    xi1 = xi*0.5-0.5; % get children gll points
    xi2 = xi*0.5+0.5;   
    
    for i=1:NP
        if(xi(i)<=0.0) % segregate parent gll points corresponding to each child
            xi1l(i) = xi(i); 
            j=i;
        else
            xi2l(i-j) = xi(i);
        end
    end

    L1l = lagrange_poly(xi1l,xi1); %compute projection matrices (can store them somewhere, they will be the same for all gathers)
    L2l = lagrange_poly(xi2l,xi2);
    w1l = mtimes(u1,L1l); % project
    w2l = mtimes(u2,L2l);

    for i=1:j
        u(i) = w1l(i); %construct one output vector
    end
    for i=j+1:NP
        u(i) = w2l(i-j);
    end
end