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

function [u,u1,u2]=scatter_gather_1(uin,u1in,u2in,ngl,iter,is,jac)

    NP = ngl;
%     if(mod(iter,2))
        NPx = 2*NP-3;
%     else
%         NPx = 2*NP-2;
%     end

%Compute LGL Points
    [xgl,wgl]=legendre_gauss_lobatto(NP);

    %Compute Legendre Cardinal functions and derivatives
    [psi,dpsi,xnq,wq] = lagrange_basis(NP,NP,xgl);

    %Compute LGL Points
    [xgl,wgl]=legendre_gauss_lobatto(NPx);

    %Compute Legendre Cardinal functions and derivatives
    [psi,dpsi,xnq,wqx] = lagrange_basis(NPx,NPx,xgl);
    
    
    
    xi = legendre_gauss_lobatto(NP); %get the gll points (possibly stored)
    xix = legendre_gauss_lobatto(NPx); %higher order base
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
    u1(1) = uina(1);
    u2(ngl) = uin(ngl);
    
    %gather
    
    for i=1:NPx
        if(xix(i)<=0.0) % segregate parent gll points corresponding to each child
            xi1l(i) = xix(i); 
            j=i;
        else
            xi2l(i-j) = xix(i);
        end
    end

    L1l = lagrange_poly(xi1l,xi1); %compute projection matrices (can store them somewhere, they will be the same for all gathers)
    L2l = lagrange_poly(xi2l,xi2);
    
    for i=1:ngl
       u1ina(i) = u1in(i);
       u2ina(i) = u2in(i); 
    end
    
    w1l = mtimes(u1ina,L1l); % project
    w2l = mtimes(u2ina,L2l);

    for i=1:j
        ua(i) = w1l(i); %construct one output vector
    end
    for i=j+1:NPx
        ua(i) = w2l(i-j);
    end
    
    ss1 = 0;
    for i=1:NPx
       ss1 = ss1 + wqx(i)*ua(i);
    end
    ss1a=0;
    for i=1:ngl
       ss1a = ss1a + wq(i)/2 *(u1in(i)+u2in(i));
    end
    ss1a=ss1-ss1a;
    
    %last projection (can be integral, if needed)
    L = lagrange_poly(xi,xix);
    u = mtimes(ua,L);
    u(1) = u1in(1);
    u(ngl) = u2in(ngl);
    
    
end