%---------------------------------------------------------------------%
%Written by M.A. Kopera
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [fe,ip,jc,kc,js,je,ks,ke,pi,pj] = get_pointers(ngl)

    fe(1,1) = 1; %joins face side and children elements
    fe(1,2) = 2;
    fe(2,1) = 4;
    fe(2,2) = 3;
    fe(3,1) = 1;
    fe(3,2) = 4;
    fe(4,1) = 2;
    fe(4,2) = 3;
    
    ip(1) = 2;
    ip(2) = 1;
    ip(3) = 4;
    ip(4) = 3;
    
    jc(1) = 1;
    jc(2) = ngl;
    jc(3) = ngl;
    jc(4) = 1;
    
    kc(1) = 1;
    kc(2) = 1;
    kc(3) = ngl;
    kc(4) = ngl;
    
    js(1) = 1;
    js(2) = 2;
    js(3) = 1;
    js(4) = 1;
    
    je(1) = ngl;
    je(2) = ngl;
    je(3) = ngl;
    je(4) = ngl-1;
    
    ks(1) = 1;
    ks(2) = 1;
    ks(3) = 2;
    ks(4) = 2;
    
    ke(1) = ngl;
    ke(2) = ngl;
    ke(3) = ngl;
    ke(4) = ngl;
    
    pi(1) = 1;
    pi(2) = ngl;
    pi(3) = 1;
    pi(4) = 1;
    
    pj(1) = ngl;
    pj(2) = ngl;
    pj(3) = ngl;
    pj(4) = 1;
    
end
