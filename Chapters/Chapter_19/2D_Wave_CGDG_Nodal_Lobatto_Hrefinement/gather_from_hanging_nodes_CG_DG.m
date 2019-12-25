%---------------------------------------------------------------------%
%Written by M.A. Kopera
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs_continuous = gather_from_hanging_nodes_CG_DG(rhs,rhs_continuous,intma,face,ffc,nffc,ngl)


%go over all faces and find the hanging nodes

for ifa=1:nffc
    iface = ffc(ifa);
    
    plside = 0;
   if face(iface,9)==1 && face(iface,7)==0
       parent = face(iface,5);%identify parent element
       plside = face(iface,3);%parent local face
       child1 = face(iface,6);%identify children elements 
       child2 = face(iface,8);
   elseif face(iface,9)==1 && face(iface,8)==0
        parent = face(iface,6);
        plside = face(iface,4);
        child1 = face(iface,5);
        child2 = face(iface,7);
   else
       continue
   end
   
   parin = zeros(1,ngl);
   ch1out = zeros(1,ngl);
   ch2out = zeros(1,ngl);
      
   %locate parent elements in rhs_continuous and copy data from children in
   %rhs
   switch plside
%        case 0
%            continue
       case 1
           for j=1:ngl
               parloc(j) = intma(parent,j,1);
               ch1out(j) = rhs(child1,j,ngl);
               ch2out(j) = rhs(child2,j,ngl);
           end
       case 2
           for j=1:ngl
               parloc(j) = intma(parent,ngl,j);
               ch1out(j) = rhs(child1,1,j);
               ch2out(j) = rhs(child2,1,j);
           end
       case 3           
           for j=1:ngl
               parloc(j) = intma(parent,j,ngl);
               ch1out(j) = rhs(child1,j,1);
               ch2out(j) = rhs(child2,j,1);
           end
       case 4
           for j=1:ngl
               parloc(j) = intma(parent,1,j);
               ch1out(j) = rhs(child1,ngl,j);
               ch2out(j) = rhs(child2,ngl,j);
           end
   end
   
   
   
   %gather
    parin = gather_from_children(ch1out,ch2out);
   
%    xgl = legendre_gauss_lobatto(ngl);
%    figure;
%    plot(xgl,parin,'o',(xgl-1)/2,ch1out,'gd',(xgl+1)/2,ch2out,'rx')
%    hold on;
   %distribute the data and perform DSS
    for j=2:ngl-1
       rhs_continuous(parloc(j)) = rhs_continuous(parloc(j)) + parin(j);
%        plot(xgl(j),rhs_continuous(parloc(j)),'s')
    end
%     hold off;
end

end
