%---------------------------------------------------------------------%
%Written by M.A. Kopera
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
for i=1:ngl
    for j=1:ngl
        qqp(i,j) = qp1(7,i,j);
        qqe(i,j) = qe(7,i,j);
        qq1(i,j) = qp1(26,i,j);
        qq2(i,j) = qp1(27,i,j);
        qq3(i,j) = qp1(28,i,j);
        qq4(i,j) = qp1(29,i,j);
        
        qqe1(i,j) = qe(26,i,j);
        qqe2(i,j) = qe(27,i,j);
        qqe3(i,j) = qe(28,i,j);
        qqe4(i,j) = qe(29,i,j);
        
        qqee1(i,j) = qe2(26,i,j);
        qqee2(i,j) = qe2(27,i,j);
        qqee3(i,j) = qe2(28,i,j);
        qqee4(i,j) = qe2(29,i,j);
        
        qqee(i,j) = qe1(7,i,j);
    end
end
% figure;
% surf(xgl,xgl,qqp);
% figure;
% surf(xgl*0.5-0.5,xgl*0.5-0.5,qq1);
% hold on;
% surf(xgl*0.5-0.5,xgl*0.5+0.5,qq2);
% surf(xgl*0.5+0.5,xgl*0.5+0.5,qq3);
% surf(xgl*0.5+0.5,xgl*0.5-0.5,qq4);

figure;
surf(xgl*0.5-0.5,xgl*0.5-0.5,qqe1);
hold on;
surf(xgl*0.5-0.5,xgl*0.5+0.5,qqe2);
surf(xgl*0.5+0.5,xgl*0.5+0.5,qqe3);
surf(xgl*0.5+0.5,xgl*0.5-0.5,qqe4);

% figure;
% surf(xgl,xgl,qqee);

figure;
surf(xgl*0.5-0.5,xgl*0.5-0.5,qqee1);
hold on;
surf(xgl*0.5-0.5,xgl*0.5+0.5,qqee2);
surf(xgl*0.5+0.5,xgl*0.5+0.5,qqee3);
surf(xgl*0.5+0.5,xgl*0.5-0.5,qqee4);
