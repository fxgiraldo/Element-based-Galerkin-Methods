%---------------------------------------------------------------------%
%Written by F.X. Giraldo
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [flux_q] = rusanov_flux(jac_side,nx,ny,ql,qr,ul,vl,ngl,is);

    flux_q = zeros(ngl,1);
       for l=1:ngl
              wq=jac_side(is,l);
              
              %Store Normal Vectors
              nxl=nx(is,l);
              nyl=ny(is,l);

              nxr=-nxl;
              nyr=-nyl;

              %Interpolate onto Quadrature Points
              qlq_k=ql(l);
              qrq_k=qr(l);
              u_k=ul(l);
              v_k=vl(l);

              ul_k=u_k; vl_k=v_k;
              ur_k=u_k; vr_k=v_k;

              %Compute Rusanov flux Constant
              %parent
              unl=nxl*ul_k + nyl*vl_k;
              unr=nxl*ur_k + nyl*vr_k;
              claml=abs(unl);
              clamr=abs(unr);
              clam=max(claml,clamr);
            
              %Flux Variables
              fxl=qlq_k*ul_k;
              fyl=qlq_k*vl_k;

              fxr=qrq_k*ur_k;
              fyr=qrq_k*vr_k;

              %Normal Flux Component
              flux_ql=( nxl*fxl + nyl*fyl );
              flux_qr=( nxr*fxr + nyr*fyr );
              flux_ql=flux_ql - flux_qr;

              %Dissipation Term
              diss_q=clam*(qrq_k - qlq_k);

              %Construct Rusanov Flux
              flux_q(l)=0.5*(flux_ql - diss_q);
              
          end %l   

end
