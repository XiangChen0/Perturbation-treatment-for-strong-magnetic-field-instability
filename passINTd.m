%integral for passing particles

function Izp=passINTd(omega,psi,delta,vpl,tp,z0,z)
 tic;

 z=reshape(z,length(z),1);
 
 
 %-------------------------------------------------------------------------
 %Find out the indices of z0 and z in zarray and the time requires to go
 %from z0 to z or from z to z0 for different energy wp

 phiz0=HPTL(psi,delta,z0);
 
 phiz=HPTL(psi,delta,z);
 
 vz0min=sqrt(2*phiz0+2*exp(tp(1,2)));
 
 vz0max=sqrt(2*phiz0+2*exp(tp(1,end)));
 
 vz0=vz0min+(vz0max-vz0min)*exp(-15:1e-2:0);
 
 wp=vz0.^2/2-phiz0;
 
 logwp=log(wp);
 
 [wmesh,zmesh]=meshgrid(tp(1,2:end),tp(2:end,1));
 
 %interpolation scheme of tau. Here we take advantage that
 %tau~wp^(1/2), so log(tau) is nearly in proportion with log(wp)
 tauz0=interp2(wmesh,zmesh,tp(2:end,2:end),logwp,z0);
 
 tauz=interp2(wmesh,zmesh,tp(2:end,2:end),logwp,z);
 
 tau=abs(exp(tauz0)-exp(tauz));
 
 intgdp=wp;
 
 for ii=1:length(wp)
     
     intgtau=1i*omega*dgdw(psi,delta,vpl,wp(ii))*exp(1i*omega*tau(:,ii)) ...
             .*HPTL1(psi,delta,z)./sqrt(2*wp(ii)+2*phiz);
     
     intgdp(ii)=trapz(z,intgtau);

     
 end
 %The integral Izp=\int{intgd*dfwp}
 
 Izp=trapz(vz0,intgdp);

 toc;
 