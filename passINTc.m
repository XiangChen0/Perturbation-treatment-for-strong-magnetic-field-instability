%integral for passing particles

function Izp=passINTc(omega,psi,delta,vpl,tp,z0,zm,z)
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
 
 tauzmn=interp2(wmesh,zmesh,tp(2:end,2:end),logwp,-zm);
 
 tauzmp=interp2(wmesh,zmesh,tp(2:end,2:end),logwp,+zm);
 
 Dtau=tauzmp-tauzmn;
 
 tauz=interp2(wmesh,zmesh,tp(2:end,2:end),logwp,z);
 
 tau=abs(exp(tauz0)-exp(tauz));
 
 vinf=sqrt(2*wp);
 
 intgdp=wp;
 
 for ii=1:length(wp)
     
     %integration in the range (-zm,+zm)
     intgtau=1i*omega*dgdw(psi,delta,vpl,wp(ii))*exp(1i*omega*tau(:,ii)) ...
             .*Mck(kw,z)./sqrt(2*wp(ii)+2*phiz);
         
     intgdp(ii)=trapz(z,intgtau);
     
     %integration in the range (-\infty, -zm) + (+zm,+\infty) without the
     %oscillating part
     if z0<-zm
         
         intgdp(ii)=intgdp(ii)+omega/(omega+kw*vinf(ii)/4)*exp(-1i*omega/vinf(ii)*(z0+zm)) ...
             *(Mck(kw,-zm)-Mck(kw,+zm)*exp(1i*omega*Dtau));
         
     elseif z0>zm
         
         intgdp(ii)=intgdp(ii)+omega/(-omega+kw*vinf(ii)/4)*exp(1i*omega/vinf(ii)*(z0-zm)) ...
             *(Mck(kw,-zm)*exp(1i*omega*Dtau)-Mck(kw,+zm));
         
     else
         
         intgdp(ii)=intgdp(ii) ...
             +omega/(-omega+kw*vinf(ii)/4)*exp(1i*omega*(tauz0(ii)-tauzmn(ii)))*Mck(kw,-zm) ...
             -omega/(+omega+kw*vinf(ii)/4)*exp(1i*omega*(tauzmp(ii)-tauz0(ii)))*Mck(kw,+zm) ...
             -2*omega^2/(kw^2*vinf(ii)^2/16-omega^2)*Mck(kw,z0);
         
     end
     
 end
 %The integral Izp=\int{intgd*dfwp}
 
 Izp=trapz(vz0,intgdp);

 toc;
 