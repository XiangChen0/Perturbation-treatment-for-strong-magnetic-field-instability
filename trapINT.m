function Izt=trapINT(omega,psi,delta,vpl,tt,z0)
tic;

flag=1;

if z0==0
    
    Izt=0;
    
    return;
    
elseif z0<0
    
    z0=abs(z0);
    
    flag=-1;
    
end

phiz0=HPTL(psi,delta,z0);

%--------------------------------------------------------------------------           

%Spacing of energy (or vz0). Equal space is not a good way. The principle is
%to resolve the extreme cases: wt=0 which corresponds to that the particle 
%goes very far towards the infinity; vz0=0 which corresponds to that initially 
%the particle is at a turning point.
vz0max=sqrt(2*(tt(2,end)+phiz0));

vz0min=0;

%vz01=vz0min+(vz0max-vz0min)/2*exp(-15:5e-1:0);

vz0=vz0max-(vz0max-vz0min)*exp(0:-5e-2:-15);

%vz0=[vz01 vz02];

wt=vz0.^2/2-phiz0;

%find out the value of alphaw which is needed for interpolation later
alphaw=log(wt/tt(2,2));

%position of turning point
zmax=acosh((-wt/psi).^(-1/4))*delta;

zmax(zmax==0)=1e-9;

%relative position in the trajectory of a trapped particle, defined as the
%ratio of z0/zmax
alphaz0=z0./zmax;

alphaz0(alphaz0>=1)=1-1e-6;
 
[wmesh,zmesh]=meshgrid(tt(1,2:end),tt(4:end,1));

tauz0=interp2(wmesh,zmesh,tt(4:end,2:end),alphaw,alphaz0);

tauz0=exp(tauz0)-1;

Tqt=interp2(wmesh,zmesh,tt(4:end,2:end),alphaw,1); 

Tqt=exp(Tqt)-1; 

tauz=interp2(wmesh,zmesh,tt(4:end,2:end),alphaw,tt(4:end,1));

tauz=exp(tauz)-1;

tau1=abs(tauz-tauz0);

tau2=abs(tauz+tauz0+2*Tqt);

tau3=abs(tauz+tauz0);

tau4=abs(-tauz+tauz0+2*Tqt);

z=tt(4:end,1)*zmax*(1-1e-4);

vz=sqrt(2*wt+2*HPTL(psi,delta,z));

intgdt=wt;
 
 for ii=1:length(wt)
     
     intgtau=1i*omega*HPTL1(psi,delta,z(:,ii))./vz(:,ii) ...
             *dgdw(psi,delta,vpl,wt(ii))/(1-exp(4i*omega*Tqt(ii))).* ...
             (exp(1i*omega*tau1(:,ii))+exp(1i*omega*(4*Tqt(ii)-tau1(:,ii))) ...
             +exp(1i*omega*tau2(:,ii))+exp(1i*omega*(4*Tqt(ii)-tau2(:,ii))) ...
             -exp(1i*omega*tau3(:,ii))-exp(1i*omega*(4*Tqt(ii)-tau3(:,ii))) ...
             -exp(1i*omega*tau4(:,ii))-exp(1i*omega*(4*Tqt(ii)-tau4(:,ii))));
     
     intgdt(ii)=trapz(z(:,ii),intgtau);
     
 end
 
 Izt=flag*trapz(vz0,intgdt);
 
 
 %figure;plot(vz0,intgdt);
 
 toc;
 
