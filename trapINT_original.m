function Izt=trapINT(omega,psi,delta,tt,dz,z0,z)

phiz0=HPTL(psi,delta,z0);
phizf=HPTL(psi,delta,abs(z)+dz/2);
phizb=HPTL(psi,delta,abs(z)-dz/2);
%--------------------------------------------------------------------------           
vz0max=sqrt(2*(2*tt(2,end)+phiz0));

wmid=2*max(phiz0-phizf,0);
vz0mid=sqrt(wmid);

wmin=2*max(phiz0-phizb,0);
vz0min=sqrt(wmin);
if z0==z && z0==0
     vz0min=sqrt(2*(tt(2,2)+phiz0));
end 

if vz0min>=vz0max || vz0mid>=vz0max
    Izt=0;
    return;
end
%spacing
if vz0mid==vz0min
   vz01=vz0min+(vz0max-vz0min)/2*exp(-15:1e-2:0);
   vz02=vz0max-(vz0max-vz0min)/2*exp(0:-1e-2:-15);
   vz0=[vz01 vz02];
else
    dvz01=(vz0mid-vz0min)*1e-3;
    vz01=vz0min:dvz01:vz0mid;
    vz02=vz0mid+(vz0max-vz0mid)/2*exp(-15:1e-2:0);
    vz03=vz0max-(vz0max-vz0mid)/2*exp(0:-1e-2:-15);
    vz0=[vz01 vz02 vz03];
end
%--------------------------------------------------------------------------
wt=vz0.^2/2-phiz0;
alphaw=log((tt(2,end)-wt)/(tt(2,end)-tt(2,2)));

zmax=acosh((-wt/psi).^(-1/4))*delta;
zmax(zmax==0)=1e-9;
alphaz0=z0./zmax;
alphazf=(abs(z)+dz/2)./zmax;
alphazf(alphazf>1)=1; 
tauz0=intp2(tt,alphaw,alphaz0);
tauzf=intp2(tt,alphaw,alphazf);
Tqt=intp2(tt,alphaw,1); 
tauz0=exp(tauz0)-1;
tauzf=exp(tauzf)-1;
Tqt=exp(Tqt)-1; 

if z~=0  
   alphazb=(abs(z)-dz/2)./zmax; 
   tauzb=intp2(tt,alphaw,alphazb); 
   tauzb=exp(tauzb)-1; 
   tau1=abs(tauzf-tauzb);
end
           
dgdwt=dgdw(psi,delta,vpl,wt);

if z0==z && z0==0                            
       tau1=2*tauzf;    
       tau2=2*Tqt-tauzf;
       tau3=4*Tqt-tauzf;   
       intgdt=2*dgdwt.*((exp(1i*omega*tau1)-1)./(1-exp(1i*omega*4*Tqt)) ...
              .*(exp(1i*omega*tau2)+exp(1i*omega*tau3))+exp(1i*omega*tauzf)-1);       
elseif z0==z       
       tau2=4*Tqt+tauzb-tauz0;
       tau3=2*Tqt-tauzf-tauz0;
       tau4=2*Tqt+tauzb+tauz0;
       tau5=4*Tqt-tauzf+tauz0;        
       intgdt=dgdwt.*(exp(1i*omega*tau1)-1).*(exp(1i*omega*tau2)+exp(1i*omega*tau3) ...
              +exp(1i*omega*tau4)+exp(1i*omega*tau5))./(1-exp(1i*omega*4*Tqt)) ...
              +dgdwt.*(exp(1i*omega*(tauzf-tauz0))+exp(1i*omega*(tauz0-tauzb))-2);            
elseif z==0  
       tau1=2*tauzf;
       tau2=tauz0-tauzf;    
       tau3=2*Tqt-tauzf+tauz0;
       tau4=2*Tqt-tauzf-tauz0;
       tau5=4*Tqt-tauzf-tauz0;          
       intgdt=dgdwt.*(exp(1i*omega*tau1)-1).*(exp(1i*omega*tau2)+exp(1i*omega*tau3) ...
              +exp(1i*omega*tau4)+exp(1i*omega*tau5))./(1-exp(1i*omega*4*Tqt)); 
          
elseif z>0 && z<z0
       tau2=tauz0-tauzf;    
       tau3=2*Tqt+tauzb+tauz0;
       tau4=2*Tqt-tauzf-tauz0;
       tau5=4*Tqt+tauzb-tauz0;          
       intgdt=dgdwt.*(exp(1i*omega*tau1)-1).*(exp(1i*omega*tau2)+exp(1i*omega*tau3) ...
              +exp(1i*omega*tau4)+exp(1i*omega*tau5))./(1-exp(1i*omega*4*Tqt)); 
          
elseif z<0    
       tau2=tauzb+tauz0;
       tau3=2*Tqt-tauzf+tauz0;
       tau4=2*Tqt+tauzb-tauz0;
       tau5=4*Tqt-tauzf-tauz0;                            
       intgdt=dgdwt.*(exp(1i*omega*tau1)-1).*(exp(1i*omega*tau2)+exp(1i*omega*tau3) ...
              +exp(1i*omega*tau4)+exp(1i*omega*tau5))./(1-exp(1i*omega*4*Tqt)); 
          
else   
       tau2=tauzb-tauz0;
       tau3=2*Tqt-tauzf-tauz0;
       tau4=2*Tqt+tauzb+tauz0;
       tau5=4*Tqt-tauzf+tauz0;                            
       intgdt=dgdwt.*(exp(1i*omega*tau1)-1).*(exp(1i*omega*tau2)+exp(1i*omega*tau3) ...
              +exp(1i*omega*tau4)+exp(1i*omega*tau5))./(1-exp(1i*omega*4*Tqt));
          
end
%plot(vz0,real(intgdt),'.');
Izt=trapz(vz0,intgdt);
