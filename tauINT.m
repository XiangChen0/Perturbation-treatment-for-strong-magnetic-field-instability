%Time integral
function hWz=tauINT(omega,psi,delta,Wz,z0,phi1z)
%passing particle
if Wz>=0
   %z=-32*delta:(1e-3)*delta:32*delta;
   vz_1=1./sqrt(Wz+2*HPTL(psi,delta,phi1z(1,:)));
   tz1=cumtrapz(phi1z(1,:),vz_1);
   [~,n0]=min(abs(phi1z(1,:)-z0));
   tz=abs(tz1-tz1(n0));
   ftz=exp(1i*omega*tz);
   hWz=-trapz(ftz(1:n0),phi1z(2,1:n0))+trapz(ftz(n0+1:end),phi1z(2,n0+1:end));
%trapped particle
else
    zmax=sqrt(-log(-Wz/2/psi))*sqrt(2)*delta;
    dx=phi1z(1,2)-phi1z(1,1);
    zmax=dx*floor(zmax/dx);
    z=-zmax:dx:zmax;
    z(1)=z(1)+(1e-4)*delta;
    z(end)=z(end)-(1e-4)*delta;
    vz_1=1./sqrt(Wz+2*HPTL(psi,delta,z));
    tz1=cumtrapz(z,vz_1);
    T=(tz1(end)-tz1(1))*2;
    %disp([2*pi/T,zmax]);
    tz2=tz1+T/2;
    tz=[tz1,tz2];
    [~,n0]=min(abs(z-z0));
    tz=abs(tz-tz(n0));
    tz=mod(tz,T);
    mn=length(phi1z(2,:))/2+1/2-zmax/dx;
    mp=length(phi1z(2,:))/2+1/2+zmax/dx;
    phi1trap=phi1z(2,mn:mp);
    phi1trap=[phi1trap,fliplr(phi1trap)];
    hWz=-2/T*trapz(tz,phi1trap);
    for jj=1:10
        sintz=sin(jj*2*pi/T*tz);
        hWz=hWz+2/pi/jj*omega^2/(4*pi^2*jj^2/T^2-omega^2)*trapz(sintz,phi1trap);
    end
end