function y=dgdw(psi,delta,vpl,w)

indp=find(w>=0);
indn=find(w<0);

%dgdw for passing particles
yp=-1/sqrt(2*pi)/vpl^3*exp(-w(indp)/vpl^2);
   

%dgdw for trapped particles
yn=15*sqrt(2)/2/delta^2/sqrt(psi)-1/sqrt(2*pi)/vpl^3*exp(-w(indn)/vpl^2).*erfc(sqrt(-w(indn))/vpl);
    
y=w;

y(indp)=yp;
y(indn)=yn;