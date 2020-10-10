function y=PDF(psi,delta,vpl,w)

indp=find(w>=0);
indn=find(w<0);

%dgdw for passing particles
yp=1/sqrt(2*pi)/vpl*exp(-w(indp)/vpl^2);
   

%dgdw for trapped particles
yn=16*sqrt(-2*w(indn))/pi/delta^2+15*sqrt(2)*w(indn)/2/delta^2/sqrt(psi)+ ...
   1/sqrt(2*pi)/vpl*exp(-w(indn)/vpl^2).*erfc(sqrt(-w(indn))/vpl);
    
y=w;

y(indp)=yp;
y(indn)=yn;