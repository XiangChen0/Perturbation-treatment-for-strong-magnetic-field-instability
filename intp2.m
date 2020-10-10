function tau=intp2(tt,alphaw,alphaz)
wsp=tt(1,3)-tt(1,2);
zsp=tt(5,1)-tt(4,1);

ccf=floor(alphaw/wsp)+2;
ccf(ccf>=length(tt(1,:)))=length(tt(1,:))-1;
ccc=ccf+1;

rrf=floor(alphaz/zsp)+4;
rrf(rrf>=length(tt(:,1)))=length(tt(:,1))-1;
rrc=rrf+1;

M=length(tt(:,1));
tau=((tt(1,ccc)-alphaw).*(tt(rrc,1)'-alphaz).*tt(ccf*M+rrf-M)+ ...
    (tt(1,ccc)-alphaw).*(alphaz-tt(rrf,1)').*tt(ccf*M+rrc-M)+ ...
    (alphaw-tt(1,ccf)).*(tt(rrc,1)'-alphaz).*tt(ccc*M+rrf-M)+ ...
    (alphaw-tt(1,ccf)).*(alphaz-tt(rrf,1)').*tt(ccc*M+rrc-M))/ ...
    zsp/wsp;

