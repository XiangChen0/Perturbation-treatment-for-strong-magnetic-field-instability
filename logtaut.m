%Calculation of tau(z,w) for trapped particles(w<0) for any certain initial
%position z0>0;

function t=logtaut(psi,delta,alphaw,alphaz)

%Discretization (w,z) space (<w<0,-inf<z<+inf)
%zarray=zarray(zarray>=0);
%{
warray=zeros(9,15+log10(alphaw));
for ii=15:-1:1-log10(alphaw)
    warray(:,16-ii)=psi*10^(-ii):psi*10^(-ii):9*psi*10^(-ii);
end
warray=reshape(warray,1,length(warray(1,:))*length(warray(:,1)));
warray=-fliplr(warray);
wex=-psi*(1-alphaw):psi*alphaw:min(warray);
warray=[wex warray];
%}

wspc=0:-alphaw:-36;

wmin=-psi*(1-1e-6);

wmax=0;

%w1=wmin+(wmax-wmin)/2*exp(-10:alphaw:0);
warray=wmax-(wmax-wmin)*exp(wspc);

%warray=[w1 w2];

zmax=acosh((-warray/psi).^(-1/4))*delta;

zspc=0:alphaz:1-alphaz;

zmesh=zspc'*zmax;

wmesh=ones(length(zspc),1)*warray;

%The inverse of vz which is the integrand of tau-integral
vzinv=1./sqrt(2*(wmesh+HPTL(psi,delta,zmesh)));

%Find out tau(z,w)
tau=wmesh;

for ii=1:length(warray)
    
    tau(:,ii)=cumtrapz(zmesh(:,ii),vzinv(:,ii));
    
end

%First order correction to the sojourn time near the turning points
taucr=sqrt(2*alphaz*zmax./abs(HPTL1(psi,delta,zmax)));

%Quarter period of each orbit with energy w<0
Tqt=(tau(end,:)-tau(1,:))+taucr;

%Saving up w,z together with tau into matrix t

tau=[log(tau+1);log(Tqt+1)];

t=zeros(length(tau(:,1))+3,length(tau(1,:))+1);

t(1,2:end)=wspc;

t(2,2:end)=warray;

t(3,2:end)=zmax;

t(4:end,1)=(0:alphaz:1)';

t(4:end,2:end)=tau;

