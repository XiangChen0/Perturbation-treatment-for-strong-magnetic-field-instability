%Calculation of tau(z,w) for passing particles(w>0) for any certain initial
%position z0>0;
function t=logtaup(psi,delta,zarray)

%Discretization of w(>0)
warray=zeros(9,12);

for ii=15:-1:4
    warray(:,16-ii)=psi*10^(-ii):psi*10^(-ii):9*psi*10^(-ii);
end

warray=reshape(warray,1,length(warray(1,:))*length(warray(:,1)));

wex=psi*1e-3:psi*1e-3:5*psi;

warray=[warray wex];

%The grid points in [w,z] space where tau is found
[wmesh,zmesh]=meshgrid(warray,zarray);

%The integrand: 1/vz
vzinv=1./sqrt(2*(wmesh+HPTL(psi,delta,zmesh)));

%Find out tau(z,w)=int{dz/vz}
tau=wmesh;

for ii=1:length(warray)
    
    tau(:,ii)=cumtrapz(zmesh(:,ii),vzinv(:,ii));
    
end

%Saving up w,z together with tau into the matrix t
t=zeros(length(tau(:,1))+1,length(tau(1,:))+1);

t(2:end,1)=zarray;

t(1,2:end)=log(warray);

t(2:end,2:end)=log(tau+1);
