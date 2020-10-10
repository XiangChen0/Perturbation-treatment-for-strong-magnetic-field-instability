Iz=zeros(length(x),1);
%figure;
%hold on;
tic;
for jj=1:length(x)
    Iz(jj)=trapINT(omega,psi,delta,dpdft,tt,dx,x(500),x(jj));
    %WzINT(omega,psi,delta,vpl,dpdft,tp,tt,wr,dx,x(450),x(jj));
    disp(jj);
end
%figure;
plot(x,real(Iz)/dx,'.');
hold on;
plot(x,imag(Iz)/dx,'.');
grid on;
hold off;
toc;

%{
tic;
z=-10:0.1:10;
Iz=z;
for ii=1:length(z)
    Iz(ii)=WzINT(0.1+0.005i,1,2.5,z(ii),phizb);
end
plot(z,real(Iz));
hold on;
plot(z,imag(Iz));
plot(phizb(1,:),phizb(2,:));
grid on;
toc;
%}

%{
%Solving the matrix inversion
load z_Iz.mat;
k=2*pi/30;
z=Iz(1,:);
Fz=Iz(2,:);
N=length(z)+2;
dx=abs(z(2)-z(1));
z=[z(1)-dx z z(end)+dx];
A=ones(N);
A=(tril(A,-1)-tril(A,-2)+tril(A,1)-tril(A,0))/dx^2;
for ii=2:N-1
    A(ii,ii)=-2/dx^2-z(ii)^2/delta^4+3/delta^2-k^2;
end
A(1,1)=-1/dx^2;
A(end,end)=-1/dx^2;
b=[0;Fz';0];
phi_new=A\b;
hold on;
plot(z,real(phi_new));
grid on;
%}



phu=x;
for ii=1:length(x)
    w=NewMode(x(ii),x/delta);
    phu(ii)=trapz(x,w.*Jz');
end
figure;plot(x,real(phu));
hold on;plot(x,imag(phu));