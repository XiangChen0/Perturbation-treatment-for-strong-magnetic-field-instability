tic;

kp=1.2;
omega=0.1+0.01i;

%direct integration
k1=0:0.005:10;
k2=10:0.1:400;
k=[k1 k2];

z=1;

z0=(-30:0.1:30)/4;
M=z0;

for ii=1:length(z0)

%modified plasma . dispersion function
intg=2*real(Mck(k,z).*Mck(-k,z0(ii)))./gm(kp,omega,k)./(k.^2+1)./(k.^2+4)./(k.^2+9)./(k.^2+16)./(k.^2+25);

M(ii)=trapz(k,intg);

end

figure;
plot(z0,real(M));
hold on;
plot(z0,imag(M));
grid on;
xlim([min(z0),max(z0)]);
toc;

%{
%standard plasma dispersion function
intg=Mck(k,z).*Mck(-k,z0(ii))./gm(kp,omega,k)./(k.^2+1)./(k.^2+4)./(k.^2+9)./(k.^2+16)./(k.^2+25);
intg(49001)=(intg(49000)+intg(49002))/2;
M1(ii)=trapz(k,intg);
%}