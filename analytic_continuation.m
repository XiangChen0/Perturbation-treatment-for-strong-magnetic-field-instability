kp=1.2;
omega=0.1+0.01i;

dkx=0.001;
dky=0.001;
k0=-10:dkx:10;



figure;
hold on;

f=gm(kp,omega,abs(k0));
f1=(circshift(f,-1)-circshift(f,1))/2/dkx;

N=50;
for ii=1:N
    
    plot(k0,imag(f));
    g=f-1i*dky*f1;
    g1=(circshift(g,-1)-circshift(g,1))/2/dkx;
    f=g;
    f1=g1;

end
grid on;
xlim([-5 5]);

%{
figure;
plot(k0,real(g));
hold on;
plot(k0,imag(g));
grid on;
xlabel('k');
legend('Re[f]','Im[f]');


figure;
plot(k0,real(f));
hold on;
plot(k0,real(f1));
%}