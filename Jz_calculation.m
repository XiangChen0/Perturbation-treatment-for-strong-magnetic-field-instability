Jz=x;
for ii=1:length(Jz)
    Jz(ii)=passINT(omega,psi,delta,vpl,tp,x(ii),xin)+trapINT(omega,psi,delta,vpl,tt,x(ii));
end
figure;
plot(x,Jz);