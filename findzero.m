%This function solves the exact value of zeros of Gamma(k) in complex plane using gradient
%descent method. k0 is the intial guess. We know Gamma(k)=0 has totally 4
%solutions, their approximate values are +kw=kp*omega, -kw=-kp*omega, +4i
%and -4i. Taking k0 to be these approximate value, it takes very few steps
%to converges to the exact zeros.
%when kp=1.2 and omega=0.1+0.01i, the four zeros are 
%+0.103584037047468+0.011589315395968i--> +0.103553346741586 + 0.011647821758705i
%-0.103833345959402-0.009260440366076i--> -0.103553340293473 - 0.011647842301169i
%-0.022706912535874+4.404662136373913i--> +0.023269445366375 + 3.944783882416238i
%-0.023269455238880-3.944783892658416i--> -0.023269455238880 - 3.944783892658416i
function y=findzero(k0,kp,omega)

x=[real(k0),imag(k0)];
dx=1e-8;

n=20;
S=zeros(n,2);
lambda=zeros(n,1);

for jj=1:n
    
    f=gm(kp,omega,x(1)+x(2)*1i);
    
    if(abs(f)<1e-6)
        break;
    end
    
    f1=gm(kp,omega,x(1)+x(2)*1i+dx);
    f2=gm(kp,omega,x(1)+x(2)*1i+dx*1i);

    S(jj,1)=(abs(f1^2)-abs(f^2))/dx;
    S(jj,2)=(abs(f2^2)-abs(f^2))/dx;
    
    if jj==1
        lambda(jj)=1e-6;
    else
        lambda(jj)=(x-xp)*(S(jj,:)-S(jj-1,:))'/(norm(S(jj,:)-S(jj-1,:))^2);
    end
    
    xp=x;
    
    x=x-lambda(jj)*S(jj,:);
    
    
end

y=x(1)+x(2)*1i;

end