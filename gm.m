%Gamma(k)
function y=gm(k, omega, x)

    y=x;
    
    xe0=x(abs(x)<=1e-3);
    
    ye0=k^2+(1-1/omega^2)*xe0.^2-3*xe0.^4/16/omega^4;
    
    y(abs(x)<=1e-3)=ye0;
    
    
    xn0=x(abs(x)>1e-3);
    
    vn0=omega*2*sqrt(2)./xn0;
    
    
    %ind=find(imag(vn0)<=0);
    
    %alpha=imag(xn0(ind))*real(omega)./real(xn0(ind))/imag(omega);
    %we(ind)=we1(ind).*alpha+we2(ind).*(1-alpha);
    
    %vn0(ind)=conj(vn0(ind));
    we=vn0.*zetaf(vn0);
    %we(ind)=conj(we(ind));
    
    
    yn0=k^2+16+xn0.^2+16*we;
    
    y(abs(x)>1e-3)=yn0;
    
end


