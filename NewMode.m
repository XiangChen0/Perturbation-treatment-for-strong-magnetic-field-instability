%The combinition of all modes that could couple with MD4
function w=NewMode(kp,omega,z,x)

  
  indp=find(x>=z);
  indn=find(x<=z);
  
  kwp=findzero(kp*omega,kp,omega);
  kwn=findzero(-kp*omega,kp,omega);
  ksp=findzero(+4i,kp,omega);
  ksn=findzero(-4i,kp,omega);
  ampwp=1/gm1(kp,omega,kwp(1)+1i*kwp(2));
  ampwn=1/gm1(kp,omega,kwn(1)+1i*kwn(2));
  ampsp=1/gm1(kp,omega,ksp(1)+1i*ksp(2));
  ampsn=1/gm1(kp,omega,ksn(1)+1i*ksn(2));
  
  wp=(1/12-2*pi/gm(kp,omega,2i))*MD2(z).*MD2(x(indp)) ...
     -2*pi/gm(kp,omega,4i)*MD4(z).*MD4(x(indp)) ...
     +ampwp*exp(1i*kwp*(x(indp)-z)).*Fcal(kwp,z,x(indp)) ...
     +ampsp*exp(1i*ksp*(x(indp)-z)).*Fcal(ksp,z,x(indp));
  %wp= %wp=-2*pi/gm(kp,omega,4i)*MD4(z).*MD4(x(indp));
  wn=(1/12-2*pi/gm(kp,omega,2i))*MD2(z).*MD2(x(indn)) ...
     -2*pi/gm(kp,omega,4i)*MD4(z).*MD4(x(indn)) ...
     +ampwn*exp(1i*kwn*(z-x(indn))).*Fcal(kwn,x(indn),z) ...
     +ampsn*exp(1i*ksn*(z-x(indn))).*Fcal(ksn,x(indn),z);
  %wn=%wn=-2*pi/gm(kp,omega,4i)*MD4(z).*MD4(x(indn));
   
  w=x;

  w(indp)=wp;
  w(indn)=wn;
      
end


function y=Fcal(k, x, y)
 
    y=(Pk(k,x)+1i*k*Qk(k,x)).*(Pk(k,y)-1i*k*Qk(k,y))/(k^2+1)/(k^2+4)/(k^2+9)/ ...
        (k^2+16)/(k^2+25);
    
end

%Mck=exp(ikx)*(Pk+ikQk)
function y=Pk(k, x)

    y=-15*tanh(x).*(k^4+(28*sech(x).^2-15)*k^2+63*sech(x).^4-56*sech(x).^2+8);

end


function y=Qk(k, x)

    y=k^4+(105*sech(x).^2-85)*k^2+945*sech(x).^4-1155*sech(x).^2+274;

end

%Discrete mode MD2 with eigenvalue k=2i
function y=MD2(x)
    
     y=sqrt(105)/4*(2*sech(x).^2-3*sech(x).^4).*tanh(x);
     
end



%Discrete mode MD4 with eigenvalue k=4i
function y=MD4(x)
    
     y=3*sqrt(70)/8*sech(x).^4.*tanh(x);
     
end



%----------------------------Old version-----------------------------------
%{
%The combinition of all modes that could couple with MD4
function w=NewMode(z,x)

  indp=find(x>=-z);
  indn=find(x<=-z);

  wp=1/12*MD2(z)*MD2(x(indp))-35*pi/32*MD2(z)*MD2(x(indp))-2147*pi/512*MD4(z)*MD4(x(indp)) ...
     +pi*NM4p(z)*MD4(x(indp))+pi*NM4p(x(indp))*MD4(z)-315*pi/128*MD4(z)*(z+x(indp)).*MD4(x(indp));
 %
  wn=1/12*MD2(z)*MD2(x(indn))-35*pi/32*MD2(z)*MD2(x(indn))-2147*pi/512*MD4(z)*MD4(x(indn)) ...
     -pi*NM4n(z)*MD4(x(indn))-pi*MD4(z)*NM4n(x(indn))+315*pi/128*MD4(z)*(z+x(indn)).*MD4(x(indn));
 %
  w=x;

  w(indp)=wp;
  w(indn)=wn;
      
end

%New mode resulting from the second order poles at k=-4 (p=S+A)
function y=NM4p(x)

    y=exp(-4*x).*(945*sech(x).^4-6195*sech(x).^2+5634-60*tanh(x).*(56*sech(x).^2-94))/384;

end

%New mode resulting from the second order poles at k=-4 (n=S-A)
function y=NM4n(x)

     y=exp(4*x).*(945*sech(x).^4-6195*sech(x).^2+5634+60*tanh(x).*(56*sech(x).^2-94))/384;


end
%}
