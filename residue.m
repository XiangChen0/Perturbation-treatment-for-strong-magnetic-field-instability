%residues calculation

if z<=z0
    
    kw=-0.103833345959402-0.009260440366076i;
    ks=-0.023269455238880-3.944783892658416i;
    k1=-1i;k2=-2i;k3=-3i;k4=-4i;k5=-5i;
    wind=-1;
else
    
    kw=+0.103584037047468+0.011589315395968i;
    ks=-0.022706912535874+4.404662136373913i;
    k1=+1i;k2=+2i;k3=+3i;k4=+4i;k5=+5i;
    wind=+1;
end

%residues
Rs=Mck(ks,z)*Mck(-ks,z0)/gm1(kp,omega,ks)/(ks^2+1)/(ks^2+4)/(ks^2+9)/(ks^2+16)/(ks^2+25);
Rw=Mck(kw,z)*Mck(-kw,z0)/gm1(kp,omega,kw)/(kw^2+1)/(kw^2+4)/(kw^2+9)/(kw^2+16)/(kw^2+25);

R1=Mck(k1,z)*Mck(-k1,z0)/gm(kp,omega,k1)/(2*k1)/(k1^2+4)/(k1^2+9)/(k1^2+16)/(k1^2+25);
R2=Mck(k2,z)*Mck(-k2,z0)/gm(kp,omega,k2)/(k2^2+1)/(2*k2)/(k2^2+9)/(k2^2+16)/(k2^2+25);
R3=Mck(k3,z)*Mck(-k3,z0)/gm(kp,omega,k3)/(k3^2+1)/(k3^2+4)/(2*k3)/(k3^2+16)/(k3^2+25);
R4=Mck(k4,z)*Mck(-k4,z0)/gm(kp,omega,k4)/(k4^2+1)/(k4^2+4)/(k4^2+9)/(2*k4)/(k4^2+25);
R5=Mck(k5,z)*Mck(-k5,z0)/gm(kp,omega,k5)/(k5^2+1)/(k5^2+4)/(k5^2+9)/(k5^2+16)/(2*k5);

%winding number -1 if enclose the half lower plane
M=wind*2i*pi*(Rs+Rw+R1+R2+R3+R4+R5);
%------------------------------------------------


%direct integration
k1=-400:0.01:-10;
k2=-10:0.001:10;
k3=10:0.01:400;
k=[k1 k2 k3];
intg=Mck(k,z).*Mck(-k,z0)./gm(kp,omega,abs(k))./(k.^2+1)./(k.^2+4)./(k.^2+9)./(k.^2+16)./(k.^2+25);
intg(49001)=(intg(49000)+intg(49002))/2;
M1=trapz(k,intg);

%add=-64i*sqrt(2*pi)*omega./k.*exp(-8*omega^2./k.^2);



