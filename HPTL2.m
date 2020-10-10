function y=HPTL2(psi,delta,z)

y=-4*psi/delta^2*((sech(z/delta)).^6-4*(tanh(z/delta)).^2.*(sech(z/delta)).^4);