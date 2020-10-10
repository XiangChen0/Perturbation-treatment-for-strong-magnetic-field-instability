function y=HPTL1(psi,delta,z)

y=-4*psi/delta*tanh(z/delta).*(sech(z/delta)).^4;