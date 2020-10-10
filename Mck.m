%continuous eigenfunction
function y=Mck(k,x)

P=-15*tanh(x)*(k.^4+(28*sech(x)^2-15)*k.^2+63*sech(x)^4-56*sech(x)^2+8);

Q=k.^4+(105*sech(x)^2-85)*k.^2+945*sech(x)^4-1155*sech(x)^2+274;

y=(P+1i*k.*Q).*exp(1i*k*x);