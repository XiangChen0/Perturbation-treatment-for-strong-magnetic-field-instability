tic;
%frequency of the perturbation
omega=0.08+0.004i;
%electron hole width
delta=4;
%electron hole height
psi=0.7;
%parallel thermal velocity (vpl=delta/4 to ensure dg/dw is finite at w=0)
vpl=delta/4;
%spacing
dx=1e-1;
%the interval of integration(assuming \phi_1=0 outside of this region)
xmax=30;
x=-xmax:dx:xmax;

N=length(x);
%--------------------------------------------------------------------------
%Construction of the matrix A s.t. A*\phi_1=0
A=ones(N);
A=(tril(A,-1)-tril(A,-2)+tril(A,1)-tril(A,0))/dx^2;
A(1:N+1:N^2)=-2/dx^2-(16-30*(sech(x/delta)).^2)/delta^2;
%caculation of tau for passing particles
xin=x(1)-dx:dx/2:x(end)+dx;
tp=logtaup(psi,delta,xin);
%caculation of tau for trapped particles
%spacing factor of energy for trapped particles
alphaw=1e-2;
%spacing factor of z for trapped particles
alphaz=5e-4;
tt=logtaut(psi,delta,alphaw,alphaz);

B=zeros(N);
for ii=(N+1)/2:N
    disp([ii toc]);
    for jj=1:N  
        B(ii,jj)=passINT(omega,psi,delta,vpl,tp,dx,x(ii),x(jj))+ ...
                 trapINT(omega,psi,delta,tt,dx,x(ii),x(jj));
        B(N+1-ii,N+1-jj)=B(ii,jj);
        %B(jj,ii)=B(ii,jj);
        %B(N+1-jj,N+1-ii)=B(ii,jj);
    end
end
A=A+B;
%Boundary conditions
A(1,1)=1/dx^2;
A(1,2)=-A(1,1);
A(end,end)=A(1,1);
A(end,end-1)=-A(1,1);
%--------------------------------------------------------------------------
%Find out the eigenvalue lambda that has the smallest positive real part and
%its imaginary part is negligible compared to its real part. In other words,
%the smallest and the closest to the positive real axis eigenvalue. 
%The eigenvector that the eigenvalue corresponds to is the solution we are
%looking for.
[V,D]=eig(A);
lm=diag(D);
lm(real(lm)<0)=nan;
[lmin,lind]=min(abs(lm));
plot(x,real(V(:,lind)));
grid on;
xlabel('x');
ylabel('\phi_1');

time=toc;

save A2;

