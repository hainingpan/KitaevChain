function re=dos(mu,t,Delta,n,omega,delta)
ham=H(mu,t,Delta,n);
K=(omega+1i*delta)*speye(2*n)-ham;
Kv=eigs(K,40,0,'Tolerance',1e-5,'MaxIterations',20000);
re=-sum(imag(1./Kv))/pi;
end