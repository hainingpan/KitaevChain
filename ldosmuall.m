function re=ldosmuall(mulist,t,Delta,phi,n,omega,delta)
% ldos vs x at disorder in mu
% ham=Hmu(mulist,t,Delta,n);
ham=Hmu_D(mulist,t,Delta,phi,n);
G=inv(full((omega+1i*delta)*speye(2*n)-ham)); 
Gdiag=diag(G);
re=-imag(transpose(sum(reshape(Gdiag,[],2),2)))/pi;
end