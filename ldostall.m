function re=ldostall(mu,tlist,Delta,n,omega,delta)
% ldos vs x at disorder in t
ham=Ht(mu,tlist,Delta,n);
G=inv(full((omega+1i*delta)*speye(2*n)-ham)); 
Gdiag=diag(G);
re=-imag(transpose(sum(reshape(Gdiag,[],2),2)))/pi;
end