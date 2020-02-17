%% ldos vs x
function re=ldosall(mu,t,Delta,n,omega,delta)
ham=H(mu,t,Delta,n);
G=inv(full((omega+1i*delta)*speye(2*n)-ham)); 
Gdiag=diag(G);
re=-imag(transpose(sum(reshape(Gdiag,[],2),2)))/pi;
end