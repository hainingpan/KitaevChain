function Ham=Hmu_D(mulist,t,Delta,phi,n)
tauz=diag([1,-1]);
tauy=[0,-1i;1i,0];
band1sm=(spdiags([ones(n,1)],[1],n,n));
bandm1sm=(spdiags([ones(n,1)],[-1],n,n));
band_pos=band1sm+bandm1sm;
band_neg=band1sm-bandm1sm;

% Ham=-kron(tauz,mudiag)-kron(t*tauz+1i*Delta*exp(1i*phi)*tauy,band1sm)-kron(t*tauz-1i*Delta*exp(-1i*phi)*tauy,bandm1sm);

mublk=-spdiags(mulist,[0],n,n);
tblk=-t*band_pos;
Deltablk=-Delta*exp(1i*phi)*band_neg;
Ham=[mublk+tblk,-Deltablk;conj(Deltablk),-mublk-tblk];
end