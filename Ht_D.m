function Ham=Ht_D(mu,tlist,Delta,phi,n)
tauz=diag([1,-1]);
tauy=[0,-1i;1i,0];
band1sm=(spdiags([ones(n,1)],[1],n,n));
bandm1sm=(spdiags([ones(n,1)],[-1],n,n));
band_pos=band1sm+bandm1sm;
band_neg=band1sm-bandm1sm;

% Ham=-mu*kron(tauz,speye(n))-(kron(1i*Delta*exp(1i*phi)*tauy,band1sm)-kron(1i*Delta*exp(-1i*phi)*tauy,bandm1sm))-kron(tauz,tlistmat);
mublk=-mu*speye(n);
tblk=-spdiags([[0;tlist],[tlist;0]],[1,-1],n,n);
Deltablk=-Delta*exp(1i*phi)*band_neg;
Ham=[mublk+tblk,-Deltablk;conj(Deltablk),-mublk-tblk];
end