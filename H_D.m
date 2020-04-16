function Ham=H_D(mu,t,Delta,phi,n)
%Kitaev chain in class D, where Delta is complex phi
tauz=diag([1,-1]);
tauy=[0,-1i;1i,0];
band1sm=(spdiags([ones(n,1)],[1],n,n));
bandm1sm=(spdiags([ones(n,1)],[-1],n,n));
band_pos=band1sm+bandm1sm;
band_neg=band1sm-bandm1sm;

mublk=-mu*speye(n);
tblk=-t*band_pos;
Deltablk=-Delta*exp(1i*phi)*band_neg;
Ham=[mublk+tblk,-Deltablk;conj(Deltablk),-mublk-tblk];
end