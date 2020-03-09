function Ham=Ht(mu,tlist,Delta,n)
tauz=diag([1,-1]);
tauy=[0,-1i;1i,0];
band1sm=(spdiags([ones(n,1)],[1],n,n));
bandm1sm=(spdiags([ones(n,1)],[-1],n,n));

tlistmat=spdiags([[0;tlist],[tlist;0]],[1,-1],n,n);
Ham=-mu*kron(tauz,speye(n))-kron(1i*Delta*tauy,band1sm-bandm1sm)-kron(tauz,tlistmat);
end