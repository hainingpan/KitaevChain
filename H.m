function Ham=H(mu,t,Delta,n)
tauz=diag([1,-1]);
tauy=[0,-1i;1i,0];
band1sm=(spdiags([ones(n,1)],[1],n,n));
bandm1sm=(spdiags([ones(n,1)],[-1],n,n));

Ham=-mu*kron(tauz,speye(n))-kron(t*tauz+1i*Delta*tauy,band1sm)-kron(t*tauz-1i*Delta*tauy,bandm1sm);
end