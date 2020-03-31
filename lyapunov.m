function gamma=lyapunov(mulist,tlist,Delta)
%for an N-length wire
%mulist=[mu_1,mu_2,...,mu_N]
%tlist=[t_0,t_1,..,t_N]=[t_1,t_1,...,t_{N-1},t_{N-1}];
%Deltalist=[Delta_0,Delta_1,..,t_N]=[Delta_1,Delta_1,...,Delta_{N-1},Delta_{N-1}];
%Transfer matrix is A=prod(A_n), n=1..N
%A11=-mu_n/(t_n+Delta_n)
%A12=(Delta_{n-1}-t_{n-1})/(Delta_n+t_n)
%A21=1;A22=0;
A11=-mulist./(tlist(2:end)+Delta(2:end));
A12=-(tlist(1:end-1)-Delta(1:end-1))./(tlist(2:end)+Delta(2:end));
A=eye(2);
N=length(mulist);
for i=1:N
    An=[A11(i),A12(i);1,0];
    A=A*An;
end
lambda=eig(A);
gamma=(log(max(abs(lambda)))/length(mulist));
end
