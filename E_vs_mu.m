t=1;
Delta=1;
mulist=-3:0.1:3;
n=100;
nv=10;
% nv=2*n;
en=zeros(nv,length(mulist));
for i=1:length(mulist)
    mu=mulist(i);
    ham=H(mu,t,Delta,n);
    eigval=eigs(ham,nv,1e-3,'Tolerance',1e-5,'MaxIterations',20000);
%     eigval=eig(ham);
    en(:,i)=sort(eigval(1:nv));
end
figure;
plot(mulist,en);