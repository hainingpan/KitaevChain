t=1;
Delta=1;
mulist=linspace(-3,3,60);
n=100;
nv=10;
% nv=2*n;
en=zeros(nv,length(mulist));
for i=1:length(mulist)
    mu=mulist(i);
    ham=H(mu,t,Delta,n);
%     eigval=eigs(ham,nv,1e-40,'Tolerance',1e-5,'MaxIterations',20000);
%     en(:,i)=sort(eigval(1:nv));

    eigval=eig(ham);
    en(:,i)=sort(eigval);
end
figure;
plot(mulist,en);