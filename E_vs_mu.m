function E_vs_mu(t,Delta,tVar,tlist,n)
% t=1;
% Delta=0.2;
mulist=linspace(0,4*t,201);
% n=100;
% nv=10;
% nv=2*n;
en=zeros(2*n,length(mulist));
if length(tlist)==1
    tlist=t+tVar*randn(n-1,1);
end
parfor i=1:length(mulist)
    mu=mulist(i);
    ham=Ht(mu,tlist,Delta,n);
%     ham=H(mu,t,Delta,n);
%     eigval=eigs(ham,nv,1e-40,'Tolerance',1e-5,'MaxIterations',20000);
%     en(:,i)=sort(eigval(1:nv));
    eigval=eig(ham);
    en(:,i)=sort(eigval);
end
figure;
plot(mulist/Delta,en'/Delta,'k');
xlim([0,mulist(end)/Delta]);
ylim([-2,2]);
xlabel('\mu/\Delta');
ylabel('E/\Delta');
line([t/Delta*2,t/Delta*2],[-2,2],'Color','r');
title(strcat('Energy [t/\Delta=',num2str(t/Delta),',\sigma_t/t=',num2str(tVar/t),']'));
fn_t=strcat('t',num2str(t));
fn_Delta=strcat('D',num2str(Delta));
fn_tVar=strcat('tVar',num2str(tVar));
fn=strcat(fn_t,fn_Delta,fn_tVar);
saveas(gcf,strcat(fn,'.png'));
save(strcat(fn,'_',num2str(mulist(end)),'.dat'),'en','-ascii');
end