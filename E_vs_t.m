function E_vs_t(mu,Delta,muVar,mulist,n)
% mu=1;
% Delta=0.2;
tlist=linspace(0,2*mu,201);
% n=100;
% nv=10;
% nv=2*n;
en=zeros(2*n,length(tlist));
if length(mulist)==1
    mulist=mu+muVar*randn(n,1);E
end
parfor i=1:length(tlist)
    t=tlist(i);
    ham=Hmu(mulist,t,Delta,n);
%     ham=H(mu,t,Delta,n);
%     eigval=eigs(ham,nv,1e-40,'Tolerance',1e-5,'MaxIterations',20000);
%     en(:,i)=sort(eigval(1:nv));
    eigval=eig(ham);
    en(:,i)=sort(eigval);
end
figure;
plot(tlist/Delta,en'/Delta,'k');
xlim([0,tlist(end)/Delta]);
ylim([-(2*tlist(1)+mu)/Delta,(2*tlist(1)+mu)/Delta]);
xlabel('t/\Delta');
ylabel('E/\Delta');
line([mu/Delta/2,mu/Delta/2],[-(2*tlist(1)+mu)/Delta,(2*tlist(1)+mu)/Delta],'Color','r');
title(strcat('Energy at [\mu/\Delta=',num2str(mu/Delta),',\sigma_\mu/\mu=',num2str(muVar/mu),']'));
fn_mu=strcat('m',num2str(mu));
fn_Delta=strcat('D',num2str(Delta));
fn_muVar=strcat('muVar',num2str(muVar));
fn=strcat(fn_mu,fn_Delta,fn_muVar);
saveas(gcf,strcat(fn,'.png'));
save(strcat(fn,'_',num2str(tlist(end)),'.dat'),'en','-ascii');
end