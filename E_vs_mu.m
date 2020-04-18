function [en,Pf,LE,tqpt_lb,tqpt_ub]=E_vs_mu(t,Delta,phi,tVar,tlist,n)
% t=1;
% Delta=0.2;
mulist=linspace(0,4*t,201);
% n=100;
% nv=10;
% nv=2*n;
en=zeros(2*n,length(mulist));
if length(tlist)==1
    tlist=t+tVar*randn(n-1,1);
%     tlist=t+tVar*(2*rand(n-1,1)-1);
end
Pf=zeros(1,length(mulist));
LE=zeros(1,length(mulist));
sy=[1,-1i;1,1i]/sqrt(2);
syt=kron(sy,eye(n));
% hamiltonian=@(mu) Ht(mu,tlist,Delta,n);
hamiltonian=@(mu) Ht_D(mu,tlist,Delta,phi,n);

for i=1:length(mulist)
    mu=mulist(i);
%     ham=Ht(mu,tlist,Delta,n);  
    ham=hamiltonian(mu);
%     ham=H(mu,t,Delta,n);
%     eigval=eigs(ham,nv,1e-40,'Tolerance',1e-5,'MaxIterations',20000);
%     en(:,i)=sort(eigval(1:nv));
    eigval=eig(full(ham));
    en(:,i)=sort(eigval);
    Pf(i)=sign(real(pfaffian_LTL(syt'*ham*syt)));    
    LE(i)=lyapunov(mu*ones(1,n),[tlist(1),tlist.',tlist(end)],Delta*ones(1,n+1));
end
figure;
plot(mulist/Delta,en'/Delta,'k');
xlim([0,mulist(end)/Delta]);
ylim([-2,2]);
xlabel('\mu/\Delta');
ylabel('E/\Delta');
xline(t/Delta*2,'r');
title(strcat('Energy [t/\Delta=',num2str(t/Delta),',\sigma_t/t=',num2str(tVar/t),']'));
fn_t=strcat('t',num2str(t));
fn_Delta=strcat('D',num2str(Delta),'phi',num2str(phi));
fn_tVar=strcat('tVar',num2str(tVar));
fn=strcat(fn_t,fn_Delta,fn_tVar);
saveas(gcf,strcat(fn,'.png'));

figure;
plot(mulist/Delta,Pf,'k');
xlim([0,mulist(end)/Delta]);
ylim([-1.2,1.2]);
xlabel('\mu/\Delta');
ylabel('Pf');
xline(t/Delta*2,'r');
title(strcat('Pfaffian at [t/\Delta=',num2str(t/Delta),',\sigma_t/t=',num2str(tVar/t),']'));
fn_t=strcat('t',num2str(t));
fn_Delta=strcat('D',num2str(Delta));
fn_tVar=strcat('tVar',num2str(tVar));
fn=strcat(fn_t,fn_Delta,fn_tVar);
saveas(gcf,strcat(fn,'_Pf.png'));

figure;
plot(mulist/Delta,LE,'k');
xlim([0,mulist(end)/Delta]);
xlabel('\mu/\Delta');
ylabel('LE');
yline(0);
xline(t/Delta*2,'r');
lb=min([find(LE>0,1),length(mulist)]);
ub=min([find(LE<0,1,'last'),length(mulist)]);
if lb<length(tlist)
    tqpt_lb=interp1(LE(lb-1:lb),tlist(lb-1:lb),0);
else
    tqpt_lb=[];
end
if ub<length(tlist)
    tqpt_ub=interp1(LE(ub:ub+1),tlist(ub:ub+1),0);
else
    tqpt_ub=[];
end
title(strcat('Lyapunov exponent [t/\Delta=',num2str(t/Delta),',\sigma_t/t=',num2str(tVar/t),'] TQPT\in[',num2str(tqpt_lb/Delta),',',num2str(tqpt_ub/Delta),']'));
fn_t=strcat('t',num2str(t));
fn_Delta=strcat('D',num2str(Delta));
fn_tVar=strcat('tVar',num2str(tVar));
fn=strcat(fn_t,fn_Delta,fn_tVar);
saveas(gcf,strcat(fn,'_LE.png'));
end