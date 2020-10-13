function [en,Pf,LE,tqpt_lb,tqpt_ub]=E_vs_t(mu,Delta,phi,muVar,mulist,n)
% mu=1;
% Delta=0.2;
if muVar<2
    tlist=linspace(0,mu,401);
else
    tlist=linspace(0,8*mu,8*401);
end
% n=100;
% nv=10;
% nv=2*n;
en=zeros(2*n,length(tlist));
if length(mulist)==1
    mulist=mu+muVar*randn(n,1);
%     mulist=mu+muVar*(2*rand(n,1)-1);
end
Pf=zeros(1,length(tlist));
LE=zeros(1,length(tlist));
sy=[1,-1i;1,1i]/sqrt(2);
syt=kron(sy,eye(n));
% hamiltonian=@(t) Hmu(mulist,t,Delta,n);
hamiltonian=@(t) Hmu_D(mulist,t,Delta,phi,n);

for i=1:length(tlist)
    t=tlist(i);
    ham=hamiltonian(t);
%     ham=H(mu,t,Delta,n);
%     eigval=eigs(ham,nv,1e-40,'Tolerance',1e-5,'MaxIterations',20000);
%     en(:,i)=sort(eigval(1:nv));
    eigval=eig(full(ham));
    en(:,i)=sort(eigval);
    Pf(i)=sign(real(pfaffian_LTL(syt'*ham*syt)));    
    LE(i)=lyapunov(mulist,t*ones(1,n+1),Delta*ones(1,n+1));
end
figure;
plot(tlist/Delta,en'/Delta,'k');
xlim([0,tlist(end)/Delta]);
ylim([-(2*tlist(1)+mu)/Delta,(2*tlist(1)+mu)/Delta]);
xlabel('t/\Delta');
ylabel('E/\Delta');
xline(mu/Delta/2,'r');
title(strcat('Energy at [\mu/\Delta=',num2str(mu/Delta),',\sigma_\mu/\mu=',num2str(muVar/mu),']'));
fn_mu=strcat('m',num2str(mu));
fn_Delta=strcat('D',num2str(Delta),'phi',num2str(phi));
fn_muVar=strcat('muVar',num2str(muVar));
fn_L=strcat('L',num2str(n));
fn=strcat(fn_mu,fn_Delta,fn_muVar,fn_L);
saveas(gcf,strcat(fn,'.png'));

figure;
plot(tlist/Delta,Pf);
xlim([0,tlist(end)/Delta]);
ylim([-1.2,1.2]);
xlabel('t/\Delta');
ylabel('Pf');
xline(mu/Delta/2,'r');
title(strcat('Pfaffian at [\mu/\Delta=',num2str(mu/Delta),',\sigma_\mu/\mu=',num2str(muVar/mu),']'));
fn_mu=strcat('m',num2str(mu));
fn_Delta=strcat('D',num2str(Delta),'phi',num2str(phi));
fn_muVar=strcat('muVar',num2str(muVar));
fn=strcat(fn_mu,fn_Delta,fn_muVar);
saveas(gcf,strcat(fn,'_Pf.png'));

figure;
plot(tlist/Delta,LE);
xlim([0,tlist(end)/Delta]);
xlabel('t/\Delta');
ylabel('LE');
yline(0);
xline(mu/Delta/2,'r');
lb=min([find(LE<0,1),length(tlist)]);
ub=min([find(LE>0,1,'last'),length(tlist)]);
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
title(strcat('Lyapunov exponent at [\mu/\Delta=',num2str(mu/Delta),',\sigma_\mu/\mu=',num2str(muVar/mu),'] TQPT\in[',num2str(tqpt_lb/Delta),',',num2str(tqpt_ub/Delta),']'));
fn_mu=strcat('m',num2str(mu));
fn_Delta=strcat('D',num2str(Delta),'phi',num2str(phi));
fn_muVar=strcat('muVar',num2str(muVar));
fn=strcat(fn_mu,fn_Delta,fn_muVar);
saveas(gcf,strcat(fn,'_LE.png'));
end