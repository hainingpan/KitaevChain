function LDOS_vs_t(mu,Delta,muVar,mulist,n)
% mu=1;
% Delta=.2;
delta=1e-3;
tlist=linspace(0,2*mu,201);
energylist=linspace(-mu,mu,201);
% n=100;
ldosmap=zeros(length(tlist),length(energylist),n);
lent=length(tlist);
lenenergy=length(energylist);
% muVar=1;
if length(mulist)==1
    mulist=mu+muVar*randn(n,1);
end

parfor i=1:lent
    for j=1:lenenergy
    t=tlist(i);
    energy=energylist(j);
    ldosmap(i,j,:)=ldosmuall(mulist,t,Delta,n,energy,delta);
    end
end

E_vs_t(mu,Delta,muVar,mulist,n);
figure;
LDOS_L=(squeeze(ldosmap(:,:,1)))';
surf(tlist/Delta,energylist/Delta,LDOS_L,'edgecolor','none');
view(2);
colorbar;colormap hot;
xlabel('t/\Delta');
ylabel('E/\Delta');
title(strcat('LDOS on the left end [\mu/\Delta=',num2str(mu/Delta),',\sigma_\mu/\mu=',num2str(muVar/mu),']'));
fn_mu=strcat('m',num2str(mu));
fn_Delta=strcat('D',num2str(Delta));
fn_muVar=strcat('muVar',num2str(muVar));
fn=strcat(fn_mu,fn_Delta,fn_muVar,'_LDOS_L');
saveas(gcf,strcat(fn,'.png'));
save(strcat(fn,'.dat'),'LDOS_L','-ascii');

figure;
LDOS_M=(squeeze(ldosmap(:,:,floor(n/2))))';
surf(tlist/Delta,energylist/Delta,LDOS_M,'edgecolor','none');
view(2);
colorbar;colormap hot;
xlabel('t/\Delta');
ylabel('E/\Delta');
title(strcat('LDOS in the bulk [\mu/\Delta=',num2str(mu/Delta),',\sigma_\mu/\mu=',num2str(muVar/mu),']'));
fn=strcat(fn_mu,fn_Delta,fn_muVar,'_LDOS_M');
saveas(gcf,strcat(fn,'.png'));
save(strcat(fn,'.dat'),'LDOS_M','-ascii');


figure;
LDOS_R=(squeeze(ldosmap(:,:,end)))';
surf(tlist/Delta,energylist/Delta,LDOS_R,'edgecolor','none');
view(2);
colorbar;colormap hot;
xlabel('t/\Delta');
ylabel('E/\Delta');
title(strcat('LDOS on the right end [\mu/\Delta=',num2str(mu/Delta),',\sigma_\mu/\mu=',num2str(muVar/mu),']'));
fn=strcat(fn_mu,fn_Delta,fn_muVar,'_LDOS_R');
saveas(gcf,strcat(fn,'.png'));
save(strcat(fn,'.dat'),'LDOS_R','-ascii');


figure;
DOS=(mean(ldosmap,3))';
surf(tlist/Delta,energylist/Delta,DOS,'edgecolor','none');
view(2);
colorbar;colormap hot;
xlabel('t/\Delta');
ylabel('E/\Delta');
title(strcat('DOS [\mu/\Delta=',num2str(mu/Delta),',\sigma_\mu/\mu=',num2str(muVar/mu),']'));
fn=strcat(fn_mu,fn_Delta,fn_muVar,'_DOS');
saveas(gcf,strcat(fn,'.png'));
save(strcat(fn,'.dat'),'DOS','-ascii');

fn=strcat(fn_mu,fn_Delta,fn_muVar,'_mulist');
save(strcat(fn,'.dat'),'mulist','-ascii');

