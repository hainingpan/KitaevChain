function LDOS_vs_mu(t,Delta,tVar,tlist,n)
% t=1;
% Delta=0.2;
delta=1e-3;
mulist=linspace(0,4*t,201);
energylist=linspace(-2*Delta,2*Delta,201);
% n=100;
ldosmap=zeros(length(mulist),length(energylist),n);
lenmu=length(mulist);
lenenergy=length(energylist);

% tVar=1;
if length(tVar)==1
    tlist=t+tVar*randn(n-1,1);
    save('tlist.dat','tlist','-ascii');
end
parfor i=1:lenmu
    for j=1:lenenergy
    mu=mulist(i);
    energy=energylist(j);
    ldosmap(i,j,:)=ldostall(mu,tlist,Delta,n,energy,delta);
    end
end


E_vs_mu(t,Delta,tVar,tlist,n);

figure;
LDOS_L=(squeeze(ldosmap(:,:,1)))';
surf(mulist/Delta,energylist/Delta,LDOS_L,'edgecolor','none');
view(2);
colorbar;colormap hot;
xlabel('\mu/\Delta');
ylabel('E/\Delta');
xlim([0,mulist(end)/Delta]);
ylim([-2,2]);
title(strcat('LDOS on the left end [t/\Delta=',num2str(t/Delta),',\sigma_t/t=',num2str(tVar/t),']'));
fn_t=strcat('t',num2str(t));
fn_Delta=strcat('D',num2str(Delta));
fn_tVar=strcat('muVar',num2str(tVar));
fn=strcat(fn_t,fn_Delta,fn_tVar,'_LDOS_L');
saveas(gcf,strcat(fn,'.png'));
save(strcat(fn,'.dat'),'LDOS_L','-ascii');


figure;
LDOS_M=(squeeze(ldosmap(:,:,floor(n/2))))';
surf(mulist/Delta,energylist/Delta,LDOS_M,'edgecolor','none');
view(2);
colorbar;colormap hot;
xlabel('\mu/\Delta');
ylabel('E/\Delta');
xlim([0,mulist(end)/Delta]);
ylim([-2,2]);
title(strcat('LDOS in the bulk [t/\Delta=',num2str(t/Delta),',\sigma_t/t=',num2str(tVar/t),']'));
fn_t=strcat('t',num2str(t));
fn_Delta=strcat('D',num2str(Delta));
fn_tVar=strcat('muVar',num2str(tVar));
fn=strcat(fn_t,fn_Delta,fn_tVar,'_LDOS_M');
saveas(gcf,strcat(fn,'.png'));
save(strcat(fn,'.dat'),'LDOS_M','-ascii');


figure;
LDOS_R=(squeeze(ldosmap(:,:,end)))';
surf(mulist/Delta,energylist/Delta,LDOS_R,'edgecolor','none');
view(2);
colorbar;colormap hot;
xlabel('\mu/\Delta');
ylabel('E/\Delta');
xlim([0,mulist(end)/Delta]);
ylim([-2,2]);
title(strcat('LDOS on the right end [t/\Delta=',num2str(t/Delta),',\sigma_t/t=',num2str(tVar/t),']'));
fn_t=strcat('t',num2str(t));
fn_Delta=strcat('D',num2str(Delta));
fn_tVar=strcat('muVar',num2str(tVar));
fn=strcat(fn_t,fn_Delta,fn_tVar,'_LDOS_R');
saveas(gcf,strcat(fn,'.png'));
save(strcat(fn,'.dat'),'LDOS_R','-ascii');


figure;
DOS=(mean(ldosmap,3))';
surf(mulist/Delta,energylist/Delta,DOS,'edgecolor','none');
view(2);
colorbar;colormap hot;
xlabel('\mu/\Delta');
ylabel('E/\Delta');
xlim([0,mulist(end)/Delta]);
ylim([-2,2]);
title(strcat('DOS [t/\Delta=',num2str(t/Delta),',\sigma_t/t=',num2str(tVar/t),']'));
fn_t=strcat('t',num2str(t));
fn_Delta=strcat('D',num2str(Delta));
fn_tVar=strcat('muVar',num2str(tVar));
fn=strcat(fn_t,fn_Delta,fn_tVar,'_DOS');
saveas(gcf,strcat(fn,'.png'));
save(strcat(fn,'.dat'),'DOS','-ascii');
