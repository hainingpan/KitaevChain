t=1;
Delta=1;
delta=1e-3;
mulist=linspace(-3,3,400);
energylist=linspace(-2*t,2*t,400);
n=100;
ldosmap=zeros(length(mulist),length(energylist),n);
lenmu=length(mulist);
lenenergy=length(energylist);
parfor i=1:lenmu
    for j=1:lenenergy
    mu=mulist(i);
    energy=energylist(j);
    ldosmap(i,j,:)=ldosall(mu,t,Delta,n,energy,delta);
    end
end

figure;
surf(mulist,energylist,squeeze(ldosmap(:,:,1))','edgecolor','none');
view(2);
colorbar;colormap hot;
xlabel('\mu/t');
ylabel('E/t');
title('LDOS at the left end');

figure;
surf(mulist,energylist,squeeze(ldosmap(:,:,50))','edgecolor','none');
view(2);
colorbar;colormap hot;
xlabel('\mu/t');
ylabel('E/t');
title('LDOS at the bulk');

figure;
surf(mulist,energylist,squeeze(ldosmap(:,:,end))','edgecolor','none');
view(2);
colorbar;colormap hot;
xlabel('\mu/t');
ylabel('E/t');
title('LDOS at the right end');