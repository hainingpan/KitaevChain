function amp=wf(mu,Delta,muVar,mulist,n)
t=1;
if length(mulist)==1
    mulist=mu+muVar*randn(n,1);
end
ham=Hmu(mulist,t,Delta,n);
[vec,val]=eig(full(ham));
val=diag(val);
[val_s,val_I]=sort(val);
wf_p=vec(:,val_I(n+1));
wf_e=vec(:,val_I(n));

amp=abs(wf_p(1:n)).^2+abs(wf_p(n+1:2*n)).^2+abs(wf_e(1:n)).^2+abs(wf_e(n+1:2*n)).^2;
figure;
plot(1:n,amp);
title(strcat("Energy:",num2str(val_s(n+1))));
end
