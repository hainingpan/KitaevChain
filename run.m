% LDOS_vs_t(1,0.2,0,0.2,0,20)

dirlist=1:100;
muVarlist=[.4,.6,.8,1,1.5,2];
addpath(pwd)

for dir=dirlist
    mkdir(num2str(dir))
    cd(num2str(dir))
    for muVar=muVarlist
        fprintf("%i %0.2f\n",dir,muVar)
        LDOS_vs_t(1,0.2,0,muVar,0,20)
    end
    cd ..
end
        