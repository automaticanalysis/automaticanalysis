a=Aacc();
a.aacc_crypto.setKey('rhodricusacklori');



figure(10);
tots=[];
for reps=1:1
    for m=1:2:4
        fprintf('Rep %d threads %d\n',reps,m);
        tic
        a.MAXTHREADS=m;
%        a.uploadDicomDirectory('/home/rhodri/Series_001_CBU_Localiser');
        a.uploadDicomDirectory('/home/rhodri/Series_002_CBU_MPRAGE');
        elapsed=toc;
        tots=[tots;m elapsed];
        scatter(tots(:,1),tots(:,2));
        drawnow;
    end;
end;