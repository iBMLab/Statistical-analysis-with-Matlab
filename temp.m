nfactor=2;
nlevel=2;
ntrial=nlevel^nfactor*20;
DX=zeros(ntrial,nfactor);
for in=1:nfactor
    kk=nlevel^in;
    tkk=reshape(repmat(1:kk,[ntrial/kk,1]),[ntrial,1]);
    DX(:,in)=nlevel-mod(tkk,nlevel);
end
DX=nominal(DX);
tbl = array2table(DX);
%%
tbl.y=random('norm',0,2,[ntrial,1]);
mdl1 = fitlme(tbl,'y~DX1*DX2');anova(mdl1)
mdl2 = fitlme(tbl,'y~DX1*DX2','DummyvarCoding','effect');anova(mdl2)
mdl = fitlm(tbl,'y~DX1*DX2');anova(mdl)
%% 
nboot=1000;
Type2=zeros(nboot,1);
rng(1)
rmmean=1;
for iboot=1:nboot
    disp(iboot)
    y=random('norm',0,2,[ntrial,1]);
    if rmmeam==1
        
    end
    tbl.y=y;
    mdl = fitlm(tbl,'y~DX1*DX2*DX3');
    tbl2 = anova(mdl,'component',3);
    Type2(iboot,1)=sum(tbl2.pValue<.05);
end