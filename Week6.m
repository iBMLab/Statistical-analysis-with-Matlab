%% GLM - regession model
rand('state',0);
rho=[.9,.6,.3,0];
slope=[1,.4,-.2,-.8];
irr=2;
iss=1;
n = 200;% N trials
iplot=0;
mu1=[0 0];
Sigma1=[1 0; 0 1];
R1=chol(Sigma1);
Z1 = repmat(mu1,n,1) + randn(n,2)*R1;
X=Z1(:,1);
% Gaussian copula
% mu=[0 0];
Sigma=[1 rho(irr); rho(irr) 1];
% Z = mvnrnd(mu,Sigma,n);
% U = normcdf(Z,0,1);
R=chol(Sigma);
Z1(:,2)=randn(n,1);
U = Z1*R;

iplot=iplot+1;

% introduce slope
U(:,2)=U(:,2).*slope(iss);
% introduce constant
constant=min(U(:,2));
cst=constant*sign(constant);
U(:,2)=U(:,2)+cst;
%
ytmp=U(:,2);
%
figure;
scatterhist(X,ytmp,'kernel','on','Location','SouthEast',...
    'Direction','out','LineStyle',{'-'},...
    'LineWidth',1.5,'MarkerSize',10,'Marker','.');
hold on
[rho1,pval]=corr(X,ytmp);
covm=cov(X,ytmp);
rho2=covm(2,1)./(sqrt(covm(1,1))*sqrt(covm(2,2)));
[beta,bci,rho2,RINT,STATS] = regress(ytmp,[X,ones(length(X),1)]);

plot(X,X.*beta(1)+beta(2),'r','LineWidth',2) 
% title(['Y_' num2str(iplot) '=' num2str(slope(iss)) '*X+c_' num2str(iplot) '; (rho_' num2str(iplot) '=' num2str(rho(irr)) ')']);
% axis([-2.2 2.2 0 5.5])
xlabel('Predictor');
ylabel('Response');

% introduce subject information
Ns=5;
sbjidx=reshape(repmat(1:Ns,n/Ns,1),n,1);
figure;
scatterhist(X,ytmp,'Group',sbjidx,'kernel','on','Location','SouthEast',...
    'Direction','out','LineStyle',{'-'},...
    'LineWidth',1.5,'MarkerSize',10,'Marker','.');

hold on

plot(X,X.*beta(1)+beta(2),'r','LineWidth',2)
% title(['Y_' num2str(iplot) '=' num2str(slope(iss)) '*X+c_' num2str(iplot) '; (rho_' num2str(iplot) '=' num2str(rho(irr)) ')']);
% axis([-2.2 2.2 0 5.5])
for is=1:Ns
    X2=X(sbjidx==is);
    Y2=ytmp(sbjidx==is);
    [betat,bci,rho2,RINT,STATS] = regress(Y2,[X2,ones(length(X2),1)]);
    plot(X,X.*betat(1)+betat(2),'--')
end
xlabel('Predictor');
ylabel('Response');
%% introduce intercept for subject
figure;
SbjItp=random('normal',0,2,Ns,1);
ytmpnew=ytmp;
for is=1:Ns
    ytmpnew(sbjidx==is)=ytmp(sbjidx==is)+SbjItp(is);
end
scatterhist(X,ytmpnew,'Group',sbjidx,'kernel','on','Location','SouthEast',...
    'Direction','out','LineStyle',{'-'},...
    'LineWidth',1.5,'MarkerSize',10,'Marker','.');
hold on
[beta,bci,rho2,RINT,STATS] = regress(ytmpnew,[X,ones(length(X),1)]);

plot(X,X.*beta(1)+beta(2),'r','LineWidth',2) 
% title(['Y_' num2str(iplot) '=' num2str(slope(iss)) '*X+c_' num2str(iplot) '; (rho_' num2str(iplot) '=' num2str(rho(irr)) ')']);
% axis([-2.2 2.2 0 5.5])
for is=1:Ns
    X2=X(sbjidx==is);
    Y2=ytmpnew(sbjidx==is);
    [betat,bci,rho2,RINT,STATS] = regress(Y2,[X2,ones(length(X2),1)]);
    plot(X,X.*betat(1)+betat(2),'--')
end
xlabel('Predictor');
ylabel('Response');
%% introduced subject related intercept & slope
Sbjslp=random('normal',slope(iss),.9,1,Ns);
Yorg=U(:,2);
Yorg=Yorg.*reshape(repmat(Sbjslp,n/Ns,1),n,1);
% introduce constant
% constant=min(Yorg);
% cst=constant*sign(constant);
% Yorg=Yorg+cst;
%
ytmp2t=Yorg;
figure;
% SbjItp=random('normal',0,2,Ns,1);
ytmp2=ytmp2t;
for is=1:Ns
    ytmp2(sbjidx==is)=ytmp2t(sbjidx==is)+SbjItp(is);
end
scatterhist(X,ytmp2,'Group',sbjidx,'kernel','on','Location','SouthEast',...
    'Direction','out','LineStyle',{'-'},...
    'LineWidth',1.5,'MarkerSize',10,'Marker','.');
hold on
[beta,bci,rho2,RINT,STATS] = regress(ytmp2,[X,ones(length(X),1)]);

plot(X,X.*beta(1)+beta(2),'r','LineWidth',2) 
% title(['Y_' num2str(iplot) '=' num2str(slope(iss)) '*X+c_' num2str(iplot) '; (rho_' num2str(iplot) '=' num2str(rho(irr)) ')']);
% axis([-2.2 2.2 0 5.5])
for is=1:Ns
    X2=X(sbjidx==is);
    Y2=ytmp2(sbjidx==is);
    [betat,bci,rho2,RINT,STATS] = regress(Y2,[X2,ones(length(X2),1)]);
    plot(X,X.*betat(1)+betat(2),'--')
end

xlabel('Predictor');
ylabel('Response');

%% fixed and random slope model (for categorical predictors)
C1m=2;
C2m=5;
figure;hold on
for is=1:10
    matmp=zeros(100,2);
    matmp(1:50,1)=random('normal',10,.1,50,1);
    matmp(51:100,1)=random('normal',20,.1,50,1);
    matmp(1:50,2)=random('normal',C1m,1,50,1);
    matmp(51:100,2)=random('normal',C2m+randn,1+randn,50,1);
    plot(matmp(:,1),matmp(:,2),'.','MarkerSize',10);
    plot([10,20],[mean(matmp(1:50,2),1),mean(matmp(51:100,2),1)],'lineWidth',2)
end
xlim([5,25])
%% HLM LMM rmANOVA comparison (Mixed model)
Ns=20;
Condi1=2;
Condi2=3;
Nt=10;
beta0=[1 2 3 1 2 3]';% [g1_c1 g1_c2 g1_c3, g2_c1 g2_c2 g2_c3]
%
DXss=kron(eye(6),ones(Nt,1));
DXall=repmat(DXss,[Ns,1]);
% sbjidx=reshape(repmat(1:Ns,[Nt*Condi1*Condi2,1]),[Ns*Nt*Condi1*Condi2,1]);
sbjidx=reshape(repmat(1:(Ns*2),[Nt*Condi2,1]),[Ns*Nt*Condi1*Condi2,1]);
sbjintercept=random('Normal',0,1.5,Ns,1);
sbjintM=zeros(size(sbjidx));
for is=1:Ns
    sbjintM(sbjidx==is)=sbjintercept(is);
end
Y=random('Normal',0,5,Ns*Condi1*Condi2*Nt,1)+DXall*beta0+sbjintM;

%
% HLM
beta1=zeros(Ns,Condi2);
DXss2=DXss(1:Nt*Condi2,1:Condi2);
for is=1:Ns
    Ytmp=Y(sbjidx==is);
    DXtmp=[DXss2,ones(length(DXss2),1)];
    % DXtmp=DXss2;
    invD=pinv(DXtmp);
    betatmp=invD*Ytmp;
    beta1(is,:)=betatmp(1:Condi2);
end
DX2tmp=repmat(1:Condi2,Ns,1);DX2tmp(2:2:Ns,:)=DX2tmp(2:2:Ns,:)+Condi2;
Y2=beta1(:);DX2=DX2tmp(:);
DXg=dummyvar(DX2);
beta2=DXg\Y2;
mse=mean((DXg*beta2-Y2).^2);
covhlm=mse*((DXg'*DXg)^-1);

%
% LMM
lmm=fitlmematrix(DXall,Y,ones(length(Y),1),sbjidx,'Dummyvarcoding','effect');%,'CovariancePattern','Diagonal');
betalmm=lmm.Coefficients.Estimate;
covlmm=lmm.CoefficientCovariance;

% LMM 2
tbl=dataset;
tbl.Y=Y;
tbl.Condi1=nominal(2-double(sum(DXall(:,1:3),2)>0));
tbl.Condi2=nominal(double(sum(DXall(:,[1 4]),2)>0)+double(sum(DXall(:,[2 5]),2)>0)*2+double(sum(DXall(:,[3 6]),2)>0)*3);
tbl.sbjidx=sbjidx;
lmm2=fitlme(tbl,'Y~Condi1*Condi2+(1|sbjidx)','Dummyvarcoding','effect');
anova(lmm2)

barall=zeros(2,3);
for ic1=1:2
    for ic2=1:3
        indx=double(tbl.Condi1)==ic1&double(tbl.Condi2)==ic2;
        barall(ic1,ic2)=mean(Y(indx));
    end
end
figure;
subplot(1,3,1)
bar(mean(barall,2))

subplot(1,3,2)
bar(mean(barall,1))

subplot(1,3,3)
bar(barall)

% ANOVA
dfehlm=length(Y2)-rank(DXg);
C=limo_OrthogContrasts([Condi1,Condi2]);
clear pl Fl df1l df2l ph Fh df1h df2h
for ic=1:3
c=C{ic};
[pl(ic),Fl(ic),df1l(ic),df2l(ic)]=coefTest(lmm,c);
df1h(ic)=rank(c);
df2h(ic)=dfehlm;%(Ns-1)*rank(c);
Fh(ic)=((c*beta2)'*((c*covhlm*c')^-1)*(c*beta2))./rank(c);
ph(ic)=1-fcdf(Fh(ic),df1h(ic),df2h(ic));
end
FtableLME=dataset(Fl',df1l',df2l',pl','Varnames',{'Fstat','DF1','DF2','pValue'},'ObsNames',{'Condi1','Condi2','Interation'})
FtableHLM=dataset(Fh',df1h',df2h',ph','Varnames',{'Fstat','DF1','DF2','pValue'},'ObsNames',{'Condi1','Condi2','Interation'})

% rmANOVA
rmdatamat=zeros(Ns,Nt*Condi2);
ii=0;
for is=1:Ns
    ii=ii+1;
    indx=find(tbl.sbjidx==is);
    rmdatamat(ii,:)=tbl.Y(indx);
end
within=[tbl(indx,[3])];
varname=cell(size(within,1),1);
for iv=1:size(within,1)
    varname{iv,:}=strjoin({'y',num2str(iv)},'');
end
varname{iv+1,:}='Condi1';
rmdatamat(2:2:Ns,end+1)=1;
rmdatamat(:,end)=rmdatamat(:,end)+1;
between=mat2dataset(rmdatamat,'Varnames',varname);
between.Condi1=ordinal(between.Condi1);
rm = fitrm(between,string([ 'y1-y' num2str(iv) ' ~ Condi1'] ),'WithinDesign',within,'WithinModel','Condi2');
anovatbl = anova(rm);anovatbl(2,:)
ranovatbl = ranova(rm)
