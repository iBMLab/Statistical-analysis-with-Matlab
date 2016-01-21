% Week 3
% This week is the starting point of linear statistics: regression, t-test
% and ANOVA
%% regression demo
rand('state',0);
rho=[.9,.6,.3,0];
slope=[1,.4,-.2,-.8];
irr=1;
iss=1;
n = 200;% N trials
iplot=0;
mu1=[0 0];
Sigma1=[1 0; 0 1];
R1=chol(Sigma1);
Z1 = repmat(mu1,n,1) + randn(n,2)*R1;
X=Z1(:,1);
% Gaussian copula
figure;hold on
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
% normalized
constant=min(U(:,2));
cst=constant*sign(constant);
U(:,2)=U(:,2)+cst;
ytmp=U(:,2);
scatterhist(X,ytmp);hold on
plot(X,X.*slope(iss)+cst,'r')
title(['Y_' num2str(iplot) '=' num2str(slope(iss)) '*X+c_' num2str(iplot) '; (rho_' num2str(iplot) '=' num2str(rho(irr)) ')']);
% axis([-2.2 2.2 0 5.5])
xlabel('Rating');
ylabel('Response');

[rho1,pval]=corr(X,ytmp);
covm=cov(X,ytmp);
rho2=covm(2,1)./(sqrt(covm(1,1))*sqrt(covm(2,2)));
[beta,bci,rho2,RINT,STATS] = regress(ytmp,[X,ones(length(X),1)]);
%%
figure;

for irr=1:length(rho)
    for iss=1:length(slope)
        % Gaussian copula
        subplot(4,4,iplot);hold on
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
        % normalized
        constant=min(U(:,2));
        cst=constant*sign(constant);
        U(:,2)=U(:,2)+cst;
        
        plot(X,U(:,2),'.k');
        plot(X,X.*slope(iss)+cst,'r')
        title(['Y_' num2str(iplot) '=' num2str(slope(iss)) '*X+c_' num2str(iplot) '; (rho_' num2str(iplot) '=' num2str(rho(irr)) ')']);
        axis([-2.2 2.2 0 5.5])
    end
end
xlabel('Rating');
ylabel('Response');
%% Understand p-value: fisher's exact test
% example in http://mathworld.wolfram.com/FishersExactTest.html
x = table([5;1],[0;4],'VariableNames',{'Math_Mag','Science'},'RowNames',{'Math','Biology'})
[h,p,stats] = fishertest(x)
% example in https://en.wikipedia.org/wiki/Fisher%27s_exact_test
% x = table([1;11],[9;3],'VariableNames',{'Men','Women'},'RowNames',{'Dieting','Non_dieting'})
% [h,p,stats] = fishertest(x,'Tail','left','Alpha',0.05)
%% One-sample and paired-sample t-test
% doc ttest
rand('state',0);

nitem=20;
m1=1;
v1=1;
x = random('normal',m1,v1,[nitem,1]);
[H1,P1,CI1,STATS1] = ttest(x,m1)
[H2,P2,CI2,STATS2] = ttest(x-m1)

m2=1.2;
y = random('normal',m2,v1,[nitem,1]);
[H3,P3,CI3,STATS3] = ttest(x,y)
[H4,P4,CI4,STATS4] = ttest(x-y)
%% Two-sample t-test (or independent sample t-test)
% doc ttest2
rand('state',0);

m1=1;
v1=1;
x = random('normal',m1,v1,[nitem,1]);

m2=1.2;
v2=1.3;
y = random('normal',m2,v2,[nitem,1]);
[H5,P5,CI5,STATS5] = ttest(x,y)
[H6,P6,CI6,STATS6] = ttest2(x,y)
%% anova
% http://ch.mathworks.com/help/stats/analysis-of-variance-and-covariance.html
% one-way ANOVA
% for two-sample case, it's the same as the two-sample t-test
Y=[x,y]
[p,tbl,stats] = anova1(Y)

% but now you can extent it to multiple case (e.g., three sample case):
Y=zeros(nitem,3);
m1=1;v1=1;
m2=2;v2=1.2;
m3=1.6;v3=.9;
Y(:,1) = random('normal',m1,v1,[nitem,1]);
Y(:,2) = random('normal',m2,v2,[nitem,1]);
Y(:,3) = random('normal',m3,v3,[nitem,1]);
[p,tbl,stats] = anova1(Y)


