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
[H7,P7,CI7,STATS7] = ttest2(x,y,'Vartype','unequal')
%% more t-test (in a linear model)
rand('state',0);

m1=1;
v1=1;
nitem1=20;
x = random('normal',m1,v1,[nitem1,1]);

m2=1.2;
v2=1.3;
nitem2=25;
y = random('normal',m2,v2,[nitem2,1]);
[uH1,uP1,uCI1,uSTATS1] = ttest2(x,y)
% [uH2,uP2,uCI2,uSTATS2] = ttest2(x,y,'Vartype','unequal')

X=[ones(size(x));ones(size(y))*2];
X(:,2)=[ones(size([x;y]))];
[beta,bci,rho2,RINT,STATS] = regress([x;y],X);
%% Type 1 error and multiple comparison correction
m=1;
v=1;
nitem2=100;
ntest=10000;
Xall = random('normal',m,v,[nitem2,ntest]);
Yall = random('normal',m,v,[nitem2,ntest]);
[uH3,uP3,uCI3,uSTATS3] = ttest2(Xall,Yall,'Dim',1,'Vartype','unequal');
Type1error=sum(uH3)./ntest;
% bonferroni correction
pnew=.05./ntest;
Type1error2=sum(uP3<=pnew)./ntest;

%% anova
% http://ch.mathworks.com/help/stats/analysis-of-variance-and-covariance.html
% one-way ANOVA
% for two-sample case, it's the same as the two-sample t-test

m1=1;
v1=1;
x = random('normal',m1,v1,[nitem,1]);

m2=1.2;
v2=1.3;
y = random('normal',m2,v2,[nitem,1]);
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


