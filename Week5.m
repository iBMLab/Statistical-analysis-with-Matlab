%% anova (1 way ANOVA)
% http://ch.mathworks.com/help/stats/analysis-of-variance-and-covariance.html
% one-way ANOVA
% for two-sample case, it's the same as the two-sample t-test
% it's like a t-test, but for multiple condition
nitem=20;
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
m3=1.6;v3=.9;
Y(:,1) = random('normal',m1,v1,[nitem,1]);
Y(:,2) = random('normal',m2,v2,[nitem,1]);
Y(:,3) = random('normal',m3,v3,[nitem,1]);
[p,tbl,stats] = anova1(Y)

%% one way anova in a regression model
X=zeros(size(Y,1)*size(Y,2),size(Y,2));
X(1:size(Y,1),1)=1;
X([size(Y,1)+1:size(Y,1)*2],2)=1;
X([size(Y,1)*2+1:size(Y,1)*3],3)=1;

Y2=[Y(:,1);Y(:,2);Y(:,3)];
[beta,bci,rho2,RINT,STATS] = regress(Y2,X);

disp(['The F value from the one way ANOVA is ' num2str(tbl{2,5})])
disp(['The p value from the one way ANOVA is ' num2str(tbl{2,6})])
disp(' ')
disp(['The F value from the linear model is ' num2str(STATS(2))])
disp(['The p value from the linear model is ' num2str(STATS(3))])

%% two-way anova with balanced designs
Cond2=2;
ntime=20;
Y2way=zeros(ntime*Cond2,3);

Y2way(:,1) = random('normal',m1,v1,[nitem*Cond2,1]);
Y2way(:,2) = random('normal',m2,v2,[nitem*Cond2,1]);
Y2way(:,3) = random('normal',m3,v3,[nitem*Cond2,1]);
[p,tbl,stats] = anova2(Y2way,ntime);
%% comparison amount multiple condition
figure;
tbl
[c,m,h,nms] = multcompare(stats,'Estimate','column')
%% using linear model class in matlab
response=[Y2way(:,1);Y2way(:,2);Y2way(:,3)];
condi1=[ones(ntime*Cond2,1);ones(ntime*Cond2,1)*2;ones(ntime*Cond2,1)*3];
condi2=repmat([ones(ntime,1);ones(ntime,1)*2],[3,1]);
datM=dataset(response,condi1,condi2);
datM.condi1=nominal(datM.condi1,{'a','b','c'});
datM.condi2=nominal(datM.condi2,{'level1','level2'});

disp(datM)

lm=fitlm(datM,'response~condi1+condi2+condi1:condi2')
anova(lm)

%% linear model again
DX=kron(eye(6),ones(nitem,1));
lm2=fitlm(DX,response,'Intercept',false);
beta=lm2.Coefficients.Estimate;
covm=lm2.CoefficientCovariance;
dfehlm=length(response)-rank(DX);
C{1}=[1 1 -1 -1 0 0;.5 .5 .5 .5 -1 -1];
C{2}=[1 -1 1 -1 1 -1];
C{3}=[1 -1 -1 1 0 0;.5 -.5 .5 -.5 -1 1];
clear pl Fl df1l df2l ph Fh df1h df2h
for ic=1:3
c=C{ic};
[pl(ic),Fl(ic),df1l(ic)]=coefTest(lm2,c);
df1h(ic)=rank(c);
df2h(ic)=dfehlm;%(Ns-1)*rank(c);
Fh(ic)=((c*beta)'*((c*covm*c')^-1)*(c*beta))./rank(c);
ph(ic)=1-fcdf(Fh(ic),df1h(ic),df2h(ic));
end
Ftable1=dataset(Fl',df1l',df2h',pl','Varnames',{'Fstat','DF1','DF2','pValue'},'ObsNames',{'Condi1','Condi2','Interation'})
Ftable2=dataset(Fh',df1h',df2h',ph','Varnames',{'Fstat','DF1','DF2','pValue'},'ObsNames',{'Condi1','Condi2','Interation'})
