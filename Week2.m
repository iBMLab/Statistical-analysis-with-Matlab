% Week 2
% This week is about covariance, correlation and regression between 2
% random variables.
% record the computation of mean and variance from last week:
%%
rng(1);
nitem=1000;
% x = randn(nitem,1);% generate n random sample from a normal distribution
% these n samples are random realization of the underlying distribution
% alternatively, you can view it as a random sampling of the same
% distribution for n times

x = random('exp',1,[nitem,1]);

hist(x,100);% histogram for visualizing the data
hold on;

mean_x_1=mean(x);% compute the mean, the expected value of x
line([mean_x_1,mean_x_1],[0,40],'color','r','linewidth',1.5)
mean_x_2=sum(x)/length(x);% compute the mean, again
disp(['The mean of x is ' num2str(mean_x_1)])

var_x_1=var(x);% compute the variance, the expected value of (x-E(x))^2
line([mean_x_1-var_x_1/2,mean_x_1+var_x_1/2],[20,20],'color','r','linewidth',1.5)
var_x_2=sum((x-sum(x)/length(x)).^2)/(length(x)-1);% compute the variance, again
disp(['The var of x is ' num2str(var_x_1)])
%
%% we can simulate a demo for the central limit theorm 
% x = random('norm',0,1,[nitem,10000]);
x = random('exp',1,[nitem,10000]);
meanx=mean(x,1);
figure;
hist(meanx,100);% histogram for visualizing the data

%% different way to compute covariance 
x=randn(nitem,1);
y=randn(nitem,1);
figure;
scatterhist(x,y)
cov_xy_1 = cov(x,y);
disp(['The covariance between x and y is ' num2str(cov_xy_1(2,1))])

cov_xy_2 = sum((x-mean(x)).*(y-mean(y)))/(nitem-1);
disp(['The covariance between x and y is ' num2str(cov_xy_2)])

% similar to the above
distxy=nan(length(x),length(x));
for i=1:length(x)
    for j=1:length(x)
        if i~=j
            distxy(i,j)=(x(i)-x(j))*(y(i)-y(j));
        end
    end
end
cov_xy_3 = nansum(distxy(:))/sum(~isnan(distxy(:)))/2; % unbiased
disp(['The covariance between x and y is ' num2str(cov_xy_3)])
% imagesc(distxy)
%% correlation
corr(x,y)
cov_xy_1(2,1)./(sqrt(cov_xy_1(1,1))*sqrt(cov_xy_1(2,2)))
%% regression demo
rng(1)
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
% covm=cov(X,ytmp);
% rho2=covm(2,1)./(sqrt(covm(1,1))*sqrt(covm(2,2)));
% [beta,bci,rho2,RINT,STATS] = regress(ytmp,[X,ones(length(X),1)]);
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