%% Introduction
% Statistics is the study of the collection, analysis, interpretation,
% presentation, and organization of data (Wikipedia).
%
% Inference: deductive reasoning whenever enough information is at hand to
% permit it; inductive or plausible reasoning when – as is almost
% invariably the case in real problems – the necessary information is not
% available. But if a problem can be solved by deductive reasoning,
% probability theory is not needed for it; thus our topic is the optimal
% processing of incomplete information.
% 
% Modern statistics are based on probability theory, however, there are a
% gap/distinction between 'probability theory' and 'statistical inference'
% 
%% Let's start with some sampling problem:
rng(1);
nitem=1000;
x = randn(nitem,1);% generate n random sample from a normal distribution
% these n samples are random realization of the underlying distribution
% alternatively, you can view it as a random sampling of the same
% distribution for n times

% x = random('exp',1,[nitem,1]);

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
%% biased and unbiased estimation of variance
% recall that from the last part the variance computed by matlab is 
% var_x_2=sum((x-sum(x)/length(x)).^2)/(length(x)-1);
% so why divided by length(x)-1

% Bessel correction
% A partially nice explanation could be found:
% https://liorpachter.wordpress.com/2014/05/25/bessels-correction-and-the-dangers-of-moocs/
% variance can be rewritten as a double sum corresponding to the average
% squared distance between points and the diagonal terms of the sum are
% zero in expectation. An explanation of why this fact leads to the n-1 in
% the unbiased estimator is as follows:

% The first step is to notice that the variance of a random variable is
% equal to half of the expected squared difference of two independent
% identically distributed random variables of that type. Specifically, the
% definition of variance is:

% var(X)=E((X-miu)^2) where miu=E(X) Equivalently, var(X)=E(X^2)-miu^2 Now
% suppose that Y is another random variable identically distributed to X
% and with X,Y independent. Then
% E((X-Y)^2)=2*var(X)
% E((X-Y)^2)=E((X)^2)+E((Y)^2)-2*E(X)*E(Y)=2*E(X^2)-2*miu^2=2*var(X)

% This identity motivates a rewriting of the (uncorrected) sample variance
% s_n in a way that is computationally less efficient, but mathematically
% more insightful:
dist=[];
k=0;
for i=1:length(x)
    for j=1:length(x)
        k=k+1;
        dist(k)=(x(i)-x(j))^2;
    end
end
var_x_3=sum(dist)/length(dist)/2;
disp(['The biased estimator of var of x is ' num2str(var_x_3)])

% the above will be the same as the population variance (also the biased
% estimation of the variance)
% var_x_2=sum((x-sum(x)/length(x)).^2)/(length(x))

% it's bias because we counted i==j, the distance with itself is zero. To
% visulizaed it, we can put a if statment into the loop, or we can use nanmean:
% tips: NaN is very useful dealing with empty or missing data!!!
dist=nan(length(x),length(x));
for i=1:length(x)
    for j=1:length(x)
        if i~=j
            dist(i,j)=(x(i)-x(j))^2;
        end
    end
end
var_x_4=nansum(dist(:))/sum(~isnan(dist(:)))/2; % unbiased
disp(['The unbiased estimator of var of x is ' num2str(var_x_4)])

% even better, lets use the matlab default function:
dist=pdist2(x,x).^2;
var_x_4=sum(dist(:))/(length(x)*(length(x)-1))/2;
%% Covariance
% we actually have a very good code for computing the covariance from
% above, but we will first do it with the classical way.
y=randn(nitem,1);
figure;
scatterhist(x,y)
cov_xy_1 = cov(x,y);
disp(['The covariance between x and y is ' num2str(cov_xy_1(2,1))])

cov_xy_2 = sum((x-mean(x)).*(y-mean(y)))/(length(x)-1);
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
imagesc(distxy)
