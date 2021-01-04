function int=annoygaussint1(s1,t1,s2,t2,alpha,mu)
%lbound and ubound are the bounds of our integration
%alpha is the variance, aka standard deviation squared
%mu is the mean of our gaussian
sigma=alpha.^(1/2);
term=(pi.^(1/4)).*(sigma.^(1/2))/(2.^(5/4));
%term comes from our gaussians being square rooted
%The extra 1/2 comes from our thing being annoying
lb=(t2+s2-mu)/(2.*sigma);
ub=(t1-s1+mu)/(2.*sigma);
int=term.*(erf(ub)+erf(lb));