function int=gaussint(lbound,ubound,alpha,mu)
%lbound and ubound are the bounds of our integration
%alpha is the variance, aka standard deviation squared
%mu is the mean of our gaussian
sigma=alpha.^(1/2);
term=(pi.^(1/4)).*(sigma.^(1/2))/(2.^(1/4));
%term comes from our gaussians being square rooted
lb=(mu-lbound)/(2.*sigma);
ub=(mu-ubound)/(2.*sigma);
int=-term.*(erf(ub)-erf(lb));