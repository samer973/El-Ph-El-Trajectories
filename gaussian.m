function gauss=gaussian(s,alpha,mu)
%outputs squareroot of gaussian
%mu is mean
%alpha is the variance, aka square of standard deviation
gauss=((1/(2*pi.*alpha)).^(1/4)).*exp(-((s-mu).^2)/(4.*alpha));

end