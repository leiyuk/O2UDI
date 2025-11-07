function y = logn(miu,sigma,x)
y=1./x/log(10)/sqrt(2*pi)/sigma.*exp(-(log10(x)-miu).^2/2/sigma^2);
end

