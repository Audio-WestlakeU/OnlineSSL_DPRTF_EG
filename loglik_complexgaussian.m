function loglik = loglik_complexgaussian(x,mu,sigma)
% x : dprtf estimate, with the size of (No. of freqencies) x (No. of microphone pairs)
% mu : dprtf templete, with the size of (No. of freqencies) x (No. of microphone pairs) x (No. of candidate locations)
% sigma: variance, is set to a constant, i.e. 0.1
%
% loglik: (No. of freqencies) x (No. of microphone pairs) x (No. of candidate locations)
%

D = size(mu,3);
x = repmat(x,[1,1,D]);

% The magnitude of dprtf, i.e. the interchannel level difference, is set to
% smaller than 0.5, to implicitly use the channel that has larger level as
% the reference channel. Note that dprtf magnitude had been normalized to
% [0,1], where 0.5 corresponds to the case that the two channels have the
% same level
upindx = abs(x)>0.5;
x = x.*(1-upindx)+(x.*(1-abs(x))./(abs(x)+eps)).*upindx;
mu = mu.*(1-upindx)+(mu.*(1-abs(mu))./(abs(mu)+eps)).*upindx;

loglik = -log(pi*sigma)-abs(x-mu).^2/sigma;
