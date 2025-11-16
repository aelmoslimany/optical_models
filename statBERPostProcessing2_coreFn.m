function [BER] = statBERPostProcessing2_coreFn(xhat,x,uniqueSyms,th,noiseVar,type)
%% Pre allocate variables
BER=zeros(1,length(uniqueSyms));
nBitsPerSymbol=log2(length(uniqueSyms));


switch type

    case 'Hist2'
noiseStd=sqrt(noiseVar);
ind=x==uniqueSyms(1);
BER(1)=sum(qfunc(-(xhat(ind)-th(1))/noiseStd))/(nBitsPerSymbol*length(x));
for k=2:length(th)
    ind=x==uniqueSyms(k);
BER(k)=sum(qfunc(-(xhat(ind)-th(k))/noiseStd))/(nBitsPerSymbol*length(x)) + sum(qfunc(-(-xhat(ind)+th(k-1))/noiseStd))/(nBitsPerSymbol*length(x));
end
ind=x==uniqueSyms(end);
BER(end)=sum(qfunc(-(-xhat(ind)+th(end))/noiseStd))/(nBitsPerSymbol*length(x));

    case 'Hist'
%% Calculate Noise
noise=xhat-x;
p=2.5;
%% Statstical BER calculation
ind=x==uniqueSyms(1);
noiseTmp=noise(ind);
xmin=prctile(noiseTmp,100-p);
[mu_est, sigma_est] = truncatedGaussianMLENewton_coreFn(noiseTmp(noiseTmp>xmin), xmin);
mu_est=mu_est+uniqueSyms(1);
BER(1)=qfunc((-mu_est+th(1))/sigma_est);
for k=2:length(th)
    ind=x==uniqueSyms(k);
    noiseTmp=noise(ind);
    tmp=prctile(noiseTmp,[p 100-p]);
    xmax=tmp(1);
    xmin=tmp(2);
    [mu_est, sigma_est] = truncatedGaussianMLENewton_coreFn(noiseTmp(noiseTmp>xmin), xmin);
    mu_est=mu_est+uniqueSyms(k);
    BER(k)=qfunc((-mu_est+th(k))/sigma_est);
    [mu_est, sigma_est] = truncatedGaussianMLENewton_coreFn(-noiseTmp(noiseTmp<=xmax), -xmax);
    mu_est=-mu_est+uniqueSyms(k);
    BER(k)=BER(k)+qfunc(((-th(k-1)+mu_est))/sigma_est);
end
ind=x==uniqueSyms(end);
noiseTmp=noise(ind);
xmax=prctile(noiseTmp,9);
[mu_est, sigma_est] = truncatedGaussianMLENewton_coreFn(-noiseTmp(noiseTmp<=xmax), -xmax);
mu_est=-mu_est+uniqueSyms(end);
BER(end)=qfunc(((-th(end)+mu_est))/sigma_est);
BER=BER/(nBitsPerSymbol*length(uniqueSyms));


    case 'Real'

ind=x==uniqueSyms(1);
BER(1)=sum(xhat(ind)>th(1))/(nBitsPerSymbol*length(x));
for k=2:length(th)
    ind=x==uniqueSyms(k);
BER(k)=sum(xhat(ind)>th(k) | xhat(ind)<th(k-1))/(nBitsPerSymbol*length(x));
end
ind=x==uniqueSyms(end);
BER(end)=sum(xhat(ind)<th(end))/(nBitsPerSymbol*length(x));


end

end
