function [BER] = statBERPostProcessing_coreFn(xhat,x,whiteningFilterTaps,whiteNoiseVar,uniqueSyms,th)
%% Pre allocate variables
BER=zeros(1,length(uniqueSyms));
nBitsPerSymbol=log2(length(uniqueSyms));
%% Calculate Noise after removing the white pare
noise=xhat-x;
whiteNoise=cconv(noise,whiteningFilterTaps,length(noise));
remainingNoise=noise-whiteNoise;
xhat=x+remainingNoise;

%% Statstical BER calculation
ind=x==uniqueSyms(1);
levelsMinusTh1=th(1)-xhat(ind);
BER(1)=sum(qfunc(levelsMinusTh1/sqrt(whiteNoiseVar)))/(nBitsPerSymbol*length(x));
for k=2:length(th)
    ind=x==uniqueSyms(k);
    levelsMinusTh1=th(k)-xhat(ind);
    levelsMinusTh2=xhat(ind)-th(k-1);
    BER(k)=sum(qfunc(levelsMinusTh1/sqrt(whiteNoiseVar)))/(length(x)) + sum(qfunc(levelsMinusTh2/sqrt(whiteNoiseVar)))/(nBitsPerSymbol*length(x));
end
ind=x==uniqueSyms(end);
levelsMinusTh2=xhat(ind)-th(end);
BER(end)=sum(qfunc(levelsMinusTh2/sqrt(whiteNoiseVar)))/(nBitsPerSymbol*length(x));

end
