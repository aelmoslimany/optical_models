function [wEqOptimal,MSEInd,whiteningFilterTaps,whiteNoiseVarM] = noiseySignalAnalysisForStat_coreFn(noiseySignalM,IdealSignal,nTaps,VRMS)
%% Find Optimal Filter
Y = hankel(noiseySignalM(1:nTaps),circshift(noiseySignalM.',-(nTaps-1)));
X = hankel(IdealSignal(1:nTaps),circshift(IdealSignal.',-(nTaps-1)));
C=Y*(Y');
T=X*(Y');
W=(T/C);
Xhat=W*Y;
MSE=mean(abs(X-Xhat).^2,2);
[~,MSEInd]=min(MSE);
wUI=W(MSEInd,:);
wEqOptimal=wUI;
noiseM=Xhat(MSEInd,:)-X(MSEInd,:);
whiteningFilterTaps=1;
% whiteNoiseVarM=var(noiseM);
whiteNoiseVarM=sum(abs(VRMS*wUI).^2);
%% Apply Filter to right sampled signal
%Y = hankel(noiseySignalR(1:nTaps),circshift(noiseySignalR.',-(nTaps-1)));
%xhatR=(wUI*Y).';
%noiseR=xhatR-X(MSEInd,:);
%% Apply Filter to left sampled signal
%Y = hankel(noiseySignalL(1:nTaps),circshift(noiseySignalL.',-(nTaps-1)));
%xhatL=(wUI*Y).';
%noiseL=xhatL-X(MSEInd,:);
















% %% Compute the Autocorrelation of the noise
% RMM=autocorr(noiseM,'NumLags',4*nTaps-1);
% RMM=RMM.'*var(noiseM);
% RMM=[RMM;RMM(end-1:-1:2)];
% RMMFreq=real(fft(RMM));
% [whiteNoiseVarM,ind]=min(abs(RMMFreq));
% whiteNoiseVarM=0.9*whiteNoiseVarM;
% CMMFreq=RMMFreq-whiteNoiseVarM;
% % hc=real(ifft(sqrt(CMMFreq)));
% % w=real(ifft(1./sqrt(RMMFreq)));
% % ww=cconv(w,hc,length(w));
% ww=real(ifft(sqrt(whiteNoiseVarM./RMMFreq)));
% % ww=real(ifft(1./sqrt(RMMFreq)));
% % www=real(ifft(sqrt(RMMFreq-whiteNoiseVarM)));
%
% 
% %% Compute the whitening filter
% % RMM= toeplitz(RMM,[RMM(1); RMM(end:-1:2)]);
% % W=pinv(sqrtm(RMM));
% % w=W(:,1);
% whiteningFilterTaps=[ww(1:end/2);zeros(length(noiseM)-length(ww),1);ww(end/2+1:end)];
% % zz=[www(1:end/2);zeros(length(noiseM)-length(www),1);www(end/2+1:end)];
% %% Apply the whitening filter to the noise
% whiteNoiseM=cconv(noiseM,whiteningFilterTaps,length(noiseM));
% % coloredNoiseM2=cconv(coloredNoiseM,zz,length(noiseM));
% % whiteNoiseM=noiseM-coloredNoiseM2;
% %whiteNoiseR=cconv(noiseR,whiteningFilterTaps,length(noiseR));
% %whiteNoiseL=cconv(noiseL,whiteningFilterTaps,length(noiseL));
% %% Compute the variance of the white noise
% whiteNoiseVarMTmp=var(whiteNoiseM);
% [a,lags]=autocorr(whiteNoiseM,'NumLags',4*nTaps-1);
% a=[a.';a(end-1:-1:2).']*var(whiteNoiseM);
% aFreq=real(fft(a));
% %whiteNoiseVarR=var(whiteNoiseR);
% %whiteNoiseVarL=var(whiteNoiseL);
end
