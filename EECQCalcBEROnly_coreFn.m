function [BER] = EECQCalcBEROnly_coreFn(Vin,SSPRQseq,samplesPerUI,VRMS,nextIterInit,BERType)
%% Initialization
nTaps=64;

%% Delay correction
Vin=nextIterInit.sign*circshift(Vin,-nextIterInit.delay);


%% Filtering the optical power before and after fiber
VinFiltered=real(ifft(fft(Vin).*nextIterInit.ifftshiftBesselThomson4thResponse));
VinFilteredNoDC=VinFiltered-mean(VinFiltered);
VinFilteredReshaped=reshape(VinFilteredNoDC,samplesPerUI,[]).';
SSPRQseqSym=reshape(repmat(SSPRQseq.',samplesPerUI,1),[],1);
inputVolt=reshape(SSPRQseqSym,samplesPerUI,[]).';
%% 
wUI=nextIterInit.wUI;
NMSEInd=nextIterInit.NMSEInd;
uniqueSyms=sort(unique(SSPRQseq),'ascend');
th=(uniqueSyms(1:end-1)+uniqueSyms(2:end))/2;
ell=samplesPerUI/20;
i1=nextIterInit.samplingInd;
Vpp5Percent=nextIterInit.Vpp5Percent;
whiteningFilterTaps=nextIterInit.whiteningFilterTaps;
whiteNoiseVarM=nextIterInit.whiteNoiseVarM;
%% Left Calc
yL=VinFilteredReshaped(:,mod(i1-ell-1,samplesPerUI)+1)/(Vpp5Percent);% + VRMS*randn(size(VinFilteredReshaped(:,i1)));
YL = hankel(yL(1:nTaps),circshift(yL.',-(nTaps-1)));
xhatL=(wUI*YL).';
xL=inputVolt(:,mod(i1-ell-1,samplesPerUI)+1);
XL = hankel(xL(1:nTaps),circshift(xL.',-(nTaps-1)));
if(i1-ell<1)
    xL=XL(NMSEInd+1,:).';
else
    xL=XL(NMSEInd,:).';
end

BERL=statBERPostProcessing2_coreFn(xhatL,xL,uniqueSyms,th,whiteNoiseVarM,BERType);


%% Right Calc
yR=VinFilteredReshaped(:,mod(i1+ell-1,samplesPerUI)+1)/(Vpp5Percent);% + VRMS*randn(size(VinFilteredReshaped(:,i1)));
YR = hankel(yR(1:nTaps),circshift(yR.',-(nTaps-1)));
xhatR=(wUI*YR).';
xR=inputVolt(:,mod(i1+ell-1,samplesPerUI)+1);
XR = hankel(xR(1:nTaps),circshift(xR.',-(nTaps-1)));
if(i1+ell>samplesPerUI)
    xR=XR(NMSEInd-1,:).';
else
    xR=XR(NMSEInd,:).';
end

BERR=statBERPostProcessing2_coreFn(xhatR,xR,uniqueSyms,th,whiteNoiseVarM,BERType);

%% SNR and BER calc
BER=(BERL+BERR)./2;

end
