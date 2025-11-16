function [EECQout,Ceqout,RLM,SNRdB,BER,nextIterInit] = EECQCalc2_coreFn(Vin,SSPRQseq,samplesPerUI,VRMS,BERType)
%% Initialization
SSPRQLen=length(SSPRQseq);
N=SSPRQLen*samplesPerUI;
% nTaps=32;
nTaps=64;
SSPRQseqSym=reshape(repmat(SSPRQseq.',samplesPerUI,1),[],1);

%% Freq and Discrete Time Axis
n=(0:N-1).';
fOverfr=linspace(-samplesPerUI,samplesPerUI,N+1);
fOverfr=fOverfr(1:end-1);
fPostiveInd=fOverfr>=0;
fPostiveByTsym=fOverfr(fPostiveInd)/2;

%% OMA measuring times
SSPRQseqMax=max(SSPRQseq);
SSPRQseqMin=min(SSPRQseq);
ind=SSPRQseq==SSPRQseqMax;
placesThatHasOf7consective3=conv(ind,ones(7,1));
P3ind=find(placesThatHasOf7consective3==7);
T3min=(P3ind-7+2.5)*samplesPerUI;
T3max=(P3ind-7+4.5)*samplesPerUI;
T3ind=n>=T3min(1) & n<T3max(1);
for k=2:length(T3min)
    T3ind=T3ind | (n>=T3min(k) & n<T3max(k));
end
ind=SSPRQseq==SSPRQseqMin;
placesThatHasOf6consective0=conv(ind,ones(6,1));
P0ind=find(placesThatHasOf6consective0==6);
T0min=(P0ind-6+2)*samplesPerUI;
T0max=(P0ind-6+4)*samplesPerUI;
T0ind=n>=T0min(1) & n<T0max(1);
for k=2:length(T0max)
    T0ind=T0ind | (n>=T0min(k) & n<T0max(k));
end

%% Mean Est and Delay Estimate
Vin=Vin-mean(Vin);
timeShift=-100*samplesPerUI:100*samplesPerUI;
% Correct the delay after Fiber
tmp=xcorr(Vin,SSPRQseqSym,100*samplesPerUI);
[~,ind]=max(abs(tmp));
delay=timeShift(ind);
Vin=sign(tmp(ind))*circshift(Vin,-delay);
nextIterInit.delay=delay;
nextIterInit.sign=sign(tmp(ind));

%% Compute Vpp
V3=mean(Vin(T3ind));
V0=mean(Vin(T0ind));
Vpp=V3-V0; % dBm
%% fourth-order Bessel-Thomson
y=1j*2.114*fOverfr;
besselThomson4thResponse=(105./(105+105*y+45*y.^2+10*y.^3+y.^4)).';
NoisePowerSpectrum=abs(besselThomson4thResponse(fPostiveInd)).^2;
NoisePowerSpectrum=(NoisePowerSpectrum.')./(sum(NoisePowerSpectrum));

%% Filtering the optical power before and after fiber
VinFiltered=real(ifft(fft(Vin).*ifftshift(besselThomson4thResponse)));
VinFilteredNoDC=VinFiltered-mean(VinFiltered);
nextIterInit.ifftshiftBesselThomson4thResponse=ifftshift(besselThomson4thResponse);
%% Iterations to get best TDECQ
EECQAcceptableRange=(0:0.1:60).';
const1=1./(db2mag(EECQAcceptableRange)*6*3.414);
VinFilteredReshaped=reshape(VinFilteredNoDC,samplesPerUI,[]).';
inputVolt=reshape(SSPRQseqSym,samplesPerUI,[]).';
NMSEValue=zeros(samplesPerUI,1);
leftHistInd=1:(samplesPerUI/2+1);
rightHistInd=leftHistInd+samplesPerUI/(10);
EECQ=inf*ones(samplesPerUI,1);
LeftSamplingBestTDECQ=zeros(samplesPerUI,1);
firFilterTaps=inf*ones(samplesPerUI,nTaps);
CeqdB=zeros(samplesPerUI,1);
for i1=1:samplesPerUI
    Y = hankel(VinFilteredReshaped(1:nTaps,i1),circshift(VinFilteredReshaped(:,i1).',-(nTaps-1)));
    X = hankel(inputVolt(1:nTaps,i1),circshift(inputVolt(:,i1).',-(nTaps-1)));
    C=Y*(Y');
    T=X*(Y');
    W=(T/C);
    Xhat=W*Y;
    NMSE=mean(abs(X-Xhat).^2,2);
    [NMSEValue(i1),NMSEInd]=min(NMSE);
    wUI=W(NMSEInd,:);
    wUI=wUI./sum((wUI));
    w=upsample(wUI,samplesPerUI,i1-1);
    firFilterTaps(i1,:)=wUI(end:-1:1);
    VinEqualized=cconv(circshift(VinFiltered,-(nTaps-NMSEInd)*samplesPerUI+1),w(end:-1:1).',length(VinFiltered));
    VDC=mean(VinEqualized);
    V3=mean(VinEqualized(T3ind));
    V0=mean(VinEqualized(T0ind));
    Vpp=V3-V0;
    th1=VDC-Vpp/3;
    th2=VDC;
    th3=VDC+Vpp/3;
    binsLimit=[V0-Vpp/2,V3+Vpp/2];
    Deltay=Vpp/300;
    yi=linspace(V0-Vpp/2+Deltay,V3+Vpp/2-Deltay,600);
    Heq=wUI*exp(-1j*2*pi*((0:length(wUI)-1).')*fPostiveByTsym);
    Ceq=sqrt(sum(NoisePowerSpectrum.*abs(Heq).^2));
    CeqdB(i1)=pow2db(Ceq);
    sigmaGRange=const1.*Vpp;
    Gth1=exp((-((yi-th1)./(Ceq.*sigmaGRange)).^2)/2).*Deltay./(Ceq.*sigmaGRange*sqrt(2*pi));
    Gth2=exp((-((yi-th2)./(Ceq.*sigmaGRange)).^2)/2).*Deltay./(Ceq.*sigmaGRange*sqrt(2*pi));
    Gth3=exp((-((yi-th3)./(Ceq.*sigmaGRange)).^2)/2).*Deltay./(Ceq.*sigmaGRange*sqrt(2*pi));
    %
    opticalPowerEqualizedShifted=circshift(VinEqualized,-samplesPerUI/2);
    TDECQL=inf*ones(length(leftHistInd),1);
    TDECQR=inf*ones(length(leftHistInd),1);
    for i2=1:length(leftHistInd)
        LHSData=opticalPowerEqualizedShifted(leftHistInd(i2):samplesPerUI:end);
        tmp=histcounts(LHSData,"BinLimits",binsLimit,"BinWidth",Deltay,"Normalization","probability");
        if(length(tmp)>600)
            tmp(600)=sum(tmp(600:end));
            tmp=tmp(1:600);
        end
        CFL1=[flip(cumsum(tmp(200:-1:1))) cumsum(tmp(201:end))];
        CFL2=[flip(cumsum(tmp(300:-1:1))) cumsum(tmp(301:end))];
        CFL3=[flip(cumsum(tmp(400:-1:1))) cumsum(tmp(401:end))];
        SERL1=CFL1*(Gth1.');
        SERL2=CFL2*(Gth2.');
        SERL3=CFL3*(Gth3.');
        SERL=(SERL1+SERL2+SERL3);
        indL=SERL<4.8e-4;
        if(any(indL))
            TDECQTmp=EECQAcceptableRange(indL);
            TDECQL(i2)=TDECQTmp(1);
        else
            TDECQL(i2)=inf;
        end

        RHSData=opticalPowerEqualizedShifted(rightHistInd(i2):samplesPerUI:end);
        tmp=histcounts(RHSData,"BinLimits",binsLimit,"BinWidth",Deltay,"Normalization","probability");
        if(length(tmp)>600)
            tmp(600)=sum(tmp(600:end));
            tmp=tmp(1:600);
        end
        CFR1=[flip(cumsum(tmp(200:-1:1))) cumsum(tmp(201:end))];
        CFR2=[flip(cumsum(tmp(300:-1:1))) cumsum(tmp(301:end))];
        CFR3=[flip(cumsum(tmp(400:-1:1))) cumsum(tmp(401:end))];
        SERR1=CFR1*(Gth1.');
        SERR2=CFR2*(Gth2.');
        SERR3=CFR3*(Gth3.');
        SERR=(SERR1+SERR2+SERR3);
        indR=SERR<4.8e-4;
        if(any(indR))
            TDECQTmp=EECQAcceptableRange(indR);
            TDECQR(i2)=TDECQTmp(1);
        else
            TDECQR(i2)=inf;
        end
    end
    TDECQTmp=max([TDECQL TDECQR],[],2);
    [EECQ(i1),ind]=min(TDECQTmp);
    LeftSamplingBestTDECQ(i1)=leftHistInd(ind);
end
[EECQout,i1]=min(EECQ);
Ceqout=CeqdB(i1);
nextIterInit.samplingInd=i1;
%% Compute BER for a given VRMS
ell=samplesPerUI/20;
y=VinFilteredReshaped(:,i1);
tmp=prctile(y,[2.5 97.5]);
Vpp5Percent=tmp(2)-tmp(1);
y=VinFilteredReshaped(:,i1)/(Vpp5Percent) + VRMS*randn(size(VinFilteredReshaped(:,i1)));
yR=VinFilteredReshaped(:,mod(i1+ell-1,samplesPerUI)+1)/(Vpp5Percent) + VRMS*randn(size(VinFilteredReshaped(:,i1)));
yL=VinFilteredReshaped(:,mod(i1-ell-1,samplesPerUI)+1)/(Vpp5Percent) + VRMS*randn(size(VinFilteredReshaped(:,i1)));
nextIterInit.Vpp5Percent=Vpp5Percent;
x=inputVolt(:,i1);
[wUI,NMSEInd,whiteningFilterTaps,whiteNoiseVarM] = noiseySignalAnalysisForStat_coreFn(y,x,nTaps,VRMS);

nextIterInit.wUI=wUI;
nextIterInit.whiteningFilterTaps=whiteningFilterTaps;
nextIterInit.whiteNoiseVarM=whiteNoiseVarM;
nextIterInit.NMSEInd=NMSEInd;
uniqueSyms=sort(unique(x),'ascend');
th=(uniqueSyms(1:end-1)+uniqueSyms(2:end))/2;

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
NMSEL=mean(whiteNoiseVarM+abs(xL-xhatL).^2)./(mean(abs(xL).^2)*1.060507730482917); % the factor 1.060507730482917 is because the stress code is not balanced

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
NMSER=mean(whiteNoiseVarM+abs(xR-xhatR).^2)./(mean(abs(xR).^2)*1.060507730482917); % the factor 1.060507730482917 is because the stress code is not balanced

BERR=statBERPostProcessing2_coreFn(xhatR,xR,uniqueSyms,th,whiteNoiseVarM,BERType);

%% SNR and BER calc
SNRdB=pow2db(2./(NMSEL+NMSER));
BER=(BERL+BERR)./2;
%% Check if the FIR taps are within the constraints
MinTDECQFIRTaps=firFilterTaps(i1,:);
MinTDECQFIRTaps=MinTDECQFIRTaps./sum(MinTDECQFIRTaps);
[maxTapValue,indOfMaxTap]=max(MinTDECQFIRTaps);
Y = hankel(VinFilteredReshaped(1:nTaps,i1),circshift(VinFilteredReshaped(:,i1).',-(nTaps-1)));
X = hankel(inputVolt(1:nTaps,i1),circshift(inputVolt(:,i1).',-(nTaps-1)));
C=Y*(Y');
T=X*(Y');
W=(T/C);
Xhat=W*Y;
NMSE=mean(abs(X-Xhat).^2,2)./mean(abs(X).^2,2);
[NMSEValue(i1),NMSEInd]=min(NMSE);
wUI=W(NMSEInd,:);
wUI=wUI./sum((wUI));
w=upsample(wUI,samplesPerUI,i1-1);
VinEqualized=cconv(circshift(VinFiltered,-(nTaps-NMSEInd)*samplesPerUI+1),w(end:-1:1).',length(VinFilteredNoDC));

%% Delay Estimate
tmp=xcorr(VinEqualized,SSPRQseqSym,100*samplesPerUI);
[~,ind]=max(tmp);
delay=timeShift(ind);
VinEqualized=circshift(VinEqualized,-delay);

%% RLM Computation
VinEqualizedUndersampled = VinEqualized(1:samplesPerUI:end);
timeShift=-100*samplesPerUI:100*samplesPerUI;
% Correct the delay before Fiber
tmp=xcorr(VinEqualizedUndersampled,SSPRQseq,100*samplesPerUI);
[~,ind]=max(tmp);
delay=timeShift(ind);
VinEqualizedUndersampled=circshift(VinEqualizedUndersampled,-delay);

ind1=SSPRQseq>2/3;
ind0p3=SSPRQseq<2/3 & SSPRQseq>0;
indn0p3=SSPRQseq>-2/3 & SSPRQseq<0;
indn1=SSPRQseq<-2/3;
L1=mean(VinEqualizedUndersampled(ind1));
L0p3=mean(VinEqualizedUndersampled(ind0p3));
Ln0p3=mean(VinEqualizedUndersampled(indn0p3));
Ln1=mean(VinEqualizedUndersampled(indn1));
RLM=100*(min([abs(L1-L0p3) abs(L0p3-Ln0p3) abs(Ln0p3-Ln1)])/(abs(L1-Ln1)/3));

end
