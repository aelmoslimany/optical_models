function [TECQout,TDECQout,OMAdBm,ER,Ceqout,RLM] = TDECQCalc_coreFn(opticalPowerBeforeFiber,opticalPowerAfterFiber,SSPRQseq,samplesPerUI,EnablePlotEyeDiagram)
%% Initialization
SSPRQLen=length(SSPRQseq);
N=SSPRQLen*samplesPerUI;
nTaps=15;
SSPRQseqSym=reshape(repmat(SSPRQseq.',samplesPerUI,1),[],1);

if ~exist('EnablePlotEyeDiagram','var')
    EnablePlotEyeDiagram=1;
end
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
T3min=(P3ind(1)-7+2.5)*samplesPerUI;
T3max=(P3ind(1)-7+4.5)*samplesPerUI;
T3ind=n>=T3min & n<T3max;
ind=SSPRQseq==SSPRQseqMin;
placesThatHasOf6consective0=conv(ind,ones(6,1));
P0ind=find(placesThatHasOf6consective0==6);
T0min=(P0ind(1)-6+2)*samplesPerUI;
T0max=(P0ind(1)-6+4)*samplesPerUI;
T0ind=n>=T0min & n<T0max;

%% Delay Estimate
timeShift=-100*samplesPerUI:100*samplesPerUI;
% Correct the delay before Fiber
tmp=xcorr(opticalPowerBeforeFiber,SSPRQseqSym,100*samplesPerUI);
[~,ind]=max(abs(tmp));
delay=timeShift(ind);
if(sign(tmp(ind))==-1)
opticalPowerBeforeFiber=circshift(opticalPowerBeforeFiber,-delay);
tmpMean=mean(opticalPowerBeforeFiber);
opticalPowerBeforeFiber=-1*(opticalPowerBeforeFiber-tmpMean)+tmpMean;
else
opticalPowerBeforeFiber=circshift(opticalPowerBeforeFiber,-delay);
end

% Correct the delay after Fiber
tmp=xcorr(opticalPowerAfterFiber,SSPRQseqSym,100*samplesPerUI);
[~,ind]=max(abs(tmp));
delay=timeShift(ind);
if(sign(tmp(ind))==-1)
opticalPowerAfterFiber=circshift(opticalPowerAfterFiber,-delay);
tmpMean=mean(opticalPowerAfterFiber);
opticalPowerAfterFiber=-1*(opticalPowerAfterFiber-tmpMean)+tmpMean;
else
opticalPowerAfterFiber=circshift(opticalPowerAfterFiber,-delay);
end

%% Pavg, OMA, and ER
Pavg=pow2db(mean(opticalPowerBeforeFiber))+30;
P3=mean(opticalPowerBeforeFiber(T3ind));
P0=mean(opticalPowerBeforeFiber(T0ind));
OMAdBm=pow2db(P3-P0)+30; % dBm
ER=pow2db(P3./P0); % dB

%% fourth-order Bessel-Thomson
y=1j*2.114*fOverfr;
besselThomson4thResponse=(105./(105+105*y+45*y.^2+10*y.^3+y.^4)).';
NoisePowerSpectrum=abs(besselThomson4thResponse(fPostiveInd)).^2;
NoisePowerSpectrum=(NoisePowerSpectrum.')./(sum(NoisePowerSpectrum));

%% Filtering the optical power before and after fiber
opticalPowerBeforeFiberFiltered=real(ifft(fft(opticalPowerBeforeFiber).*ifftshift(besselThomson4thResponse)));
opticalPowerAfterFiberFiltered=real(ifft(fft(opticalPowerAfterFiber).*ifftshift(besselThomson4thResponse)));

%% Interpolate if the samplesPerUI is not multiple of 20
if(samplesPerUI>=20)
    samplesPerUITmp=floor(samplesPerUI/20)*20;
else
    samplesPerUITmp=20;
end

% Interpolation of the optical power before the fiber
dataInTmp=[opticalPowerBeforeFiberFiltered;opticalPowerBeforeFiberFiltered(1:1000)];
timeOld=(0:length(dataInTmp)-1).';
time=(0:length(SSPRQseq)*samplesPerUITmp-1).'*samplesPerUI/samplesPerUITmp;
opticalPowerBeforeFiberFiltered=interp1(timeOld,dataInTmp,time,"cubic");
% Interpolation of the optical power before the fiber
dataInTmp=[opticalPowerAfterFiberFiltered;opticalPowerAfterFiberFiltered(1:1000)];
opticalPowerAfterFiberFiltered=interp1(timeOld,dataInTmp,time,"cubic");

%% Update the frequency axis, Ts, and samplesPerUI and prepare the input seq sym
samplesPerUI=samplesPerUITmp;
N=length(SSPRQseq)*samplesPerUI;
n=(0:N-1).';
fOverfr=linspace(-samplesPerUI,samplesPerUI,N+1).';
fOverfr=fOverfr(1:end-1);
fPostiveInd=fOverfr>=0;
fPostiveByTsym=fOverfr(fPostiveInd).'/2;

%% Update OMA measuring times
SSPRQseqMax=max(SSPRQseq);
SSPRQseqMin=min(SSPRQseq);
ind=SSPRQseq==SSPRQseqMax;
placesThatHasOf7consective3=conv(ind,ones(7,1));
P3ind=find(placesThatHasOf7consective3==7);
T3min=(P3ind(1)-7+2.5)*samplesPerUI;
T3max=(P3ind(1)-7+4.5)*samplesPerUI;
T3ind=n>=T3min & n<T3max;
ind=SSPRQseq==SSPRQseqMin;
placesThatHasOf6consective0=conv(ind,ones(6,1));
P0ind=find(placesThatHasOf6consective0==6);
T0min=(P0ind(1)-6+2)*samplesPerUI;
T0max=(P0ind(1)-6+4)*samplesPerUI;
T0ind=n>=T0min & n<T0max;

%% Convolution of the input sequence with rectangular pulse
SSPRQseqSym=reshape(repmat(SSPRQseq.',samplesPerUI,1),[],1);

%% Redefine fourth-order Bessel-Thomson
y=1j*2.114*fOverfr;
besselThomson4thResponse=(105./(105+105*y+45*y.^2+10*y.^3+y.^4)).';
NoisePowerSpectrum=abs(besselThomson4thResponse(fPostiveInd)).^2;
NoisePowerSpectrum=(NoisePowerSpectrum)./(sum(NoisePowerSpectrum));

%% Remove DC 
opticalPowerBeforeFiberFilteredNoDC=opticalPowerBeforeFiberFiltered-mean(opticalPowerBeforeFiberFiltered);
opticalPowerAfterFiberFilteredNoDC=opticalPowerAfterFiberFiltered-mean(opticalPowerAfterFiberFiltered);

%% Iterations to get best TDECQ
TDECQAcceptableRange=(0:0.1:30).';
const1=1./(db2pow(TDECQAcceptableRange)*6*3.414);
opticalPowerFilteredReshaped=reshape(opticalPowerAfterFiberFilteredNoDC,samplesPerUI,[]).';
inputVolt=reshape(SSPRQseqSym,samplesPerUI,[]).';
NMSEValue=zeros(samplesPerUI,1);
leftHistInd=1:(samplesPerUI/2+1);
rightHistInd=leftHistInd+samplesPerUI/(10);
TDECQ=inf*ones(samplesPerUI,1);
LeftSamplingBestTDECQ=zeros(samplesPerUI,1);
CeqdB=zeros(samplesPerUI,1);
for i1=1:samplesPerUI
    Y = hankel(opticalPowerFilteredReshaped(1:nTaps,i1),circshift(opticalPowerFilteredReshaped(:,i1).',-(nTaps-1)));
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
    opticalPowerEqualized=cconv(circshift(opticalPowerAfterFiberFiltered,-(nTaps-NMSEInd)*samplesPerUI+1),w(end:-1:1).',length(opticalPowerAfterFiberFiltered));
    averagePower=mean(opticalPowerEqualized);
    P3=mean(opticalPowerEqualized(T3ind));
    P0=mean(opticalPowerEqualized(T0ind));
    OMA=P3-P0;
    th1=averagePower-OMA/3;
    th2=averagePower;
    th3=averagePower+OMA/3;
    binsLimit=[P0-OMA/2,P3+OMA/2];
    Deltay=OMA/300;
    yi=linspace(P0-OMA/2+Deltay,P3+OMA/2-Deltay,600);
    Heq=wUI*exp(-1j*2*pi*((0:length(wUI)-1).')*fPostiveByTsym);
    Ceq=sqrt(sum(NoisePowerSpectrum.*abs(Heq).^2));
    CeqdB(i1)=pow2db(Ceq);
    sigmaGRange=const1.*OMA;
    Gth1=exp((-((yi-th1)./(Ceq.*sigmaGRange)).^2)/2).*Deltay./(Ceq.*sigmaGRange*sqrt(2*pi));
    Gth2=exp((-((yi-th2)./(Ceq.*sigmaGRange)).^2)/2).*Deltay./(Ceq.*sigmaGRange*sqrt(2*pi));
    Gth3=exp((-((yi-th3)./(Ceq.*sigmaGRange)).^2)/2).*Deltay./(Ceq.*sigmaGRange*sqrt(2*pi));
    %
    opticalPowerEqualizedShifted=circshift(opticalPowerEqualized,-samplesPerUI/2);
    TDECQL=inf*ones(length(leftHistInd),1);
    TDECQR=inf*ones(length(leftHistInd),1);
    for i2=1:length(leftHistInd)
        LHSData=opticalPowerEqualizedShifted(leftHistInd(i2):samplesPerUI:end);
        % tmp=histogram(LHSData,"BinLimits",binsLimit,"BinWidth",OMA/300,"Normalization","probability");
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
            TDECQTmp=TDECQAcceptableRange(indL);
            TDECQL(i2)=TDECQTmp(1);
        else
            TDECQL(i2)=inf;
        end

        RHSData=opticalPowerEqualizedShifted(rightHistInd(i2):samplesPerUI:end);
        % tmp=histogram(RHSData,"BinLimits",binsLimit,"BinWidth",OMA/300,"Normalization","probability");
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
            TDECQTmp=TDECQAcceptableRange(indR);
            TDECQR(i2)=TDECQTmp(1);
        else
            TDECQR(i2)=inf;
        end
    end
    TDECQTmp=max([TDECQL TDECQR],[],2);
    [TDECQ(i1),ind]=min(TDECQTmp);
    LeftSamplingBestTDECQ(i1)=leftHistInd(ind);
end
[TDECQout,i1]=min(TDECQ);
Ceqout=CeqdB(i1);

%% Iterations to get best TECQ
TDECQAcceptableRange=(0:0.1:30).';
const1=1./(db2pow(TDECQAcceptableRange)*6*3.414);
opticalPowerFilteredReshaped=reshape(opticalPowerBeforeFiberFilteredNoDC,samplesPerUI,[]).';
inputVolt=reshape(SSPRQseqSym,samplesPerUI,[]).';
NMSEValue=zeros(samplesPerUI,1);
leftHistInd=1:(samplesPerUI/2+1);
rightHistInd=leftHistInd+samplesPerUI/(10);
TDECQ=inf*ones(samplesPerUI,1);
LeftSamplingBestTDECQ=zeros(samplesPerUI,1);
firFilterTaps=inf*ones(samplesPerUI,nTaps);
CeqdB=zeros(samplesPerUI,1);
for i1=1:samplesPerUI
    Y = hankel(opticalPowerFilteredReshaped(1:nTaps,i1),circshift(opticalPowerFilteredReshaped(:,i1).',-(nTaps-1)));
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
    opticalPowerEqualized=cconv(circshift(opticalPowerBeforeFiberFiltered,-(nTaps-NMSEInd)*samplesPerUI+1),w(end:-1:1).',length(opticalPowerBeforeFiberFiltered));
    averagePower=mean(opticalPowerEqualized);
    P3=mean(opticalPowerEqualized(T3ind));
    P0=mean(opticalPowerEqualized(T0ind));
    OMA=P3-P0;
    th1=averagePower-OMA/3;
    th2=averagePower;
    th3=averagePower+OMA/3;
    binsLimit=[P0-OMA/2,P3+OMA/2];
    Deltay=OMA/300;
    yi=linspace(P0-OMA/2+Deltay,P3+OMA/2-Deltay,600);
    Heq=wUI*exp(-1j*2*pi*((0:length(wUI)-1).')*fPostiveByTsym);
    Ceq=sqrt(sum(NoisePowerSpectrum.*abs(Heq).^2));
    CeqdB(i1)=pow2db(Ceq);
    sigmaGRange=const1.*OMA;
    Gth1=exp((-((yi-th1)./(Ceq.*sigmaGRange)).^2)/2).*Deltay./(Ceq.*sigmaGRange*sqrt(2*pi));
    Gth2=exp((-((yi-th2)./(Ceq.*sigmaGRange)).^2)/2).*Deltay./(Ceq.*sigmaGRange*sqrt(2*pi));
    Gth3=exp((-((yi-th3)./(Ceq.*sigmaGRange)).^2)/2).*Deltay./(Ceq.*sigmaGRange*sqrt(2*pi));
    %
    opticalPowerEqualizedShifted=circshift(opticalPowerEqualized,-samplesPerUI/2);
    TDECQL=inf*ones(length(leftHistInd),1);
    TDECQR=inf*ones(length(leftHistInd),1);
    for i2=1:length(leftHistInd)
        LHSData=opticalPowerEqualizedShifted(leftHistInd(i2):samplesPerUI:end);
        % tmp=histogram(LHSData,"BinLimits",binsLimit,"BinWidth",OMA/300,"Normalization","probability");
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
            TDECQTmp=TDECQAcceptableRange(indL);
            TDECQL(i2)=TDECQTmp(1);
        else
            TDECQL(i2)=inf;
        end

        RHSData=opticalPowerEqualizedShifted(rightHistInd(i2):samplesPerUI:end);
        % tmp=histogram(RHSData,"BinLimits",binsLimit,"BinWidth",OMA/300,"Normalization","probability");
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
            TDECQTmp=TDECQAcceptableRange(indR);
            TDECQR(i2)=TDECQTmp(1);
        else
            TDECQR(i2)=inf;
        end
    end
    TDECQTmp=max([TDECQL TDECQR],[],2);
    [TDECQ(i1),ind]=min(TDECQTmp);
    LeftSamplingBestTDECQ(i1)=leftHistInd(ind);
end
[TECQout,i1]=min(TDECQ);

%% Check if the FIR taps are within the constraints
MinTDECQFIRTaps=firFilterTaps(i1,:);
MinTDECQFIRTaps=MinTDECQFIRTaps./sum(MinTDECQFIRTaps);
[maxTapValue,indOfMaxTap]=max(MinTDECQFIRTaps);
Y = hankel(opticalPowerFilteredReshaped(1:nTaps,i1),circshift(opticalPowerFilteredReshaped(:,i1).',-(nTaps-1)));
X = hankel(inputVolt(1:nTaps,i1),circshift(inputVolt(:,i1).',-(nTaps-1)));
C=Y*(Y');
T=X*(Y');
W=(T/C);
Xhat=W*Y;
NMSE=mean(abs(X-Xhat).^2,2)./mean(abs(Xhat).^2,2);
[NMSEValue(i1),NMSEInd]=min(NMSE);
wUI=W(NMSEInd,:);
wUI=wUI./sum((wUI));
w=upsample(wUI,samplesPerUI,i1-1);
opticalPowerEqualized=cconv(circshift(opticalPowerBeforeFiberFiltered,-(nTaps-NMSEInd)*samplesPerUI+1),w(end:-1:1).',length(opticalPowerBeforeFiberFilteredNoDC));

%% Delay Estimate
tmp=xcorr(opticalPowerEqualized,SSPRQseqSym,100*samplesPerUI);
[~,ind]=max(tmp);
delay=timeShift(ind);
opticalPowerEqualized=circshift(opticalPowerEqualized,-delay);

%% Eye diagram
if EnablePlotEyeDiagram
figure
eyeDiagram = eyeDiagramSI;
eyeDiagram.SymbolTime =1;
eyeDiagram.SampleInterval = 1 / samplesPerUI;
eyeDiagram.Modulation = 4;
eyeDiagram(opticalPowerAfterFiber(1:100000));
plot(eyeDiagram);
title('Before Ref Rx Equalizer');
xlabel('UI')
axis square

figure
eyeDiagram = eyeDiagramSI;
eyeDiagram.SymbolTime =1;
eyeDiagram.SampleInterval = 1 / samplesPerUI;
eyeDiagram.Modulation = 4;
eyeDiagram(opticalPowerEqualized(1:100000));
plot(eyeDiagram);
xlabel('UI')
title('After Ref Rx Equalizer');
axis square
end
%% RLM Computation
opticalPowerEqualizedUndersampled = opticalPowerEqualized(1:samplesPerUI:end);
timeShift=-100*samplesPerUI:100*samplesPerUI;
% Correct the delay before Fiber
tmp=xcorr(opticalPowerEqualizedUndersampled,SSPRQseq,100*samplesPerUI);
[~,ind]=max(tmp);
delay=timeShift(ind);
opticalPowerEqualizedUndersampled=circshift(opticalPowerEqualizedUndersampled,-delay);

ind1=SSPRQseq>2/3;
ind0p3=SSPRQseq<2/3 & SSPRQseq>0;
indn0p3=SSPRQseq>-2/3 & SSPRQseq<0;
indn1=SSPRQseq<-2/3;
L1=mean(opticalPowerEqualizedUndersampled(ind1));
L0p3=mean(opticalPowerEqualizedUndersampled(ind0p3));
Ln0p3=mean(opticalPowerEqualizedUndersampled(indn0p3));
Ln1=mean(opticalPowerEqualizedUndersampled(indn1));
RLM=100*(min([abs(L1-L0p3) abs(L0p3-Ln0p3) abs(Ln0p3-Ln1)])/(abs(L1-Ln1)/3));

end
