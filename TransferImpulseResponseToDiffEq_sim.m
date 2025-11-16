%% Cleaning
clear
close all
clc
%% Initialization
Init.Ts=10e-12/240;%100e-15;
Init.IRFile='Z:\MTKOpticalModel\OpticalModelProj\DataFiles\TxAFE\TxAFEIR3.csv';
Init.Tsym=1/(100e9);
Init.NPole=14;
Init.Oversampling=round(Init.Tsym./Init.Ts);
Init.OmegaHighPass=2*pi*[47e3;48e3;49e3;50e3];
%% Read the IR
IR=readmatrix(Init.IRFile);
time=(0:length(IR)-1)*Init.Ts;
freq=linspace(-1/(2*Init.Ts),1/(2*Init.Ts),1e6+1);
freq=freq(1:end-1).';
s=2j*pi*freq;
IRFreq=fftshift(fft(IR,1e6));
indFreqPos=freq>=0;
freqPos=freq(indFreqPos);
IRFreqPos=IRFreq(indFreqPos);
%% Fit to differential equaiton
[fit,errdB] = rationalfit(freqPos,IRFreqPos,NPoles=Init.NPole);
IRFreqFitted=sum(fit.C.'./(s-fit.A.'),2);
ind1stOrder=imag(fit.C)==0;
omegaPoleTs=zeros(sum(ind1stOrder),1);
omegaPoleTs(:,1)=-fit.A(ind1stOrder).*Init.Ts;
A=-fit.C(ind1stOrder)./fit.A(ind1stOrder);
alphaTs=zeros(sum(~ind1stOrder),1);
betaTs=zeros(sum(~ind1stOrder),1);
alphaTs(:,1)=-real(fit.A(~ind1stOrder)).*Init.Ts;
betaTs(:,1)=-imag(fit.A(~ind1stOrder)).*Init.Ts;
AR=-real(fit.C(~ind1stOrder)./fit.A(~ind1stOrder));
AI=-imag(fit.C(~ind1stOrder)./fit.A(~ind1stOrder));
alphaTs=alphaTs(1:2:end,:);
betaTs=betaTs(1:2:end,:);
AR=AR(1:2:end);
AI=AI(1:2:end);
polesData=[length(omegaPoleTs);A;omegaPoleTs./Init.Ts;length(AR);AR;AI;alphaTs./Init.Ts;betaTs./Init.Ts];
writematrix(polesData,'Z:\MTKOpticalModel\OpticalModelProj\DataFiles\TxAFE\TxAFERationalTF3.csv');
%% Visualization
figure
plot(time,IR)
figure
plot(freq,mag2db(abs(IRFreq)))
hold on
plot(freq,mag2db(abs(IRFreqFitted)))
%% Stress Code
SSPRQseq = importSSPRQ_coreFn("D:\OneDrive - InfiniLink\SerDesModeling\HighSpeedCommunication\channelsData\sequences\SSPRQ_sequence.csv");
SSPRQseqSym = (SSPRQseq-1.5)/1.5;
SSPRQseqSym = reshape(repmat(SSPRQseqSym.',Init.Oversampling,1),[],1);
SSPRQseqSymUndersampled = (SSPRQseq-1.5)/1.5;
freq=linspace(-1/(2*Init.Ts),1/(2*Init.Ts),length(SSPRQseqSym)+1);
freq=freq(1:end-1).';
time = (0:length(SSPRQseqSym)-1)*Init.Ts;
%% Apply the filtering
y1=cconv(SSPRQseqSym,IR,length(SSPRQseqSym));
y1=real(ifft(ifftshift(fftshift(fft(y1)).*exp(1j*pi*freq*Init.Ts))));
y2=linearFiltering_coreFn([SSPRQseqSym;SSPRQseqSym],A,omegaPoleTs,AR,AI,alphaTs,betaTs);
y2=y2(end/2+1:end);
figure
plot(time,y1)
hold on
plot(time,y2)
%% Compute SQNR
SQNRdB=pow2db(var(y1)./var(y1-y2))
%% HighPass
omegaPoleHighTs=Init.OmegaHighPass.*Init.Ts;
AHigh=zeros(length(omegaPoleHighTs),1);
for i=1:length(Init.OmegaHighPass)
    tmp=Init.OmegaHighPass;
    tmp(i)=[];
    AHigh(i)=(prod((1-Init.OmegaHighPass(i)./Init.OmegaHighPass))-prod(-Init.OmegaHighPass(i)./Init.OmegaHighPass))./prod(1-Init.OmegaHighPass(i)./tmp);
end
ADC=1;
[y3] = highpassFiltering_coreFn([y2;y2;y2;y2],AHigh,omegaPoleHighTs,ADC);
y3=y3(3*end/4+1:end);
figure
plot(time,y2)
hold on
plot(time,y3)
figure
plot(time,y3-y2)
%%
[y4] = highpassFiltering_coreFn(ones(4e7,1),AHigh,omegaPoleHighTs,ADC);
figure
plot((0:length(y4)-1)*Init.Ts,y4)
