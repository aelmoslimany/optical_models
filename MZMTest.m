%% Cleaning
clear
close all
clc
%% Initialization
Init.Tsym=1/(100e9);
Init.BaudRate=1/Init.Tsym;
Init.Oversampling=20;
Init.Ts=Init.Tsym/Init.Oversampling;
Init.Tr2080=2e-12;
Init.Vpp=1;
Init.Offset=0;
filePath='D:\OneDrive - InfiniLink\MAIN\Engineering\Projects\MTK_OpticalEngineModeling\Deliveries\ThirdDelivery\Config\OpticalModelConfig1.ini';
%% Stress Code
SSPRQseq = importSSPRQ_coreFn("D:\OneDrive - InfiniLink\MAIN\Engineering\Projects\MTK_OpticalEngineModeling\Deliveries\ThirdDelivery\Sequences\SSPRQ_sequence.csv");
SSPRQseqSym = (Init.Vpp/2)*(SSPRQseq-1.5)/1.5 + Init.Offset;
SSPRQseqSym = reshape(repmat(SSPRQseqSym.',Init.Oversampling,1),[],1);
SSPRQseqSym=repmat(SSPRQseqSym,1,1);
SSPRQseqSymUndersampled = (SSPRQseq-1.5)/1.5;
freq=linspace(-1/(2*Init.Ts),1/(2*Init.Ts),length(SSPRQseqSym)+1);
freq=freq(1:end-1);
time = (0:length(SSPRQseqSym)-1)*Init.Ts;
%% Gaussian pulse
txDriverTF=exp(-2*(pi*freq*Init.Tr2080/1.6832).^2);
txDriverImpulseResponse= ifft(ifftshift(txDriverTF));
txDriverImpulseResponseCausalTmp=circshift(txDriverImpulseResponse,2*round(Init.Tr2080/Init.Ts));
[txDriverImpulseResponseCausal,~]=enforceCausality(txDriverImpulseResponseCausalTmp,1e-13,1e-15);
txDriverImpulseResponseCausal=txDriverImpulseResponseCausal./sum(txDriverImpulseResponseCausal);
gaussianTF=fftshift(fft(txDriverImpulseResponseCausal)).';
figure
plot(freq,mag2db(abs(gaussianTF)))
%% adding rise time
SSPRQSeqSymWithRiseTime = real(ifft(ifftshift(fftshift(fft(SSPRQseqSym)).*gaussianTF)));
figure
plot(time,SSPRQseqSym)
hold on
plot(time,SSPRQSeqSymWithRiseTime)
%% MZM Model Paramaters
HeaterPhase=pi/2;
Vbias=-1;
Zsource = 1*36;
ZLoad = 1*36;
ModulatorLength = 5e-3;
ArmMismatch = 0;%
N_sections = 10;%Number of sections for modelling RF
Wavelength = 1.311e-6;
Opticalpower = 10^(7/10)*1e-3;%Input laser power
%% MZM Model
Vin=SSPRQSeqSymWithRiseTime.';
[Outputpower,OutputField] = MZMmodelfn(Vin,Init.Ts,Vbias,ModulatorLength,N_sections,Zsource,ZLoad,ArmMismatch,Opticalpower,HeaterPhase,Wavelength);


