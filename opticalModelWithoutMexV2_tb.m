%% Cleaning
clear
close all
clc
%% Initialization
Init.Tsym=1/(70e9);
Init.BaudRate=1/Init.Tsym;
Init.Oversampling=20;
Init.Ts=Init.Tsym/Init.Oversampling;
Init.Tr2080=2e-12;
Init.Vpp=1;
Init.Offset=0;
filePath='D:\OneDrive - InfiniLink\MAIN\Engineering\Projects\MTK_OpticalEngineModeling\Deliveries\ThirdDelivery\Config\OpticalModelConfigOE.ini';
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

%% Optical Model Paramaters
OpticalModelInit = readIniFile(filePath);

%% MATLAB Implementation of the optical model
tic;
[Vout,TECQ,TDECQ,OMA,ER,Ceq,RLM] = opticalModelV2_coreFn(SSPRQSeqSymWithRiseTime,Init.Ts,Init.Oversampling,OpticalModelInit,SSPRQseqSymUndersampled);
toc
%% C++ Implementation of the optical model
tic;
writematrix(SSPRQSeqSymWithRiseTime,'./InputFile.csv');
command = sprintf('%s %s %s', 'opticalModelMainV2', sprintf('"%s"', filePath),num2str(Init.Ts));
[status] = system(command);
VoutC=readmatrix(SSPRQSeqSymWithRiseTime,'./OutputFile.csv');
toc
%% Check
NMSE=sum(abs(VoutC-Vout).^2)./sum(abs(Vout).^2)
