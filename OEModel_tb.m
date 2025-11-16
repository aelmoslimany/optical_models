%% Cleaning
clear
close all
clc
%% Initialization
Init.Tsym=1/(100e9);
Init.BaudRate=1/Init.Tsym;
Init.Oversampling=24;
Init.OversamplingIdeal=120;
Init.Ts=Init.Tsym/Init.Oversampling;
Init.TsIdeal=Init.Tsym/Init.OversamplingIdeal;
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
%% Ideal Optical Mod
[pout] = powerModulator_coreFn(SSPRQSeqSymWithRiseTime,Init.Vpp,db2pow(15),db2pow(4.2-30));
%% Optical Model Paramaters
% TxAFE
Init.TxAFE.Gain=1;
Init.TxAFE.nStages=3;
Init.TxAFE.ClippingRatio=2;
Init.TxAFE.VppReq=1;
Init.TxAFE.NonlinearCoeffCell = {[1.01, -0.05366, -0.06297];[0.9906, 0.6077, -6.079, 7.45, -2.619];[0.9984,-0.05206]};
Init.TxAFE.LinearIRCell=cell(Init.TxAFE.nStages,1);
for i=1:Init.TxAFE.nStages
tmp=readmatrix(['D:\OneDrive - InfiniLink\MAIN\Engineering\Projects\MTK_OpticalEngineModeling\Deliveries\ThirdDelivery\ImpulseResponseFiles\TxAFE\TxAFEIR' num2str(i) '.csv']);
Init.TxAFE.LinearIRCell{i}=tmp;
end

% Laser
Init.Laser.Wavelength = 1304; % nm
Init.Laser.Power = 10; % dBm
Init.Laser.RIN = -inf; % dBc/Hz
Init.Laser.CouplingLosses = 0; % dB

% E/O
Init.EO.Type = "EML";
Init.EO.IRFile = "D:\OneDrive - InfiniLink\MAIN\Engineering\Projects\MTK_OpticalEngineModeling\Deliveries\ThirdDelivery\ImpulseResponseFiles\EML\EML1.csv";
Init.EO.NLPolyCoeff = [2.7341  ,  5.4253  ,  6.9409  ,  4.8986,    1.4972  ,  0.1623];
Init.EO.Vbias = -1.6;
Init.EO.CouplingLosses = 0;
Init.EO.LinearIR=readmatrix(Init.EO.IRFile);

% Fiber
Init.Fiber.Length = 2; % Km 
Init.Fiber.Losses = 4; % dB
Init.Fiber.DispersionSlope = 0.092; % ps/nm^2.Km
Init.Fiber.ZeroDispersionWavelength = 1324; % nm
Init.Fiber.DGD = 2.3; % ps
Init.Fiber.FastAxisRatio = 0.5;

% O/E
Init.OE.Type = "GFPD";
Init.OE.DCbias=-1;
Init.OE.PD_width = 0.6e-6;%%PD width (m)
Init.OE.Temperature = 80;%%Temperature (C)
Init.OE.ZinIRFile='D:\OneDrive - InfiniLink\MAIN\Engineering\Projects\MTK_OpticalEngineModeling\Deliveries\ThirdDelivery\ImpulseResponseFiles\PD\ZinIR.csv';
Init.OE.Responsivity = 0.9;
Init.OE.ZinIR=readmatrix(Init.OE.ZinIRFile);
Init.OE.EnableShotNoise= false;
Init.OE.EnableDarkNoise= false;
tmp=readmatrix('D:\OneDrive - InfiniLink\MAIN\Engineering\Projects\MTK_OpticalEngineModeling\Deliveries\ThirdDelivery\ImpulseResponseFiles\PD\ZinPDEParams.csv');
n1stOrder=tmp(1);
k=2;
Init.OE.A=tmp(k:k+n1stOrder-1);
k=k+n1stOrder;
Init.OE.omegaPoleTs=tmp(k:k+n1stOrder-1)*Init.Ts;
k=k+n1stOrder;
n2ndOrder=tmp(k);
k=k+1;
Init.OE.AR=tmp(k:k+n2ndOrder-1);
k=k+n2ndOrder;
Init.OE.AI=tmp(k:k+n2ndOrder-1);
k=k+n2ndOrder;
Init.OE.alphaTs=tmp(k:k+n2ndOrder-1)*Init.Ts;
k=k+n2ndOrder;
Init.OE.betaTs=tmp(k:k+n2ndOrder-1)*Init.Ts;

% Rx AFE
Init.RxAFE.EnableTIANoise= true;
Init.RxAFE.Gain=1;
Init.RxAFE.ClippingRatio=2;
Init.RxAFE.VppReq=0.36;
Init.RxAFE.nStages=3;
Init.RxAFE.LinearIRCell{1}=readmatrix('D:\OneDrive - InfiniLink\MAIN\Engineering\Projects\MTK_OpticalEngineModeling\Deliveries\ThirdDelivery\ImpulseResponseFiles\RxAFE\RxAFEIR1.csv');
Init.RxAFE.LinearIRCell{3}=readmatrix('D:\OneDrive - InfiniLink\MAIN\Engineering\Projects\MTK_OpticalEngineModeling\Deliveries\ThirdDelivery\ImpulseResponseFiles\RxAFE\RxAFEIR3.csv');
Init.RxAFE.NonlinearCoeffCell={[1.0,-31.3371, 6847.755],[1.0,-0.03431,-4.82],[1.003,-0.1188, 0.1453]};
Init.RxAFE.OmegaTs=[0.0060, 0.0; 0.0020, 0.0];
Init.RxAFE.A=[0.4173; 0.0207];
Init.RxAFE.AlphaTs=[0.0225, 0.0; 0.0126, 0.00060822; 0.0074, -0.0022; 0.0123, 0.00011063; 0.0169, -0.00020688];
Init.RxAFE.BetaTs=[-0.0911, 0.0; -0.0410, -0.0050; -0.0410, -0.00099149; -0.0228, 0.00051085; -0.0145, 0.00054284];
Init.RxAFE.AR=[0.0022;2.2477;-0.7260;-12.1032;12.7080];
Init.RxAFE.AI=[-0.0039;0.5276;-2.2890;7.8365;-13.4718];
%% MATLAB Implementation of the OE
% OpticalElectricFieldFastAxisOut=repmat(sqrt(SSPRQSeqSymWithRiseTime+0.53),1,1)*4e-2+1j*eps;
% OpticalElectricFieldSlowAxisOut=OpticalElectricFieldFastAxisOut;
OpticalElectricFieldFastAxisOut=sqrt(pout/2)+1j*eps;
OpticalElectricFieldSlowAxisOut=sqrt(pout/2)+1j*eps;
tic;
[Iout] = OE_coreFn(OpticalElectricFieldFastAxisOut,OpticalElectricFieldSlowAxisOut,Init.Ts,Init.Oversampling,Init.OE,Init.Laser.Wavelength,Init.OE.ZinIR,Init.OE.EnableShotNoise,Init.OE.EnableDarkNoise);
toc
%% C++ Implementation of the fiber
% tic;
% [IoutC] = OEV2_mex(OpticalElectricFieldFastAxisOut,OpticalElectricFieldSlowAxisOut,filePath,Init.Ts,Init.Oversampling);
% toc
%% Check
% NMSE=sum(abs(IoutC-Iout).^2)./sum(abs(Iout).^2)

%% Check the oversampling ratio
SSPRQseq = importSSPRQ_coreFn("D:\OneDrive - InfiniLink\MAIN\Engineering\Projects\MTK_OpticalEngineModeling\Deliveries\ThirdDelivery\Sequences\SSPRQ_sequence.csv");
SSPRQseqSym = (Init.Vpp/2)*(SSPRQseq-1.5)/1.5 + Init.Offset;
SSPRQseqSym = reshape(repmat(SSPRQseqSym.',Init.OversamplingIdeal,1),[],1);
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
%%
OpticalElectricFieldFastAxisOut=repmat(sqrt(SSPRQSeqSymWithRiseTime+0.53),1,1)*4e-2+1j*eps;
OpticalElectricFieldSlowAxisOut=OpticalElectricFieldFastAxisOut;
tic;
[IoutIdeal] = OEV2_mex(OpticalElectricFieldFastAxisOut,OpticalElectricFieldSlowAxisOut,filePath,Init.TsIdeal,Init.OversamplingIdeal);
toc
ratio=round(Init.OversamplingIdeal./Init.Oversampling);

figure
hold on
plot((0:length(IoutC)-1)*ratio,IoutC)
plot((0:length(IoutIdeal)-1),IoutIdeal)

x=zeros(size(IoutIdeal));
for i=1:ratio
OpticalElectricFieldFastAxisOutTmp=OpticalElectricFieldFastAxisOut(i:ratio:end);
OpticalElectricFieldSlowAxisOutTmp=OpticalElectricFieldSlowAxisOut(i:ratio:end);
[IoutCTmp] = OEV2_mex(OpticalElectricFieldFastAxisOutTmp,OpticalElectricFieldSlowAxisOutTmp,filePath,Init.Ts,Init.Oversampling);
x(i:ratio:end)=IoutCTmp;
end




L=(-100:1:100)*ratio;
for i=1:length(L)
Ioutshifted= real(ifft(ifftshift(fftshift(fft(IoutIdeal)).*exp(1j*L(i)*2*pi*Init.TsIdeal.*freq.'))));
NMSEOversampling(i)=sum(abs(Ioutshifted(1000:end-1000)-x(1000:end-1000)).^2)./sum(abs(x(1000:end-1000)).^2);
end
figure
plot(pow2db(NMSEOversampling))




