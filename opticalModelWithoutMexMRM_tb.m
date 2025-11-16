%% Cleaning
clear
close all
clc

%% Initialization
Init.Tsym=1/(100e9); % Defines the symbol period
Init.BaudRate=1/Init.Tsym; % Sets the baud rate based on the symbol period
Init.Oversampling=20; % Specifies an oversampling factor.
Init.Ts=Init.Tsym/Init.Oversampling; % Defines the sampling time interval based on the symbol period and oversampling factor.
Init.Tr2080=2e-12; % Indicates a rise/fall time in seconds.
Init.Vpp=2; % Sets the peak-to-peak voltage of the input signal.
Init.Offset=0; % Defines the DC offset.
Init.EnableTDECQCalc = true; % Enables the calculation of Transmitter and Dispersion Eye Closure Quaternary (TDECQ) for performance analysis.
Init.TDECQNumFold=3; % the number of TDECQ folds, improving the statistical reliability of the TDECQ calculation.
Init.EnableEECQCalc = true; % Enables the calculation of EECQ.
Init.VRMS = 2e-3/0.3; % The RMS voltage of the noise normalized to signal with Vpp = 1V
Init.EnablePlotEyeDiagram=true; % Plot the eyediagram
Init.BERType = 'Hist2';

%% Add the MATLAB folder path
addpath('../../Code/MATLAB/FunctionsFiles/');
addpath('../../Code/MATLAB/SimFiles/')
%% Read the envirnoment variable of the root path
% Get the value of the environment variable OPTICALMODELPROJROOT
projectRoot = getenv('OPTICALMODELPROJROOT');

if isempty(projectRoot)
    error('Environment variable OPTICALMODELPROJROOT is not set.');
end

% Construct the input and output file paths
inputFilePath = fullfile(projectRoot, 'TMP', 'inputData.csv');
outputFilePath = fullfile(projectRoot, 'TMP', 'outputData.csv');

% Display the file paths (for debugging purposes)
disp(['Input file path: ', inputFilePath]);
disp(['Output file path: ', outputFilePath]);


%% Stress Code
stressCodePath=fullfile(projectRoot, 'Sequences/SSPRQ_sequence.csv');
SSPRQseq = readmatrix(stressCodePath);
SSPRQseqSym = (Init.Vpp/2)*(SSPRQseq-1.5)/1.5 + Init.Offset;
SSPRQseqSym = reshape(repmat(SSPRQseqSym.',Init.Oversampling,1),[],1);
SSPRQseqSym = repmat(SSPRQseqSym,Init.TDECQNumFold,1);
SSPRQseqSymUndersampled = (SSPRQseq-1.5)/1.5;
SSPRQseqSymUndersampled = repmat(SSPRQseqSymUndersampled,Init.TDECQNumFold,1);
freq=linspace(-1/(2*Init.Ts),1/(2*Init.Ts),length(SSPRQseqSym)+1);
freq=freq(1:end-1);
time = (0:length(SSPRQseqSym)-1)*Init.Ts;

%% Gaussian pulse
txDriverTF=exp(-2*(pi*freq*Init.Tr2080/1.6832).^2);
txDriverImpulseResponse= ifft(ifftshift(txDriverTF));
txDriverImpulseResponseCausalTmp=circshift(txDriverImpulseResponse,2*round(Init.Tr2080/Init.Ts));
[txDriverImpulseResponseCausal,~]=enforceCausality_coreFn(txDriverImpulseResponseCausalTmp,1e-13,1e-15);
txDriverImpulseResponseCausal=txDriverImpulseResponseCausal./sum(txDriverImpulseResponseCausal);
gaussianTF=fftshift(fft(txDriverImpulseResponseCausal)).';

%% adding rise time
SSPRQSeqSymWithRiseTime = real(ifft(ifftshift(fftshift(fft(SSPRQseqSym)).*gaussianTF)));
figure
plot(time,SSPRQseqSym)
hold on
plot(time,SSPRQSeqSymWithRiseTime)

%% Configuration File (Mitsubishi)
i=1;
Init.ConfigFileName='OpticalModelConfigEOGFMRM.ini';
configFilePath=fullfile(projectRoot, '/Tutorials/Tutorial10/',Init.ConfigFileName);
%% Optical Model Paramaters
OpticalModelInit = readIniFile_coreFn(configFilePath);

%% MATLAB Implementation of the optical model
tic;
[Vout,TECQ(i),TDECQ(i),OMA(i),ER(i),Ceq(i),RLM(i),EECQ(i),CeqRx(i),RLMRx(i),SNRdB(i),BER(i,:)] = opticalModelV2_coreFn(SSPRQSeqSymWithRiseTime,Init.Ts,Init.Oversampling,OpticalModelInit,SSPRQseqSymUndersampled,Init.EnableTDECQCalc,Init.TDECQNumFold,Init.EnablePlotEyeDiagram,Init.EnableEECQCalc,Init.VRMS,Init.BERType);
toc

%% C++ Implementation of the optical model
writematrix(SSPRQSeqSymWithRiseTime(1:end/Init.TDECQNumFold),inputFilePath);
tic;
if ispc
    executablePath=fullfile(projectRoot, '/Code/CPP/','opticalModelV2NoMex.exe');
    command = sprintf('powershell -Command \"%s %s %s %s \"', sprintf("& \'%s\' ", executablePath),num2str(Init.Ts), sprintf(" \'%s\' ", configFilePath),num2str(Init.Oversampling));
    [status] = system(command);
elseif isunix
    executablePath=fullfile(projectRoot, '/Code/CPP/','opticalModelV2NoMex');
    command = sprintf('%s %s %s %s ',  executablePath, num2str(Init.Ts), configFilePath, num2str(Init.Oversampling));
    [status] = system(command);
end
VoutC=readmatrix(outputFilePath);
toc

%% Check
L=length(VoutC);
SMNR(i)=-pow2db(sum(abs(VoutC(end/2+1:end)-Vout(L/2+1:length(VoutC))).^2)./sum(abs(Vout).^2));
BERTheoritical(i)=pamBerAnalyticalComputation_coreFn(4,db2pow(SNRdB(i)));
disp(['Signal-to-Modeling-Noise = ' num2str(SMNR(i)) ' dB'])
disp(['OMA = ' num2str(OMA(i)) ' dBm'])
disp(['ER = ' num2str(ER(i)) ' dB'])
disp(['TECQ = ' num2str(TECQ(i)) ' dB'])
disp(['TDECQ = ' num2str(TDECQ(i)) ' dB'])
disp(['Ceq = ' num2str(Ceq(i)) ' dB'])
disp(['RLM = ' num2str(RLM(i)) ' %'])
disp(['SNR = ' num2str(SNRdB(i)) ' dB'])
disp(['BER = ' num2str(sum(BER(i,:))) ' bit/sec'])
disp(['Approximate BER = ' num2str(sum(BERTheoritical(i))) ' bit/sec (Calculated using the erfc relation and SNR)'])
