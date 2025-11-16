function [chirpSignalTime1Real,chirpSignalTime1Imag,timeAxis,freqAxis,Ts,Fstep] =  fmcwGen_coreFn(BWchirp,OSF,tchirp)
%% FMCW Paramaters
fstart = 0; % Chirp starting Frequency
% tchirp = (2^16-1)*10e-12;%0.5e-7; % Chirp Period
% BWchirp = 400e9; % Bandwidth of the Chirp signal
% OSF = 50; % Oversampling Ratio
Fs=BWchirp*OSF; % Sampling Frequency
Ts=1./Fs; % Sampling Time
N=ceil(tchirp/Ts);
Fstep=1./(Ts*N); % Frequency step
timeAxis=0:Ts:(Ts*(N-1));
freqAxis=-Fs/2:Fstep:Fs/2-Fstep;
Vpp=2;
%% Translate the variable to baseband
fStart= fstart;
fEnd= fstart+BWchirp;
fCarrier=(fStart+fEnd)/2;
fStart=fStart-fCarrier;
fEnd=fEnd-fCarrier;
%% Generate the chirp signal in time and frequency
fChirp=linspace(fStart,fEnd,N);
phaseAngle1=[cumsum(2*pi*fChirp*Ts)];
chirpSignalTime1=(Vpp/2)*exp(1j*phaseAngle1);
chirpSignalFreq1=Ts*fftshift(fft(chirpSignalTime1));
chirpSignalTime1Real=real(chirpSignalTime1);
chirpSignalTime1Imag=imag(chirpSignalTime1);

%% Plot
figure(1)
plot(freqAxis,pow2db(abs(chirpSignalFreq1).^2/(N*Ts)))
