close all
gamman=1.27e9;
P=7.76e4;
Gn=5.62e3;
Rsp=1.28e12;
N=2.14e8;
gammae=1/2.2e-9;
omegar=2*pi*2.65e9;
gammar=1.92e9;
Fp=sqrt(Rsp*P);
Fn=sqrt(Rsp*P+gammae*N);
omega=2*pi*linspace(-1000e9,1000e9,1e7+1);
omega=omega(1:end-1);
% TF=((gamman+1i*omega)*Fp+(Gn*P*Fn))./((omegar+omega-1i*gammar).*(omegar-omega+1i*gammar));
TFabs=sqrt(2*Rsp*((gamman.^2 + omega.^2)+Gn^2*P^2*(1+gammae*N/(Rsp*P))-2*gamman*Gn*P)./(P*((omegar-omega).^2+gammar^2).*((omegar+omega).^2+gammar^2)));
X=hilbert(log(abs(TFabs)));
TF=TFabs.*exp(1i*(-imag(X)));
semilogx(omega/2/pi,mag2db(abs(TF)))
xlim([0.1e9,10e9])
ylim([-145,-106])

%% Fitting using rationfit
s=1j*omega;
freqPosInd=omega>=0;
[fit,errdB] = rationalfit(omega(freqPosInd)/2/pi,TF(freqPosInd),NPoles=3);
TFFreqFitted=sum(fit.C.'./(s.'-fit.A.'),2);
hold on
semilogx(omega/2/pi,mag2db(abs(TFFreqFitted)))
