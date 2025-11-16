function [so,a,OmegaTotOfVpn] = ringDynamicsWithHighPowerNonlinearityWithMEX_coreFn(inputVoltage,Ts,omegaBiasInitial,domega_obydV,gammaTot,gammaKappa,selfHeat1stOrderPolesHz,selfHeat1stOrderGains,selfHeat2ndOrderPolesHz,selfHeat2bdOrderGains,domegabybo,domegabybTPA,Pin,IL)
[omegaBiasInitial, ~, aInitial] = mrmRingInitialSoln2_coreFn(domega_obydV, mean(inputVoltage), gammaKappa, Pin, gammaTot, omegaBiasInitial, domegabybo, domegabybTPA, IL);

%% Use the initial condition
c1=-gammaTot/2+1j*omegaBiasInitial;
c2=sqrt(gammaKappa*Pin);
% 1st order LPF paramters
g1 = exp(-2*pi*selfHeat1stOrderPolesHz(:)*Ts);
g2 = (selfHeat1stOrderGains(:).*(1-exp(-2*pi*selfHeat1stOrderPolesHz(:)*Ts)));
% 2nd order LPF paramters
QR=real(selfHeat2bdOrderGains(:)).*(1-exp(-2*pi*real(selfHeat2ndOrderPolesHz(:))*Ts).*cos(2*pi*imag(selfHeat2ndOrderPolesHz(:))*Ts))-imag(selfHeat2bdOrderGains(:)).*exp(-2*pi*real(selfHeat2ndOrderPolesHz(:))*Ts).*sin(2*pi*imag(selfHeat2ndOrderPolesHz(:))*Ts);
QI=imag(selfHeat2bdOrderGains(:)).*(1-exp(-2*pi*real(selfHeat2ndOrderPolesHz(:))*Ts).*cos(2*pi*imag(selfHeat2ndOrderPolesHz(:))*Ts))+real(selfHeat2bdOrderGains(:)).*exp(-2*pi*real(selfHeat2ndOrderPolesHz(:))*Ts).*sin(2*pi*imag(selfHeat2ndOrderPolesHz(:))*Ts);
g3 = 2*exp(-2*pi*real(selfHeat2ndOrderPolesHz(:))*Ts).*cos(2*pi*imag(selfHeat2ndOrderPolesHz(:))*Ts);
g4 = -exp(-4*pi*real(selfHeat2ndOrderPolesHz(:))*Ts);
g5 = 2*QR;
g6 = -2*exp(-2*pi*real(selfHeat2ndOrderPolesHz(:))*Ts).*(QR.*cos(2*pi*imag(selfHeat2ndOrderPolesHz(:))*Ts)-QI.*sin(2*pi*imag(selfHeat2ndOrderPolesHz(:))*Ts));

omegaOfVpn=domega_obydV*inputVoltage(:);
selfHeatGain=selfHeat1stOrderGains(:);
selfHeatGainComplex=2*real(selfHeat2bdOrderGains(:));
[a,OmegaTotOfVpn] = LumpedRingWithHighPowerNonlinearity2_mex(omegaOfVpn,g1,g2,g3,g4,g5,g6,selfHeatGain,selfHeatGainComplex,domegabybo,domegabybTPA,c1,c2,Ts,aInitial);

so=sqrt(Pin)-1j*a*sqrt(gammaKappa);
end
