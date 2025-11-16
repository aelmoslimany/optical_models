function [so,a,OmegaTotOfVpn] = ringDynamicsWithHighPowerNonlinearityWithoutMEX_coreFn(inputVoltage,Ts,omegaBiasInitial,domega_obydV,gammaTot,gammaKappa,selfHeat1stOrderPolesHz,selfHeat1stOrderGains,selfHeat2ndOrderPolesHz,selfHeat2bdOrderGains,domegabybo,domegabybTPA,Pin,IL)
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

%% MRM Without Mex
omegaOfVpnLen=length(inputVoltage);
LPFLen=length(selfHeatGain);
complexLPFLen=length(selfHeatGainComplex);
OmegaTotOfVpn=zeros(omegaOfVpnLen,1);
a=zeros(omegaOfVpnLen,1);
ainit=aInitial;
bo = abs(aInitial).^2;
bTPA=abs(bo).^2;
aSq=bo;
aSqSq=bTPA;
aSqD=bo;
aSqSqD=bTPA;
y0 = zeros(LPFLen,1);
ySq0 = zeros(LPFLen,1);
for i_tmp2=1:LPFLen
    y0(i_tmp2) = bo*selfHeatGain(i_tmp2);
    ySq0(i_tmp2) = bTPA*selfHeatGain(i_tmp2);
end

z0 = zeros(complexLPFLen,1);
z1 = zeros(complexLPFLen,1);
zSq0 = zeros(complexLPFLen,1);
zSq1 = zeros(complexLPFLen,1);
for i_tmp2=1:complexLPFLen
    z0(i_tmp2) = bo*selfHeatGainComplex(i_tmp2);
    z1(i_tmp2) = bo*selfHeatGainComplex(i_tmp2);
    zSq0(i_tmp2) = bTPA*selfHeatGainComplex(i_tmp2);
    zSq1(i_tmp2) = bTPA*selfHeatGainComplex(i_tmp2);
end


for i_tmp=1:omegaOfVpnLen
    OmegaTotOfVpn(i_tmp)=omegaOfVpn(i_tmp)+domegabybo*bo+domegabybTPA*bTPA;
    f0=c1+1i*OmegaTotOfVpn(i_tmp);
    f1=exp(Ts*f0);
    f2=(-1i*c2)/f0;
    a(i_tmp) = f1*ainit + f2*(f1-1.0);
    ainit=a(i_tmp);

    % 1st order LPF
    for i_tmp2=1:LPFLen
        y0(i_tmp2) = g1(i_tmp2)*y0(i_tmp2) + g2(i_tmp2)*aSq;
        ySq0(i_tmp2) = g1(i_tmp2)*ySq0(i_tmp2) + g2(i_tmp2)*aSqSq;
    end

    % 2nd order LPF
    for i_tmp2=1:complexLPFLen
        tmp = z0(i_tmp2);
        z0(i_tmp2) = g3(i_tmp2)*z0(i_tmp2) + g4(i_tmp2)*z1(i_tmp2) + g5(i_tmp2)*aSq + g6(i_tmp2)*aSqD;
        z1(i_tmp2) = tmp;
        tmp = zSq0(i_tmp2);
        zSq0(i_tmp2) = g3(i_tmp2)*zSq0(i_tmp2) + g4(i_tmp2)*zSq1(i_tmp2) + g5(i_tmp2)*aSqSq + g6(i_tmp2)*aSqSqD;
        zSq1(i_tmp2) = tmp;
    end

    % Compute the output of the LPF
    aSqD=aSq;
    aSqSqD=aSqSq;
    aSq=abs(a(i_tmp)).^2;
    aSqSq=abs(aSq).^2;
    bo=0;
    bTPA=0;

    for i_tmp2=1:LPFLen
        bo = bo + y0(i_tmp2);
        bTPA = bTPA + ySq0(i_tmp2);
    end

    for i_tmp2=1:complexLPFLen
        bo = bo + z0(i_tmp2);
        bTPA = bTPA + zSq0(i_tmp2);
    end
end

%% Output Field
so=sqrt(Pin)-1j*a*sqrt(gammaKappa);
end
