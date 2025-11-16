function [PLaserCoupled,RinNoise] = LaserModel_coreFn(Ts,Len,Power,CouplingLosses,RIN)
PLaserCoupled=Power-CouplingLosses;
NoisePower=db2pow(PLaserCoupled+RIN-30)/(2*Ts);
PLaserCoupled=db2pow(PLaserCoupled-30);
RinNoise=0*sqrt(NoisePower)*randn(Len,1);
end
