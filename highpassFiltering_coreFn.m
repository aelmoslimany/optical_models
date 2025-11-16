function [Vout] = highpassFiltering_coreFn(Vin,A,omegaPoleTs,ADC)
Vout1stOrder=zeros(length(Vin),length(A));
Vout=zeros(size(Vin));
omegaPoleInst=omegaPoleTs(:,1);
B1=exp(-omegaPoleInst.');
for i=3:length(Vin)
    Vout1stOrder(i,:)=B1.*Vout1stOrder(i-1,:)+(A.').*(1-B1).*Vin(i);
    Vout(i)=sum(Vout1stOrder(i,:));
end
Vout=ADC*Vin-Vout;
end
