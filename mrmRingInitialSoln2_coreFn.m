function [omegaBiasInitial, aSq, aInitial, isFailed] = mrmRingInitialSoln2_coreFn(domega_obydV, inputVoltageLevels, gammaKappa, Pin, gammaTot, omegaBiasInitial, domegabybo, domegabybTPA, IL)
%% GET Initial Solution
% inputVoltageMean=domega_obydV*mean(inputVoltageLevels);
aOld1=(1j*sqrt(gammaKappa)*sqrt(Pin)./(-gammaTot/2 +1j*(omegaBiasInitial+domega_obydV*inputVoltageLevels)));
aOld=mean(aOld1);
aSq=0;
aSqOld=0;
aSqSq=0;
aSqSqOld=0;
TRegulaization=1e-6;
Omegap=1/1e-4;
mm=inf;
OmegaStep=2*pi*1e9;
wasItPos=true;
iLoop1=0;
while(mm>1e-5 && iLoop1<1e3)
    iLoop1=iLoop1+1;
    m=inf;
    iLoop=0;
    while((~(m<1e-10 && iLoop>20)) && iLoop<1e3)
        aInitial1=1j*sqrt(gammaKappa)*sqrt(Pin)./(-gammaTot/2 +1j*(omegaBiasInitial+domega_obydV*inputVoltageLevels+domegabybo*aSq+ domegabybTPA*aSqSq ));
        aInitial=mean(aInitial1);
        aSq=aSqOld+TRegulaization*Omegap*(-aSqOld+mean(abs(aOld1).^2));
        aSqOld=aSq;
        aSqSq=aSqSqOld+TRegulaization*Omegap*(-aSqSqOld+mean(abs(aOld1).^4));
        aSqSqOld=aSqSq;
        m=abs(aInitial-aOld);
        aOld1=aInitial;
        iLoop=iLoop+1;
    end
    mm=(mean(abs(sqrt(Pin)-1j*sqrt(gammaKappa)*aInitial1).^2)-IL*Pin);
    % (abs(sqrt(Pin)-1j*sqrt(gammaKappa)*aInitial).^2/Pin)
    if(mm>0)
        if(~wasItPos)
            OmegaStep=OmegaStep/2;
        end
        omegaBiasInitial=omegaBiasInitial-OmegaStep;
        wasItPos=true;
    elseif(mm<0)
        if(wasItPos)
            OmegaStep=OmegaStep/2;
        end
        omegaBiasInitial=omegaBiasInitial+OmegaStep;
        wasItPos=false;
    end
    mm=abs(mm);
end
if(mm>1e-5)
    isFailed=true;
else
    isFailed=false;
end
end
