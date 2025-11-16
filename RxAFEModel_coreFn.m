function [Vout,Gain,VoutLinear] = RxAFEModel_coreFn(Iin,AutoGain,Gain,nStages,LinearIRCell,NonlinearCoeffCell,clippingRatio,VppReq,A,omegaTs,AR,AI,alphaTs,betaTs,EnableTIANoise,Ts,Oversampling)
%% Initialization
step=200/Oversampling;
Iin=Iin-mean(Iin);
%% Auto Gain or not
if(AutoGain)
    % Initialization
    Vout=Iin;
    % % Apply only Linear Stages
    for i=1:nStages
        if(i==2)
            omegaTsTmp=omegaTs;
            omegaTsTmp(:,2)=0;
            alphaTsTmp=alphaTs;
            alphaTsTmp(:,2)=0;
            betaTsTmp=betaTs;
            betaTsTmp(:,2)=0;
            [Vout] = dynamicNonlinearity4_coreFn([Vout(end-1000+1:end);Vout],A,omegaTsTmp*step,AR,AI,alphaTsTmp*step,betaTsTmp*step);
            Vout=Vout(1000+1:end);
        else
            % Linear Stage
            tmp=LinearIRCell{i};
            tmp=step*tmp(1:step:end);
            Vout=real(ifft(fft(Vout).*fft(tmp,length(Vout))));
        end
    end
    P = prctile(Vout,clippingRatio/2);
    Gain=(-1/P)*VppReq/2;
    VoutLinear=Gain*Vout;
end
%% Initialization
if(EnableTIANoise)
    
GaindB=mag2db(Gain);
Gain2dB=72.5+GaindB;
noiseSigma=(14+20.5*(1+tanh(-0.957-0.22*(Gain2dB-56)+0.00724*(Gain2dB-56).^2)))*1e-12/sqrt(2*Ts);
Vout=Gain*Iin+Gain*noiseSigma*randn(size(Iin));
else
Vout=Gain*Iin;
end
%% Apply the stages
for i=1:nStages
    if(i==2)
        [Vout] = dynamicNonlinearity4_coreFn([Vout(end-1000+1:end);Vout],A,omegaTs*step,AR,AI,alphaTs*step,betaTs*step);
        Vout=Vout(1000+1:end);

    else
        % Linear Stage
        tmp=LinearIRCell{i};
        tmp=step*tmp(1:step:end);
        Vout=real(ifft(fft(Vout).*fft(tmp,length(Vout))));

    end
    % Nonlinear Stage
    polyCoeff=NonlinearCoeffCell{i};
    tmp = polyCoeff(1)*Vout;
    for k=2:length(polyCoeff)
        tmp= tmp + polyCoeff(k)*Vout.^(2*k-1);
    end
    Vout=tmp;

end
end
