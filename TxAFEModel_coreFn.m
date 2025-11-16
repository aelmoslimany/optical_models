function [Vout,Gain,VoutLinear] = TxAFEModel_coreFn(Vin,AutoGain,Gain,nStages,LinearIRCell,NonlinearCoeffCell,clippingRatio,VppReq,Oversampling)
%% Initialization
% step=200/Oversampling;
step=240/Oversampling;
%% Auto Gain or not
if(AutoGain)
    % Initialization
    Vout=Vin;
    % Apply only Linear Stages
    for i=1:nStages
        % Linear Stage
        tmp=LinearIRCell{i};
        tmp=step*tmp(1:step:end);
        Vout=real(ifft(fft(Vout).*fft(tmp,length(Vout))));
    end
    P = prctile(Vout,clippingRatio/2);
    Gain=(-1/P)*(VppReq/2);
    VoutLinear=Gain*Vout;
end
%% Initialization
Vout=Gain*Vin;
%% Apply the stages
for i=1:nStages
    % Linear Stage
    tmp=LinearIRCell{i};
    tmp=step*tmp(1:step:end);
    Vout=real(ifft(fft(Vout).*fft(tmp,length(Vout))));
    % Nonlinear Stage
    polyCoeff=NonlinearCoeffCell{i};
    tmp = polyCoeff(1)*Vout;
    for k=2:length(polyCoeff)
        tmp= tmp + polyCoeff(k)*Vout.^(2*k-1);
    end
    Vout=tmp;
end
end
