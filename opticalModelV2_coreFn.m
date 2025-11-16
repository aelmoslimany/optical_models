function [Vout,TECQ,TDECQ,OMA,ER,Ceq,RLM,EECQ,CeqoutRx,RLMRx,SNRdB,BER] = opticalModelV2_coreFn(Vin,Ts,Oversampling,OpticalModelInit,SSPRQseqSymUndersampled,EnableTDECQCalc,TDECQNumFold,EnableEyeDiagramPlot,EnableEECQCalc,VRMS,BERType)
%% Tx AFE
[Vout] = AFEModelV2_coreFn(Vin,Ts,Oversampling,OpticalModelInit.TXAFE);

%% E/O
[OpticalElectricFieldOut] = EOV2_coreFn(Vout,Ts,Oversampling,OpticalModelInit.EO);

%% Fiber
[OpticalElectricFieldFastAxisOut,OpticalElectricFieldSlowAxisOut] = FIBERModelV2_coreFn(OpticalElectricFieldOut,Ts,Oversampling,OpticalModelInit.FIBER);

%% TDECQ Calc
if(EnableTDECQCalc)
    if(TDECQNumFold>=3)
        startInd=length(SSPRQseqSymUndersampled)/TDECQNumFold +1;
        endInd=length(SSPRQseqSymUndersampled)*(TDECQNumFold-1)/TDECQNumFold;
        oversamplingStartInd=length(OpticalElectricFieldFastAxisOut)/TDECQNumFold +1;
        oversamplingEndInd=length(OpticalElectricFieldFastAxisOut)*(TDECQNumFold-1)/TDECQNumFold;
        [TECQ,TDECQ,OMA,ER,Ceq,RLM] = TDECQCalc_coreFn(abs(OpticalElectricFieldOut(oversamplingStartInd:oversamplingEndInd)).^2,abs(OpticalElectricFieldFastAxisOut(oversamplingStartInd:oversamplingEndInd)).^2+abs(OpticalElectricFieldSlowAxisOut(oversamplingStartInd:oversamplingEndInd)).^2,SSPRQseqSymUndersampled(startInd:endInd),Oversampling,EnableEyeDiagramPlot);
    else
        [TECQ,TDECQ,OMA,ER,Ceq,RLM] = TDECQCalc_coreFn(abs(OpticalElectricFieldOut).^2,abs(OpticalElectricFieldFastAxisOut).^2+abs(OpticalElectricFieldSlowAxisOut).^2,SSPRQseqSymUndersampled,Oversampling,EnableEyeDiagramPlot);

    end
else
    TECQ=nan;
    TDECQ=nan;
    OMA=nan;
    ER=nan;
    Ceq=nan;
    RLM=nan;
end

%% O/E
[Iout] = OEV2_coreFn(OpticalElectricFieldFastAxisOut,OpticalElectricFieldSlowAxisOut,Ts,Oversampling,OpticalModelInit.OE);

%% Rx AFE
[Vout] = AFEModelV2_coreFn(Iout,Ts,Oversampling,OpticalModelInit.RXAFE);

%% EECQ Calc
if(EnableEECQCalc)
    if(TDECQNumFold>=3)
        startInd=length(SSPRQseqSymUndersampled)/TDECQNumFold +1;
        endInd=length(SSPRQseqSymUndersampled)*(TDECQNumFold-1)/TDECQNumFold;
        oversamplingStartInd=length(OpticalElectricFieldFastAxisOut)/TDECQNumFold +1;
        oversamplingEndInd=length(OpticalElectricFieldFastAxisOut)*(TDECQNumFold-1)/TDECQNumFold;
        [EECQ,CeqoutRx,RLMRx,SNRdB,BER] = EECQCalc_coreFn(Vout(oversamplingStartInd:oversamplingEndInd),SSPRQseqSymUndersampled(startInd:endInd),Oversampling,VRMS,BERType);
    else
        [EECQ,CeqoutRx,RLMRx,SNRdB,BER] = EECQCalc_coreFn(Vout,SSPRQseqSymUndersampled,Oversampling,VRMS,BERType);

    end
else
    EECQ=nan;
    CeqoutRx=nan;
    RLMRx=nan;
    SNRdB=nan;
    BER=nan;
end

end
