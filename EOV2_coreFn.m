function [OpticalElectricFieldOut] = EOV2_coreFn(Vin,Ts,Oversampling,EO)
projectRoot = getenv('OPTICALMODELPROJROOT');


switch EO.Type

    case 'Ideal'
        if strcmp(EO.EnableAutoCalcVppIn,'true')
            Vp = prctile(Vin, [100*EO.ClippingRatio/2, 100*(1-EO.ClippingRatio/2)]);
            EO.VppIn = 2*max(Vp);
        end
        [Pout] = powerModulator_coreFn(Vin,EO.VppIn,db2pow(EO.ER),db2pow(EO.OMA-30));
        Pout(Pout<0)=0;
        OpticalElectricFieldOut=sqrt(db2pow(-1*EO.EOCouplingLosses)*Pout);

    case 'EML'
        EO.LaserPower=db2pow(EO.LaserPower-30)*db2pow(-1*EO.LaserInCouplingLosses);
        for i=1:EO.NumStages
            tag=['Stage' num2str(i) 'ParamsFile'];
            filePath=fullfile(projectRoot,EO.(tag));
            params=readmatrix(filePath);
            switch EO.StagesTypes{i}

                case 'LinearRational'
                    n1stOrder=params(1);
                    k=2;
                    EO.A=params(k:k+n1stOrder-1);
                    k=k+n1stOrder;
                    EO.omegaPoleTs=params(k:k+n1stOrder-1)*Ts;
                    k=k+n1stOrder;
                    n2ndOrder=params(k);
                    k=k+1;
                    EO.AR=params(k:k+n2ndOrder-1);
                    k=k+n2ndOrder;
                    EO.AI=params(k:k+n2ndOrder-1);
                    k=k+n2ndOrder;
                    EO.alphaTs=params(k:k+n2ndOrder-1)*Ts;
                    k=k+n2ndOrder;
                    EO.betaTs=params(k:k+n2ndOrder-1)*Ts;

                    [Vin] = linearFiltering_coreFn(Vin,EO.A,EO.omegaPoleTs,EO.AR,EO.AI,EO.alphaTs,EO.betaTs);
                case 'LinearIR'
                    % Linear Stage
                    step=Ts/params(1);
                    LinearIR=params(2:end);
                    tmp=interp1((0:length(LinearIR)-1)*params(1),LinearIR,0:Ts:(length(LinearIR)-1)*params(1),"linear");
                    tmp=step*tmp.';
                    VinTmp=[Vin;zeros(length(tmp),1)];
                    VinTmp=real(ifft(fft(VinTmp).*fft(tmp,length(VinTmp))));
                    Vin=VinTmp(1:length(Vin));

                case 'NonlinearPoly'
                    % Add Bias
                    Vin=Vin+EO.Vbias;
                    % Nonlinear Stage
                    polyCoeff=params;
                    tmp = polyCoeff(2)*Vin+polyCoeff(1);
                    for k=3:length(polyCoeff)
                        tmp= tmp + polyCoeff(k)*Vin.^(k-1);
                    end
                    Vin=0.5*(tanh(tmp)+1);
                    OpticalElectricFieldOut=sqrt(EO.LaserPower*db2pow(-EO.EOCouplingLosses).*Vin);

                case 'NonlinearPolyDynamicTemperature'
                    n1stOrder=params(1);
                    k=2;
                    EO.A=params(k:k+n1stOrder-1);
                    k=k+n1stOrder;
                    EO.omegaPoleTs=params(k:k+n1stOrder-1)*Ts;
                    k=k+n1stOrder;
                    n2ndOrder=params(k);
                    k=k+1;
                    EO.AR=params(k:k+n2ndOrder-1);
                    k=k+n2ndOrder;
                    EO.AI=params(k:k+n2ndOrder-1);
                    k=k+n2ndOrder;
                    EO.alphaTs=params(k:k+n2ndOrder-1)*Ts;
                    k=k+n2ndOrder;
                    EO.betaTs=params(k:k+n2ndOrder-1)*Ts;
                    k=k+n2ndOrder;
                    EO.Poly0=params(k:3:end);
                    EO.Poly1=params(k+1:3:end);
                    EO.Poly2=params(k+2:3:end);
                    [Pout] = polyWithHighPowerNonlinearity_coreFn(Vin,EO.A,EO.omegaPoleTs,EO.AR,EO.AI,EO.alphaTs,EO.betaTs,EO.Poly0,EO.Poly1,EO.Poly2,EO.Vbias,EO.LaserPower);
                    Pout(Pout<0)=0;
                    OpticalElectricFieldOut=sqrt(db2pow(-1*EO.EOCouplingLosses)*Pout);

                case 'NonlinearODE'
                    n1stOrder=params(1);
                    k=2;
                    EO.A=params(k:k+n1stOrder-1);
                    k=k+n1stOrder;
                    EO.omegaPoleTs=params(k:k+n1stOrder-1)*Ts;
                    k=k+n1stOrder;
                    EO.omegaPoleTs1=params(k:k+n1stOrder-1)*Ts;
                    k=k+n1stOrder;
                    n2ndOrder=params(k);
                    k=k+1;
                    EO.AR=params(k:k+n2ndOrder-1);
                    k=k+n2ndOrder;
                    EO.AI=params(k:k+n2ndOrder-1);
                    k=k+n2ndOrder;
                    EO.alphaTs=params(k:k+n2ndOrder-1)*Ts;
                    k=k+n2ndOrder;
                    EO.betaTs=params(k:k+n2ndOrder-1)*Ts;
                    k=k+n2ndOrder;
                    EO.alphaTs1=params(k:k+n2ndOrder-1)*Ts;
                    k=k+n2ndOrder;
                    EO.betaTs1=params(k:k+n2ndOrder-1)*Ts;
                    [Vin] = dynamicNonlinearity_coreFn(Vin,EO.A,EO.omegaPoleTs,EO.omegaPoleTs1,EO.AR,EO.AI,EO.alphaTs,EO.betaTs,EO.alphaTs1,EO.betaTs1);

            end
        end

    case 'MZM'
        if strcmp(EO.EnableAutoCalcVppIn,'true')
            Vp = prctile(Vin, [100*EO.ClippingRatio/2, 100*(1-EO.ClippingRatio/2)]);
            EO.VppIn = EO.VppIn./(2*max(Vp));
            Vin=Vin*EO.VppIn;
        end
        EO.LaserPower=db2pow(EO.LaserPower-30)*db2pow(-1*EO.LaserInCouplingLosses);
        Vin2=[Vin;zeros(1e4,1)];
        EO.Wavelength=EO.Wavelength*1e-9;
        EO.LaserPower=EO.LaserPower*db2pow(EO.EOCouplingLosses);
        [~,OpticalElectricFieldOut] = MZMmodelfn(Vin2,Ts,EO.Vbias,EO.MZMLength,EO.NumOfSections,EO.MZMSourceImpedance,EO.MZMLoadImpedance,EO.MZMArmMismatch,EO.LaserPower,EO.MZMHeaterPhase,EO.Wavelength);
        OpticalElectricFieldOut=OpticalElectricFieldOut(:);
        lenPre = round((50e-12 / Ts) * (EO.MZMLength/5e-3));
        OpticalElectricFieldOut=circshift(OpticalElectricFieldOut,lenPre);
        OpticalElectricFieldOut=OpticalElectricFieldOut(1:length(Vin));


    case 'MRM'
        EO.LaserPower=db2pow(EO.LaserPower-30)*db2pow(-1*EO.LaserInCouplingLosses);
        filePath=fullfile(projectRoot,EO.('TempFilterParamsFile'));
        params=readmatrix(filePath);
        n1stOrder=params(1);
        k=2;
        EO.A=params(k:k+n1stOrder-1);
        k=k+n1stOrder;
        EO.omegaPoleTs=params(k:k+n1stOrder-1)*Ts;
        k=k+n1stOrder;
        n2ndOrder=params(k);
        k=k+1;
        EO.AR=params(k:k+n2ndOrder-1);
        k=k+n2ndOrder;
        EO.AI=params(k:k+n2ndOrder-1);
        k=k+n2ndOrder;
        EO.alphaTs=params(k:k+n2ndOrder-1)*Ts;
        k=k+n2ndOrder;
        EO.betaTs=params(k:k+n2ndOrder-1)*Ts;
        selfHeat1stOrderGains=EO.A.';
        selfHeat1stOrderPolesHz=EO.omegaPoleTs.'/(2*pi*Ts);
        selfHeat2ndOrderPolesHz=(EO.alphaTs+1j*EO.betaTs).'/(2*pi*Ts);
        selfHeat2bdOrderGains=(EO.AR+1j* EO.AI).';
        RoundTripTime = 2 * pi * EO.MRMRadius * EO.ng / 3e8;
        gammaKappa = (EO.KappaSq / 100) / RoundTripTime;
        gammaLoss = (1 - 10.0^(EO.LossPerRoundTripdB / 10.0)) / RoundTripTime;
        gammaTot = gammaKappa + gammaLoss;
        Vin = sum(EO.OmegaofVpnCoeff .* (Vin+EO.Vbias).^(1:length(EO.OmegaofVpnCoeff)),2);
        [OpticalElectricFieldOut,~,~] = ringDynamicsWithHighPowerNonlinearityWithMEX_coreFn(Vin,Ts,2*pi*80e9,1,gammaTot,gammaKappa,selfHeat1stOrderPolesHz,selfHeat1stOrderGains,selfHeat2ndOrderPolesHz,selfHeat2bdOrderGains,EO.domegaobydAbsaSq,EO.domegaobydAbsaSqSq,EO.LaserPower,db2pow(-1*EO.ILTarget));
end



if strcmp(EO.EnableRinNoise,'true')
    % parse rational filtering file
    tmp=readmatrix(fullfile(projectRoot,EO.RinFilterParamsFile));
    n1stOrder=tmp(1);
    k=2;
    EO.A=tmp(k:k+n1stOrder-1);
    k=k+n1stOrder;
    EO.omegaPoleTs=tmp(k:k+n1stOrder-1)*Ts;
    k=k+n1stOrder;
    n2ndOrder=tmp(k);
    k=k+1;
    EO.AR=tmp(k:k+n2ndOrder-1);
    k=k+n2ndOrder;
    EO.AI=tmp(k:k+n2ndOrder-1);
    k=k+n2ndOrder;
    EO.alphaTs=tmp(k:k+n2ndOrder-1)*Ts;
    k=k+n2ndOrder;
    EO.betaTs=tmp(k:k+n2ndOrder-1)*Ts;

    RinNoiseUnfiltered=randn(size(Vin))/sqrt(2*Ts);
    [RinNoiseFiltered] = linearFiltering_coreFn(RinNoiseUnfiltered,EO.A,EO.omegaPoleTs,EO.AR,EO.AI,EO.alphaTs,EO.betaTs);
    tmp=1+RinNoiseFiltered;
    tmp(tmp<0)=0;
    OpticalElectricFieldOut=sqrt(tmp).*OpticalElectricFieldOut;
end


end



