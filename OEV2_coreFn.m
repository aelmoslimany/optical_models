function [Iout] = OEV2_coreFn(OpticalElectricFieldFastAxisIn,OpticalElectricFieldSlowAxisIn,Ts,Oversampling,OE)
OECouplingLossesFastAxis=db2pow(-1*OE.CouplingLossesFastAxis);
OECouplingLossesSlowAxis=db2pow(-1*OE.CouplingLossesSlowAxis);

switch OE.Type

    case 'Ideal'


        responsivityFastAxis = OE.Responsivity * OECouplingLossesFastAxis;
        responsivitySlowAxis = OE.Responsivity * OECouplingLossesSlowAxis;
        
        Iout=responsivityFastAxis*abs(OpticalElectricFieldFastAxisIn).^2+responsivitySlowAxis*abs(OpticalElectricFieldSlowAxisIn).^2;
    
        if strcmp(OE.EnableShotNoise,'true')
            Iavg = mean(Iout);
    
            % Implementation with shot noise and dark noise
            sigmaNoise = sqrt(2 * 1.602176634e-19 * abs(Iavg) * (1 / (2 * Ts)));
            Iout = Iout + sigmaNoise * randn(size(Iout));
        end

    case 'GFPD'
        opticalPower=OECouplingLossesFastAxis*abs(OpticalElectricFieldFastAxisIn).^2+OECouplingLossesSlowAxis*abs(OpticalElectricFieldSlowAxisIn).^2;
        opticalPowerTmp=opticalPower;
        opticalPowerTmp(pow2db(opticalPowerTmp)<-40) = db2pow(-40);
        Responsivity = Detector_responsivity(opticalPowerTmp,OE.DCBias,OE.PDWidth,OE.Wavelength*1e-9); %%Detector resposivity vs time
        C = TransitBandwidth(opticalPowerTmp,OE.DCBias,OE.PDWidth,OE.Temperature); %%Transit time Equivalent capacitor
        % Differential equation parameters
        N=length(opticalPower);
        timeVec=Ts*(0:N-1).';
        Q = 1.602176634e-19 ;

        if strcmp(OE.EnableShotNoise,'true')
        f = (Responsivity.*opticalPower + sqrt(2*Q*(Responsivity.*opticalPower)*(1/(2*Ts))).*randn(size(opticalPower)))./(2*C);%% Main signal + Shot noise
        else
        f = (Responsivity.*opticalPower)./(2*C);
        end
        
                   
        a = (1./(2*C));
        
        Id2 = zeros(N,1);
        Id2(1) = (Ts*(f(1))/2)./(1+a(1)*Ts/2);
        % tic;
        for i=2:N
        Id2(i)=(Id2(i-1)-(a(i-1)*Id2(i-1))*Ts/2+Ts*(f(i)+f(i-1))/2)./(1+a(i)*Ts/2);
        end
        
        
        if strcmp(OE.EnableShotNoise,'true')
        Id_dark = -1*DarkCurrent(OE.DCBias,OE.PDWidth,OE.Temperature);
        Id_dark = Id_dark + sqrt(2*Q*abs(Id_dark)*(1/(2*Ts))).*randn(size(Id2));%%Add dark current shot noise
        Id2 = Id2 + Id_dark;
        end
        
        projectRoot = getenv('OPTICALMODELPROJROOT');
        tmp=readmatrix(fullfile(projectRoot,OE.ZinIRFile));
        n1stOrder=tmp(1);
        k=2;
        OE.A=tmp(k:k+n1stOrder-1);
        k=k+n1stOrder;
        OE.omegaPoleTs=tmp(k:k+n1stOrder-1)*Ts;
        k=k+n1stOrder;
        n2ndOrder=tmp(k);
        k=k+1;
        OE.AR=tmp(k:k+n2ndOrder-1);
        k=k+n2ndOrder;
        OE.AI=tmp(k:k+n2ndOrder-1);
        k=k+n2ndOrder;
        OE.alphaTs=tmp(k:k+n2ndOrder-1)*Ts;
        k=k+n2ndOrder;
        OE.betaTs=tmp(k:k+n2ndOrder-1)*Ts;

        [Iout] = linearFiltering_coreFn(Id2,OE.A,OE.omegaPoleTs,OE.AR,OE.AI,OE.alphaTs,OE.betaTs);
        % Iout=Id2;
end

end
