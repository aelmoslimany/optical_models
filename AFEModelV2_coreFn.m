function [Vout] = AFEModelV2_coreFn(Vin,Ts,Oversampling,AFE)

switch AFE.Type

    case 'Ideal'

        % DC Cancellation implementation
        if strcmp((AFE.EnableInputDCCancellation),'true')
            VinDC = mean(Vin(1:end));
            Vin = Vin - VinDC;
        end
    
        % AGC implementation
        if strcmp((AFE.EnableAGC),'true')
            Vp = prctile(Vin, [100*AFE.AGCClippingRatio/2, 100*(1-AFE.AGCClippingRatio/2)]);
            AFE.Gain = AFE.VppReq*0.5./max(abs(Vp));
            AFE.Gain = min([AFE.Gain, AFE.GainMax]);
            AFE.Gain = max([AFE.Gain, AFE.GainMin]);
        end
    
        % Compute the standard deviation of the input referred noise
        if strcmp((AFE.EnableInputRefNoise),'true')
            GaindB = 20 * log10(AFE.Gain);
            inputRefNoiseSigma = (AFE.NoiseGainParams(1) + AFE.NoiseGainParams(2) * (tanh(AFE.NoiseGainParams(3) + AFE.NoiseGainParams(4) * GaindB + AFE.NoiseGainParams(5) * GaindB^2))) / sqrt(2 * Ts);
            Vin = Vin + inputRefNoiseSigma * randn(size(Vin));
        end
    
        % Apply the gain
        Vout = Vin * AFE.Gain;


    % Get the properties for "GeneralizedHammersteinWienerModel" model from the configuration file
    case 'GeneralizedHammersteinWienerModel'
        % DC Cancellation implementation
        projectRoot = getenv('OPTICALMODELPROJROOT');
        if strcmp((AFE.EnableInputDCCancellation),'true')
            VinDC = mean(Vin(1:end));
            Vin = Vin - VinDC;
        end

        if AFE.NumStages > 0
            if(~iscell(AFE.StagesTypes))
                AFE.StagesTypes={AFE.StagesTypes};

            end
            
            for i = 1:AFE.NumStages
                if strcmp(AFE.StagesTypes{i}, 'LinearIR')
                    % Read the impulse response and resampling it
                    paramsFileField = sprintf('Stage%dParamsFile', i);
                    InputIR=readmatrix(fullfile(projectRoot,AFE.(paramsFileField)));
                    IRTs=InputIR(1);%first is the value of original Ts
                    InputIR=InputIR(2:end);
                    InputIRSize = length(InputIR);
                    OutputIRSize = floor(InputIRSize * IRTs / Ts);
                    % t=Ts:Ts:OutputIRSize*Ts;
                    % t_orig=0:IRTs:InputIRSize*IRTs-IRTs;
                    step=Ts/IRTs;
                    TsResampling=Ts;
                    IROut=step*InputIR(1);
                    % AFE.linearIR{i} = [(step) * interp1(t_orig,InputIR, t, 'linear')]; % Linear interpolation
                    for ii=1:InputIRSize-1 
                        while(TsResampling < ii * IRTs) 
                          slope = (InputIR(ii+1) - InputIR(ii)) / (IRTs);
                          doubleValue = slope * (TsResampling - (ii - 1) * IRTs) + InputIR(ii);
                          IROut=[IROut step * doubleValue]; 
                          TsResampling = TsResampling + Ts;

                        end

                    end

                    AFE.linearIR{i}=IROut;
                elseif strcmp(AFE.StagesTypes{i}, 'LinearRational')
                    paramsFileField = sprintf('Stage%dParamsFile', i);
                    tmp=readmatrix(fullfile(projectRoot,AFE.(paramsFileField)));
                    n1stOrder=tmp(1);
                    k=2;
                    AFE.A{i}=tmp(k:k+n1stOrder-1);
                    k=k+n1stOrder;
                    AFE.omegaPoleTs{i}=tmp(k:k+n1stOrder-1)*Ts;
                    k=k+n1stOrder;
                    n2ndOrder=tmp(k);
                    k=k+1;
                    AFE.AR{i}=tmp(k:k+n2ndOrder-1);
                    k=k+n2ndOrder;
                    AFE.AI{i}=tmp(k:k+n2ndOrder-1);
                    k=k+n2ndOrder;
                    AFE.alphaTs{i}=tmp(k:k+n2ndOrder-1)*Ts;
                    k=k+n2ndOrder;
                    AFE.betaTs{i}=tmp(k:k+n2ndOrder-1)*Ts;
                 elseif strcmp(AFE.StagesTypes{i}, 'NonlinearODE')

                    paramsFileField = sprintf('Stage%dParamsFile', i);
                    tmp=readmatrix(fullfile(projectRoot,AFE.(paramsFileField)));
                    n1stOrder=tmp(1);
                    k=2;
                    AFE.A{i}=tmp(k:k+n1stOrder-1);
                    k=k+n1stOrder;
                    AFE.omegaPoleTs{i}=tmp(k:k+n1stOrder-1)*Ts;
                    k=k+n1stOrder;
                    AFE.omegaPoleTs1{i}=tmp(k:k+n1stOrder-1)*Ts;
                    k=k+n1stOrder;
                    n2ndOrder=tmp(k);
                    k=k+1;
                    AFE.AR{i}=tmp(k:k+n2ndOrder-1);
                    k=k+n2ndOrder;
                    AFE.AI{i}=tmp(k:k+n2ndOrder-1);
                    k=k+n2ndOrder;
                    AFE.alphaTs{i}=tmp(k:k+n2ndOrder-1)*Ts;
                    k=k+n2ndOrder;
                    AFE.betaTs{i}=tmp(k:k+n2ndOrder-1)*Ts;
                    k=k+n2ndOrder;
                    AFE.alphaTs1{i}=tmp(k:k+n2ndOrder-1)*Ts;
                    k=k+n2ndOrder;
                    AFE.betaTs1{i}=tmp(k:k+n2ndOrder-1)*Ts;
                elseif strcmp(AFE.StagesTypes{i}, 'LinearHPFRational')
                    paramsFileField = sprintf('Stage%dParamsFile', i);
                    tmp=readmatrix(fullfile(projectRoot,AFE.(paramsFileField)));
                    AFE.AAllPass{i}=tmp(1);
                    n1stOrder=tmp(2);
                    k=3;
                    AFE.A{i}=tmp(k:k+n1stOrder-1);
                    k=k+n1stOrder;
                    AFE.omegaPoleTs{i}=tmp(k:k+n1stOrder-1)*Ts;
                    k=k+n1stOrder;
                    n2ndOrder=tmp(k);
                    k=k+1;
                    AFE.AR{i}=tmp(k:k+n2ndOrder-1);
                    k=k+n2ndOrder;
                    AFE.AI{i}=tmp(k:k+n2ndOrder-1);
                    k=k+n2ndOrder;
                    AFE.alphaTs{i}=tmp(k:k+n2ndOrder-1)*Ts;
                    k=k+n2ndOrder;
                    AFE.betaTs{i}=tmp(k:k+n2ndOrder-1)*Ts;
                elseif strcmp(AFE.StagesTypes{i}, 'NonlinearPoly') 
                    paramsFileField = sprintf('Stage%dParamsFile', i);
                    tmp=readmatrix(fullfile(projectRoot,AFE.(paramsFileField)));
                    AFE.NonlinearPoly{i}=tmp;
                elseif strcmp(AFE.StagesTypes{i}, 'NonlinearOddPoly') 
                    paramsFileField = sprintf('Stage%dParamsFile', i);
                    tmp=readmatrix(fullfile(projectRoot,AFE.(paramsFileField)));
                    AFE.NonlinearOddPoly{i}=tmp;  
                elseif strcmp(AFE.StagesTypes{i}, 'NonlinearOddNormPoly') 
                    paramsFileField = sprintf('Stage%dParamsFile', i);
                    tmp=readmatrix(fullfile(projectRoot,AFE.(paramsFileField)));
                    AFE.NonlinearOddNormPoly{i}=tmp;      
                end
            end

            
                if strcmp((AFE.EnableAGC),'true')
                    [Vout]=ApplyLinearStages_coreFn(AFE, Vin);
                    Vp = prctile(Vout, [100*AFE.AGCClippingRatio/2, 100*(1-AFE.AGCClippingRatio/2)]);
                    AFE.Gain = AFE.VppReq*0.5./max(abs(Vp));
                    AFE.Gain = min([AFE.Gain, AFE.GainMax]);
                    AFE.Gain = max([AFE.Gain, AFE.GainMin]);
                end
            
                % Compute the standard deviation of the input referred noise
                if strcmp((AFE.EnableInputRefNoise),'true')
                    GaindB = 20 * log10(AFE.Gain);
                    inputRefNoiseSigma = (AFE.NoiseGainParams(1) + AFE.NoiseGainParams(2) * (tanh(AFE.NoiseGainParams(3) + AFE.NoiseGainParams(4) * GaindB + AFE.NoiseGainParams(5) * GaindB^2))) / sqrt(2 * Ts);
                    Vin = Vin + inputRefNoiseSigma * randn(size(Vin));
                end


                 % Apply the gain and and apply AFE
                [Vout]=ApplyStages_coreFn(AFE, Vin * AFE.Gain);

            

    
         end


end


