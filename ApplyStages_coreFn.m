function Vout = ApplyStages_coreFn(AFEModel,Vin)
    Vout=Vin;
    for k = 1:AFEModel.NumStages

        if strcmp(AFEModel.StagesTypes{k}, 'LinearIR')
            % Apply linear Filtering using Impulse Response
            VoutTmp = conv(Vout,AFEModel.linearIR{k});
            Vout=VoutTmp(1:length(Vout));

        elseif strcmp(AFEModel.StagesTypes{k}, 'LinearRational')
          
            [Vout] = linearFiltering_coreFn(Vout,AFEModel.A{k},AFEModel.omegaPoleTs{k},AFEModel.AR{k},AFEModel.AI{k},AFEModel.alphaTs{k},AFEModel.betaTs{k});
        
        elseif strcmp(AFEModel.StagesTypes{k}, 'NonlinearODE')
            Input=Vout;
            [Vout] = dynamicNonlinearity_coreFn(Input,AFEModel.A{k},AFEModel.omegaPoleTs{k},AFEModel.omegaPoleTs1{k},AFEModel.AR{k},AFEModel.AI{k},AFEModel.alphaTs{k},AFEModel.betaTs{k},AFEModel.alphaTs1{k},AFEModel.betaTs1{k});
        elseif strcmp(AFEModel.StagesTypes{k}, 'LinearHPFRational')
      
            [Vout] = linearHPFFiltering_coreFn(Vout,AFEModel.AAllPass{k},AFEModel.A{k},AFEModel.omegaPoleTs{k},AFEModel.AR{k},AFEModel.AI{k},AFEModel.alphaTs{k},AFEModel.betaTs{k});
        
        elseif strcmp(AFEModel.StagesTypes{k}, 'NonlinearOddPoly')

            % Nonlinear Stage
            polyCoeff=AFEModel.NonlinearOddPoly{k};
            tmp = polyCoeff(1)*Vout;
            for kk=2:length(polyCoeff)
                tmp= tmp + polyCoeff(kk)*Vout.^(2*kk-1);
            end
            Vout=tmp;
        elseif strcmp(AFEModel.StagesTypes{k}, 'NonlinearOddNormPoly')

            % Nonlinear Stage
            polyCoeff=[1;AFEModel.NonlinearOddNormPoly{k}];
            tmp = polyCoeff(1)*Vout;
            for kk=2:length(polyCoeff)
                tmp= tmp + polyCoeff(kk)*Vout.^(2*kk-1);
            end
            Vout=tmp;
        elseif strcmp(AFEModel.StagesTypes{k}, 'NonlinearPoly')

            % Nonlinear Stage
            polyCoeff=AFEModel.NonlinearPoly{k};
            tmp = polyCoeff(1);
            for kk=2:length(polyCoeff)
                tmp= tmp + polyCoeff(kk)*Vout.^(kk-1);
            end
            Vout=tmp;   
        end
       
    end

end




