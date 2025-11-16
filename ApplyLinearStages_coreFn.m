function Vout = ApplyLinearStages_coreFn(AFEModel,Vin)
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
            [Vout] = linearFiltering_coreFn(Input,AFEModel.A{k},AFEModel.omegaPoleTs{k},AFEModel.AR{k},AFEModel.AI{k},AFEModel.alphaTs{k},AFEModel.betaTs{k});
        elseif strcmp(AFEModel.StagesTypes{k}, 'LinearHPFRational')
      
            [Vout] = linearHPFFiltering_coreFn(Vout,AFEModel.AAllPass{k},AFEModel.A{k},AFEModel.omegaPoleTs{k},AFEModel.AR{k},AFEModel.AI{k},AFEModel.alphaTs{k},AFEModel.betaTs{k});
        
        elseif strcmp(AFEModel.StagesTypes{k}, 'NonlinearOddPoly')

            % Nonlinear Stage
            polyCoeff=AFEModel.NonlinearOddPoly{k};
            Vout = polyCoeff(1)*Vout;
            
        elseif strcmp(AFEModel.StagesTypes{k}, 'NonlinearOddNormPoly')

            % Nonlinear Stage
                   
        elseif strcmp(AFEModel.StagesTypes{k}, 'NonlinearPoly')

            % Nonlinear Stage
            polyCoeff=AFEModel.NonlinearPoly{k};
            Vout = polyCoeff(1)+polyCoeff(2)*Vout;
            
        end
       
    end

end




