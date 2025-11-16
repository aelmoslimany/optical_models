function [OpticalElectricFieldFastAxisOut,OpticalElectricFieldSlowAxisOut] = FIBERModelV2_coreFn(OpticalElectricFieldIn,Ts,Oversampling,FIBER)

switch FIBER.Type

    case 'Ideal'
    Loss=db2pow(-1*FIBER.Losses);

    % Initialization
    fastAccessFactor = sqrt(FIBER.FastAxisRatio * Loss);
    slowAccessFactor = sqrt((1 - FIBER.FastAxisRatio) * Loss);
    
    % Preallocate output arrays
    OpticalElectricFieldFastAxisOut = fastAccessFactor * OpticalElectricFieldIn;
    OpticalElectricFieldSlowAxisOut = slowAccessFactor * OpticalElectricFieldIn;

    case 'LinearFiber'
    
    % Frequency and Time Axis
    OpticalElectricFieldIn=[OpticalElectricFieldIn;sqrt(mean(abs(OpticalElectricFieldIn).^2))*ones(4000,1)];
    N=length(OpticalElectricFieldIn);
    freqAxis=(-1/(2*Ts):1/(Ts*N):1/(2*Ts)-1/(Ts*N)).';
    % Transform the input optical electric field to the frequency domain
    opticalInElectricFieldFreq=fftshift(fft(OpticalElectricFieldIn));
    % Chromatic Dispersion and IL
    IntrinsicDelay = (50e-12) * (FIBER.DispersionSlope * FIBER.Length) / (2 * 0.092) + FIBER.DGD*1e-12;
    beta2=-(FIBER.DispersionSlope*1e-12*FIBER.Wavelength.^3*1e-9*(1-(FIBER.ZeroDispersionWavelength./FIBER.Wavelength).^4)/(8*pi*3e8));
    CD=exp(-1j*beta2*FIBER.Length*((2*pi*freqAxis).^2)./2);
    opticalInElectricFieldAfterCDFreq=opticalInElectricFieldFreq.*CD*db2mag(-FIBER.Losses);
    % Differential Group Dleay 
    opticalOutElectricFieldSigFastFreq=sqrt(FIBER.FastAxisRatio)*opticalInElectricFieldAfterCDFreq.*exp(-1j*2*pi*freqAxis*IntrinsicDelay);
    opticalOutElectricFieldSigSlowFreq=sqrt(1-FIBER.FastAxisRatio)*opticalInElectricFieldAfterCDFreq.*exp(-1j*2*pi*freqAxis*(IntrinsicDelay+FIBER.DGD*1e-12));
    % Going Back to Time Domain 
    OpticalElectricFieldFastAxisOut=ifft(ifftshift(opticalOutElectricFieldSigFastFreq));
    OpticalElectricFieldFastAxisOut=OpticalElectricFieldFastAxisOut(1:end-4000);
    OpticalElectricFieldSlowAxisOut=ifft(ifftshift(opticalOutElectricFieldSigSlowFreq));
    OpticalElectricFieldSlowAxisOut=OpticalElectricFieldSlowAxisOut(1:end-4000);

%     case 'LinearODEFiber'
%         % Chromatic Dispersion and IL
%     beta2=-(FIBER.DispersionSlope*1e-12*FIBER.Wavelength.^3*1e-9*(1-(FIBER.ZeroDispersionWavelength./FIBER.Wavelength).^4)/(8*pi*3e8));
%     gamma=-1j*beta2*FIBER.Length/(2*Ts^2);
%     nTerms=140;
%     A=zeros(nTerms,1);
%     for k=0:nTerms
% if(mod(k,2)==0)
% A(k+1)=(gamma^(k/2)/factorial(k/2))*hypergeom((k+1)/2,0.5,gamma);
% else
% A(k+1)=-(k+1)*(gamma^((k+1)/2)/factorial((k+1)/2))*hypergeom(k/2 +1,1.5,gamma);
% end
%     end
%     for i=0:2:nTerms
%         for k=0:i
%             A(k)=A(k)+gamma.^(i/2)/factorial(i/2) * nchoosek(i,k);
%         end
%     end
%     OpticalElectricFieldAfterCD=conv(OpticalElectricFieldIn,A)*db2mag(-FIBER.Losses);
%     OpticalElectricFieldAfterCD=OpticalElectricFieldAfterCD(nTerms:end);



end


end
