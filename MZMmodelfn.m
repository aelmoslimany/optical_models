function [Outputpower,OutputField] = MZMmodelfn(Vin,Ts,Vbias,ModulatorLength,N_sections,Zsource,ZLoad,ArmMismatch,Opticalpower,HeaterPhase,Wavelength)
Vin=Vin(:).';
SpeedofLight = 299792458;
N=length(Vin);
freq=linspace(-1/(2*Ts),1/(2*Ts),N+1);
freq=freq(1:end-1);
freqPostive=flip(-freq(freq<=0));
omega = 2*pi*freqPostive;
VinFreq=fft(Vin);
VsourceFreq=VinFreq(1:length(freqPostive));
SectionLength=ModulatorLength/N_sections; %% from [10u:0.4m]; %%Section length
[R3,R4,R5,R6,R7,R11,R12,R13,R14,L2,L3,L6,L7,C0,C1,C3,C10,ng] = MZMPhaseShifterCircuitModel(SectionLength,Vbias);%%TL parameters SectionLength: setction length, v bias voltage
ExcessLoss = 10^(-0.2950/20);% bend and heaters Transmission
%%%%%%%%%1st section
Z1 = 1./(1./(R3+1i*omega*L2) + 1/R5) + 1./(1./(1i*omega*L3) + 1/R4);
Z2 = 0;
Z3 = 1./(1i*omega*C1);
%%ABCD
A1 = 1 + Z1./Z3;
B1 = Z1 + Z2 + (Z1.*Z2)./Z3;
C1 = 1./Z3;
D1 = 1 + Z2./Z3;
T1(1,1,:)=A1;
T1(1,2,:)=B1;
T1(2,1,:)=C1;
T1(2,2,:)=D1;
%%%%%%%%%%3rd section
Z1 = 0;
Z2 = 1./( 1/R13 + 1./(1i*omega*L7)) + 1./(1./(R12 + 1i*omega*L6) + 1/R14);
Z3 = 1./( 1/R11 + (1i*omega*C10) )  + 1./(1i*omega*C3);
%%ABCD
A3 = 1 + Z1./Z3;
B3 = Z1 + Z2 + (Z1.*Z2)./Z3;
C3 = 1./Z3;
D3 = 1 + Z2./Z3;
T3(1,1,:)=A3;
T3(1,2,:)=B3;
T3(2,1,:)=C3;
T3(2,2,:)=D3;
%%%%%%%%%%%2nd section (Junction)
Ra = R6 + R7;
La = 0;
Ca = C0;
%%ABCD
A2 = 1;
B2 = 0;
C2 = 1./(R6 + R7 + 1i*(omega*La-(1./(omega.*Ca))));
D2 = 1;
T2(2,1,:)=C2;
T2(2,2,:)=D2;
T2(1,1,:)=A2;
T2(1,2,:)=B2;
%%%%%%%%%%%%%%%%%%%%%Total ABCD%%%%%%%%%%%%%%%%%
temp(1,1,:) = ones(1,length(A3));
temp(1,2,:) =zeros(1,length(A3));
temp(2,1,:) = zeros(1,length(A3));
temp(2,2,:) = ones(1,length(A3));
for Index = 1:N_sections
    temp=pagemtimes(T3,temp);
    temp=pagemtimes(T2,temp);
    temp=pagemtimes(T1,temp);
end
A = squeeze(temp(1,1,:)).';
B = squeeze(temp(1,2,:)).';
C = squeeze(temp(2,1,:)).';
D = squeeze(temp(2,2,:)).';

IL = VsourceFreq./(A.*ZLoad + B + C.*Zsource.*ZLoad + D.*Zsource);%%Load/Termination current
VL = ZLoad.*IL;%%Load/Termination Voltage

%% ABCD to S-parameters

Vn = zeros(N_sections,length(freqPostive));
TF = zeros(N_sections,length(freqPostive)); % Vn = TF.*VL

temp(1,1,:) = ones(1,length(A3));
temp(1,2,:) =zeros(1,length(A3));
temp(2,1,:) = zeros(1,length(A3));
temp(2,2,:) = ones(1,length(A3));

for Index = 1:N_sections
    temp=pagemtimes(T3,temp);
    Vn(N_sections+1-Index,:) = squeeze(temp(1,1,:)).'.*VL + squeeze(temp(1,2,:)).'.*IL;
    TF(N_sections+1-Index,:) = squeeze(temp(1,1,:)).'*ZLoad + squeeze(temp(1,2,:)).';
    temp=pagemtimes(T2,temp);
    temp=pagemtimes(T1,temp);
end

xi = 1./(1 - (omega.^2).*La.*Ca + 1i*omega.*Ra.*Ca);
Vc = xi.*Vn;
 % Vc2=xi.*TF.*VL;
Vc2=(xi.*TF./(A.*ZLoad + B + C.*Zsource.*ZLoad + D.*Zsource)) .* VsourceFreq;
z = 0.5*SectionLength:SectionLength:(ModulatorLength-0.5*SectionLength);
PhaseMismatch=bsxfun(@times,omega,z.'.*ng/SpeedofLight);
VoltageCap = Vc.*exp(1i*PhaseMismatch);
TF2=((xi.*TF./(A.*ZLoad + B + C.*Zsource.*ZLoad + D.*Zsource)).*exp(1i*PhaseMismatch));
VoltageCap2 = TF2 .* VsourceFreq;

PhaseUpperArm = zeros(1,N);
PhaseLowerArm = zeros(1,N);
LossUpperArm = zeros(1,N);
LossLowerArm = zeros(1,N);
for Index = 1:(N_sections)
    temp=[VoltageCap(Index,:) fliplr(conj(VoltageCap(Index,2:end-1)))];
    temp=real(ifft(temp));
    [phasesegment,losssegment] = MZMPhaseShifterModel(SectionLength,temp+Vbias,Wavelength); %%Phase per section
    PhaseUpperArm = PhaseUpperArm+phasesegment;
    LossUpperArm = LossUpperArm+losssegment;
    [phasesegment,losssegment] = MZMPhaseShifterModel(SectionLength,-temp+Vbias,Wavelength); %%Phase per section
    PhaseLowerArm = PhaseLowerArm+phasesegment;
    LossLowerArm = LossLowerArm+losssegment;
end
PhaseLowerArm=(100-ArmMismatch*100).*PhaseLowerArm./100;
LossLowerArm=(100-ArmMismatch*100).*LossLowerArm./100;





ElectricField = [ sqrt(0.5*Opticalpower).*exp(1j*(PhaseUpperArm+HeaterPhase)).*exp(-LossUpperArm).*ExcessLoss ; sqrt(0.5*Opticalpower).*exp(1j*PhaseLowerArm).*exp(-LossLowerArm).*ExcessLoss];
OutputField = sqrt(0.5).*(ElectricField(1,:) + ElectricField(2,:));
Outputpower = abs(OutputField(1,:)).^2;

end
