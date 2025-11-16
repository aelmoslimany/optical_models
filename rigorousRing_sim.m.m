clear
clc
close all

n_Si=3.92;%3.92;


% Pin=7;%dBm
% R=3.8e-6;%m
% BWIL=5;%dB
% kappa2=5.1/100;%
% loss=-0.65*0.62;%dB per roundtrip (normalized for 7.5um ring)
% ME=88e-12;%pm/v (for 2V swing) [Should actually equal to 96.37 (measured) 92.3 (simulated)]
% RCBW=60;%GHz

% Pin=7;%dBm
% R=4e-6;%m
% BWIL=5;%dB
% kappa2=6.4/100;%
% loss=-0.65*0.69;%dB per roundtrip (normalized for 7.5um ring)
% ME=42e-12;%pm/v (for 2V swing) [Should actually equal to 96.37 (measured) 92.3 (simulated)]
% RCBW=82;%GHz

% %%%%%%%%%%%%%%%%GF
Pin=7;%dBm
R=7.5e-6;%m
BWIL=5;%dB
kappa2=12.69/100;%
loss=-0.65;%dB per roundtrip (normalized for 7.5um ring)
ME=90e-12;%pm/v (for 2V swing) [Should actually equal to 96.37 (measured) 92.3 (simulated)]
RCBW=82;%GHz

% Pin=7;%dBm
% R=18.48e-6;%m
% BWIL=5;%dB
% kappa2=12.69/100;%
% loss=-0.65;%dB per roundtrip (normalized for 7.5um ring)
% ME=12e-12*2;%pm/v (for 2V swing) [Should actually equal to 96.37 (measured) 92.3 (simulated)]
% RCBW=82;%GHz

% Pin=7;%dBm
% R=5e-6;%m
% BWIL=5;%dB
% kappa2=11/100;%
% loss=-0.535;%dB per roundtrip (normalized for 7.5um ring)
% ME=21e-12;%pm/v (for 2V swing) [Should actually equal to 96.37 (measured) 92.3 (simulated)]
% % ME=30.1386e-12*2;%pm/v (for 2V swing) [Should actually equal to 96.37 (measured) 92.3 (simulated)]
% RCBW=89;%GHz


% Pin=7;%dBm
% R=12e-6;%m
% BWIL=5;%dB
% kappa2=20/100;%
% loss=-0.65*12/7.5*0.45*1;%dB per roundtrip (normalized for 7.5um ring)
% ME=33e-12;%pm/v (for 2V swing) [Should actually equal to 96.37 (measured) 92.3 (simulated)]
% % ME=30.1386e-12*2;%pm/v (for 2V swing) [Should actually equal to 96.37 (measured) 92.3 (simulated)]
% RCBW=54;%GHz

% Pin=7;%dBm
% R=7.5e-6;%m
% BWIL=5;%dB
% kappa2=12.69/100;%14/100;%
% loss=-0.65;%dB
% ME=90e-12*24.7/31.17;%pm/v (for 2V swing)
% RCBW=82*132.56/90.33;%GHz

% Pin=7;%dBm
% R=15e-6;%m
% BWIL=5;%dB
% kappa2=23.2/100;%
% loss=-0.65;%dB
% ME=90e-12;%pm/v (for 2V swing)
% RCBW=82;%GHz

% Pin=7;%dBm
% R=15e-6;%m
% BWIL=4.5;%dB
% kappa2=24/100;%
% loss=-0.65;%dB
% ME=90e-12;%pm/v (for 2V swing)
% RCBW=82;%GHz

% Pin=7;%dBm
% R=7.5e-6;%m
% BWIL=6;%dB
% kappa2=11/100;%
% loss=-0.65;%dB per roundtrip
% ME=90e-12;%pm/v (for 2V swing) [Should actually equal to 96.37 (measured) 92.3 (simulated)]
% RCBW=82;%GHz

% Pin=6;%dBm
% R=20e-6;%m
% BWIL=4;%dB
% kappa2=30/100;%
% loss=-0.65;%dB per roundtrip
% ME=96e-12;%pm/v (for 2V swing) [Should actually equal to 96.37 (measured) 92.3 (simulated)]
% RCBW=82;%GHz

% load MAGcurve
% plot(lambda,pow2db(IL))
% hold on
% lambda_max=1320e-9;
% lambda_min=1308e-9;
% lambda_step=1e-14;

n_g=n_Si;

lambda=linspace(1311e-9-(1.31e-6)^2/(n_g*2*pi*R)/2,1311e-9+(1.31e-6)^2/(n_g*2*pi*R)/2,100000);
lambda=fliplr(lambda);

sigma=sqrt(1-kappa2);%sqrt(1-0.11);
A=sqrt(10^(loss*R/7.5e-6/10));%0.9483*1.015;%0.949;%0.96;

c=3e8;

% er=0;
% for sigma=0.92:0.0001:0.98
%     er=er+1;

m=0;
% for phi_initial=-(1.6:0.0001:1.7)
%     m=m+1

alpha=A;%sqrt(1-A^2);
% lambda=lambda_min:lambda_step:lambda_max;
%lambda=Data1(:,1);
% phi_initial=-1.6558-4.02662;
% phi=2*pi*R*2*pi./(lambda)*n_Si + phi_initial;
phi_initial=-2*pi*R*2*pi./(1.311e-6)*n_Si;
phi=2*pi*R*2*pi./(lambda)*n_Si + phi_initial;
Ts=(sigma-A*exp(-1j*phi))./(1-sigma*A*exp(-1j*phi));

%Ts111=sqrt((sigma.^2+A.^2 -2*sigma*A*cos(phi))./(1+sigma.^2*A.^2 -2*sigma*A*cos(phi)));


Ts_min=min(Ts);
plot(lambda,pow2db(abs(Ts).^2))
figure
lambdao=1.310788287767544e-06;
% lambdao=1.310856559783117e-06;

for i=1:length(Ts)
    if 10*log10(abs(Ts(i)).^2)<=-BWIL
    lambdai=lambda(i);
    break;
    end
end

disp( [ 'FSR = ' num2str( 1e9*(1.31e-6)^2/(n_g*2*pi*R) ) ' nm'])

disp( [ 'IL = ' num2str( 10*log10(abs(Ts(i)).^2) ) ' dB'])

ERphi1=2*pi*R*2*pi./(lambdai-ME*0.5471)*n_Si + phi_initial;
ERphi2=2*pi*R*2*pi./(lambdai+ME*0.4529)*n_Si + phi_initial;
ERTs1=(sigma-A*exp(-1j*ERphi1))./(1-sigma*A*exp(-1j*ERphi1));
ERTs2=(sigma-A*exp(-1j*ERphi2))./(1-sigma*A*exp(-1j*ERphi2));
ER = max([abs(ERTs1^2) abs(ERTs2^2)])/min([abs(ERTs1^2) abs(ERTs2^2)]);

disp(['ER = ' num2str(10*log10(ER)) ' dB'])




OMA=10^(Pin/10)*(max([abs(ERTs1^2) abs(ERTs2^2)])-min([abs(ERTs1^2) abs(ERTs2^2)]));

disp(['OMA = ' num2str(10*log10(OMA)) ' dBm'])

RLMphi1=2*pi*R*2*pi./(lambdai-ME*0.5471)*n_Si + phi_initial;
RLMphi2=2*pi*R*2*pi./(lambdai-ME*0.1679)*n_Si + phi_initial;
RLMphi3=2*pi*R*2*pi./(lambdai+ME*0.1582)*n_Si + phi_initial;
RLMphi4=2*pi*R*2*pi./(lambdai+ME*0.4529)*n_Si + phi_initial;
RLMTs1sq=abs((sigma-A*exp(-1j*RLMphi1))./(1-sigma*A*exp(-1j*RLMphi1)))^2;
RLMTs2sq=abs((sigma-A*exp(-1j*RLMphi2))./(1-sigma*A*exp(-1j*RLMphi2)))^2;
RLMTs3sq=abs((sigma-A*exp(-1j*RLMphi3))./(1-sigma*A*exp(-1j*RLMphi3)))^2;
RLMTs4sq=abs((sigma-A*exp(-1j*RLMphi4))./(1-sigma*A*exp(-1j*RLMphi4)))^2;
Vmin=(RLMTs1sq+RLMTs4sq)/2;
ES1=(RLMTs2sq-Vmin)/(RLMTs1sq-Vmin);
ES2=(RLMTs3sq-Vmin)/(RLMTs4sq-Vmin);
RLM=min([(3*ES1) (3*ES2) (2-3*ES1) (2-3*ES2)]);

VpiL=2*pi/(ERphi1-ERphi2)/(1e-2/(R*2*pi));
disp(['VpiL = ' num2str(VpiL) ' V.cm'])

disp(['RLM = ' num2str(RLM)])

for i=1:length(Ts)
    if abs(Ts(i))^2<=0.5
    lambdaQ1=lambda(i);
    Qi=i;
    break;
    end
end

for i=Qi:length(Ts)
    if abs(Ts(i))^2>0.5
    lambdaQ2=lambda(i);
    break;
    end
end

disp(['Q = ' num2str(1.31e-6/(abs(lambdaQ1-lambdaQ2)))])

finesse=(1.31e-6)^2/(n_g*2*pi*R)/(abs(lambdaQ1-lambdaQ2));

disp(['Finesse = ' num2str(finesse)])

lambdao=lambdai;
T_s=2*pi*R*n_g/c;
EOr=0;
mm=0;
for f=0.1e9:0.1e9:200e9
    mm=mm+1;
    t=0:T_s:(1/f);
    phi=2*pi*R*2*pi./(lambdao)*n_Si + phi_initial;
    Ts=(sigma-A*exp(-1j*phi))./(1-sigma*A*exp(-1j*phi));
    m=1;
    for t=t
      m=m+1;
      phi=2*pi*R*2*pi./(lambdao+1e-13*sin(2*pi*f*t))*n_Si + phi_initial;
      x=A*exp(-1j*phi);
      Ts(m)=sigma-x+Ts(m-1)*sigma*x;
    end
    EOr(mm)=max(abs(Ts.^2))-min(abs(Ts.^2));
   
end
hold on
plot((0.1e9:0.1e9:200e9)/1e9,mag2db(EOr/EOr(1)))
f=0.1e9:0.1e9:200e9;
for i=1:length(EOr)
    if mag2db(EOr(i)/EOr(1))<=-3
    freqi=f(i);
    break;
    end
end

disp(['Optical BW = ' num2str(freqi/1e9) ' GHz'])
disp(['RC BW = ' num2str(RCBW) ' GHz'])

BWtot=sqrt(1/(1/(freqi/1e9)^2+(1/RCBW)^2 ));

disp(['Total BW = ' num2str(BWtot) ' GHz'])

gammaKappa=kappa2/(2*pi*R*n_g/3e8);
gammaLoss=(1-db2pow(loss*R/7.5e-6))/(2*pi*R*n_g/3e8);

disp(['gamma_Kappa = ' num2str(gammaKappa)])
disp(['gamma_Loss = ' num2str(gammaLoss)])

    %  t=0:T_s:(1/f);
    %  figure
    % plot(t,abs(Ts(2:end)).^2)

% % 
% % % hold on
% % Pin=10^(-20/10)*1e-3;
% % % plot(Data1(:,1),(Data1(:,2)/max(Data1(:,2))))
% % %plot(Data1(:,1),10*log10(Data1(:,2)/max(Data1(:,2))))
% % 
% % hold on
% % % [a b]=min(Ts);
% % % [c d]=min(Data1(:,2));
% % % delta(m)=lambda(b)-Data1(d,1);
% % % disp([lambda(b)-Data1(d,1)])
% % % end
% % 
% % gamma_k=(1-sigma^2)*c/(2*pi*R*n_g);
% % 
% % % hold off
% % % a=( sqrt(Pin)-sqrt(Pin)*Ts )/sqrt(gamma_k);
% % % plot(lambda,abs(a).^2*4.8783e+13)
% % 
% % clear Ts
% % m=0;
% % for lambda=lambda_min:lambda_step:lambda_max
% %     m=m+1;
% %     if m==1
% %         phi=2*pi*R*2*pi./lambda*n_Si + phi_initial;
% %         Ts=(sigma-A*exp(-1j*phi))./(1-sigma*A*exp(-1j*phi));
% %     else
% %         % a=( sqrt(Pin)-sqrt(Pin*abs(Ts(m-1)).^2) )/sqrt(gamma_k);
% %         a=( sqrt(Pin)-Ts(m-1)*sqrt(Pin))/sqrt(gamma_k);
% %         phi=2*pi*R*2*pi./lambda*n_Si + phi_initial + 0.63*abs(a)^2*9e25/7.5e10*1.8e-4*2*pi*R*2*pi/lambda;
% %         % Ts(m)=min([(sigma-A*exp(-1j*phi))./(1-sigma*A*exp(-1j*phi)) Ts(m-1)]);
% %         % Ts(m)=(sigma-A*exp(-1j*phi))./(1-sigma*A*exp(-1j*phi));
% %         x=A*exp(-1j*phi);
% %         Ts(m)=sigma-x+Ts(m-1)*sigma*x;
% %     end
% %     % if abs(Ts(m)-Ts_min)<0.01
% %     %     phi=2*pi*R*2*pi./lambda*n_Si + phi_initial;
% %     %     Ts(m)=(sigma-A*exp(-1j*phi))./(1-sigma*A*exp(-1j*phi));
% %     % end
% % 
% % end
% % 
% % lambda=lambda_min:lambda_step:lambda_max;
% %  plot(lambda,(abs(Ts).^2))
% % %plot(lambda,10*log10(abs(Ts).^2))
% % 
% %  xlabel("\lambda (nm)")
% %  ylabel("Normalized Transmission")
% % % ylabel("Normalized Transmission Difference")
% % % legend('GF PDK 0dBm','Rigorous 0dBm','GF PDK 7dBm','Rigorous 7dBm')
% % set(gca,'FontSize',18)
% % % axis([1.31e-6 1.315e-6 -30 0])
% % 
% % % err(er)=sum( abs( abs(Ts).^2-Data1(:,2)/max(Data1(:,2)) ));
% % % end
% % % plot(err)
% % % disp(min(err))

