function [Vout] = dynamicNonlinearity_coreFn(Vin,A,omegaPoleTs,omegaPoleTs1,AR,AI,alphaTs,betaTs,alphaTs1,betaTs1)
Vout1stOrder=zeros(length(Vin),length(A));
Vout2ndOrder=zeros(length(Vin),length(AR));
Vout=zeros(size(Vin));
omegaPoleInst=omegaPoleTs(:,1);
alphaTsInst=alphaTs(:,1);
betaTsInst=betaTs(:,1);
B1=exp(-omegaPoleInst.');
B2=exp(-alphaTsInst.');
C3=cos(betaTsInst.');
S3=sin(betaTsInst.');
QR=(AR.').*(1-B2.*C3) - (AI.').*B2.*S3;
QI=(AI.').*(1-B2.*C3) + (AR.').*B2.*S3;

i=1;
Vout1stOrder(i,:)=(A.').*(1-B1).*Vin(i);
Vout2ndOrder(i,:)=2*QR.*Vin(i);
Vout(i)=sum(Vout1stOrder(i,:))+sum(Vout2ndOrder(i,:));

omegaPoleInst=omegaPoleTs(:,1)+omegaPoleTs1(:,1).*abs(Vout(i));
alphaTsInst=alphaTs(:,1)+alphaTs1(:,1).*abs(Vout(i));
betaTsInst=betaTs(:,1)+betaTs1(:,1).*abs(Vout(i));
B1=exp(-omegaPoleInst.');
B2=exp(-alphaTsInst.');
C3=cos(betaTsInst.');
S3=sin(betaTsInst.');
QR=(AR.').*(1-B2.*C3) - (AI.').*B2.*S3;
QI=(AI.').*(1-B2.*C3) + (AR.').*B2.*S3;

i=2;
Vout1stOrder(i,:)=B1.*Vout1stOrder(i-1,:)+(A.').*(1-B1).*Vin(i);
Vout2ndOrder(i,:)=2*B2.*C3.*Vout2ndOrder(i-1,:) + 2*QR.*Vin(i) - 2*B2.*(QR.*C3 - QI.*S3).*Vin(i-1);
Vout(i)=sum(Vout1stOrder(i,:))+sum(Vout2ndOrder(i,:));

omegaPoleInst=omegaPoleTs(:,1)+omegaPoleTs1(:,1).*abs(Vout(i));
alphaTsInst=alphaTs(:,1)+alphaTs1(:,1).*abs(Vout(i));
betaTsInst=betaTs(:,1)+betaTs1(:,1).*abs(Vout(i));
B1=exp(-omegaPoleInst.');
B2=exp(-alphaTsInst.');
C3=cos(betaTsInst.');
S3=sin(betaTsInst.');
QR=(AR.').*(1-B2.*C3) - (AI.').*B2.*S3;
QI=(AI.').*(1-B2.*C3) + (AR.').*B2.*S3;

for i=3:length(Vin)
    Vout1stOrder(i,:)=B1.*Vout1stOrder(i-1,:)+(A.').*(1-B1).*Vin(i);
    Vout2ndOrder(i,:)=2*B2.*C3.*Vout2ndOrder(i-1,:) - B2.^2 .* Vout2ndOrder(i-2,:) + 2*QR.*Vin(i) - 2*B2.*(QR.*C3 - QI.*S3).*Vin(i-1);
    Vout(i)=sum(Vout1stOrder(i,:))+sum(Vout2ndOrder(i,:));
    omegaPoleInst=omegaPoleTs(:,1)+omegaPoleTs1(:,1).*abs(Vout(i));
    alphaTsInst=alphaTs(:,1)+alphaTs1(:,1).*abs(Vout(i));
    betaTsInst=betaTs(:,1)+betaTs1(:,1).*abs(Vout(i));
    B1=exp(-omegaPoleInst.');
    B2=exp(-alphaTsInst.');
    C3=cos(betaTsInst.');
    S3=sin(betaTsInst.');
    QR=(AR.').*(1-B2.*C3) - (AI.').*B2.*S3;
    QI=(AI.').*(1-B2.*C3) + (AR.').*B2.*S3;
end
end
