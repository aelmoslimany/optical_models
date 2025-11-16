function [R3,R4,R5,R6,R7,R11,R12,R13,R14,L2,L3,L6,L7,C0,C1,C3,C10,ng] = MZMPhaseShifterCircuitModel(l,v)
rsub=0.135;%% Effective substrate resistivity [ohm*m]
temperature = 300.15; %Temperature (K)
indperl = 440e-9;
indlperl = 50e-9;
resperl = 850;
reshperl = 2000;
indrsxperl = 900e3*exp(1.8e4*(l-0.01e-3));
cmetperl = 120e-12;
cansxperl = 100e-12/2;
rsx = 1/l;
csx = 8.8542e-12*11.9*rsub*l/3;
vbi_va =  0.825;%%built in voltage of pn junction
m_va = 0.3;  %% grading coefficient (related to doping profile)
ngpnte = 3.94;%%index group TE mode
condperl_p = 195;
condperl_n = 585;
capperl  =  295e-12; %%C/m
co_mod=1;
dtemp = 0;
deltemp = temperature + dtemp - 298.15;
ind = indperl*l;
indrsx = indrsxperl*l;
cap = capperl*(1 + 0.0001*deltemp)*l*co_mod;%%co_mod-->corner effect
res = (3.0*deltemp+resperl)*l;
resg_p = 1.0/(condperl_p*l);
resg_n = 1.0/(condperl_n*l);
cmet = cmetperl*l;
cansx = cansxperl*l;
vbi_vat= vbi_va-0.0015*deltemp;
if (v < vbi_vat)
    cj= cap*exp(-m_va*log((vbi_vat-v)/vbi_va));
end
%%1st section
R3 = 0.5*res;
L2 = 0.5*ind;
R5 = 0.5*indrsx;
R4 =  0.5*reshperl*l;
L3 = 0.5*indlperl*l;
C1 = cmet;
%%Junction
C0 = cj;
R6 = resg_p*(1+ 0.001 * deltemp);
R7 = resg_n*(1+ 0.001 * deltemp);
%%3rd section
R12 = 0.5*res;
L6 = 0.5*ind;
R14 = 0.5*indrsx;
R13 = 0.5*reshperl*l;
L7 = 0.5*indlperl*l;
C3 = cansx;
R11 = rsx;
C10 = csx;
ng = ngpnte;
end
