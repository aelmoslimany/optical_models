function ID_dark = DarkCurrent(Vdc_PD,PD_width,Temperature)
Vac_M = Vdc_PD;
wi_M = PD_width;
dtemp_M = Temperature-25;
Pq = 1.60217663e-19;
Pk = 1.380649e-23;
dctol = 0.955;
temperature = 298.15;
wi = wi_M;
dtemp = dtemp_M;
Vac = Vac_M;
deltemp = temperature + dtemp - 298.15;
nomtemp = temperature + dtemp;
delwi= wi-0.7e-6;
v = Vac;
if (v < 0)
   alf = -0.184 - 0.24e6*delwi;
   bet =1.3e-3 - 3e10*delwi*delwi;
   voinv = alf + bet*deltemp;

   if (wi==0.5e-6)
       curr = -16e-8*exp(0.035*deltemp)*sinh(1*v*1.09*voinv);
   else
       curr = -16e-8*exp(0.043*deltemp)*sinh(1*v*voinv);
   end

else
   neff = 1.53 - 1.8e-3*deltemp + 0.2e-4*deltemp*deltemp;
   curr = 4e-9*exp(0.048*deltemp)*(exp(Pq*v/(neff*Pk*nomtemp))-1.0);
end
dkcur = 1.0*dctol*curr;
dkcur_M = dkcur;
ID_dark = dkcur_M;
end
