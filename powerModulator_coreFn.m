function [pout] = powerModulator_coreFn(vin,vpp,ER,OMA)
PDC=(OMA/2)*(ER+1)/(ER-1);
pout=OMA*(vin/vpp)+PDC;
end
