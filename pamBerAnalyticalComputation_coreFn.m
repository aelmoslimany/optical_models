function [ber] = pamBerAnalyticalComputation_coreFn(pamOrder,SNR)
ber=(2*(pamOrder-1)/(pamOrder*log2(pamOrder))).*qfunc(sqrt(3*(SNR)./(pamOrder.^2-1)));
end

