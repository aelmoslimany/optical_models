function [impulseResponse,causalityCorrectiondB]=enforceCausality_coreFn(impulseResponse,threshold1,threshold2)
% Causality is imposed using the Alternating Projections Method. See also:
% Quatieri and Oppenheim, "Iterative Techniques for Minimum Phase Signal
% Reconstruction from Phase or Magnitude", IEEE Trans. ASSP-29, December
% 1981 (http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1163714)

%%
absImpulseResponseFreq=abs(fft(impulseResponse));
impulseResponseLen = length(impulseResponse);
originalImpulseResponse=impulseResponse;
% Assumption: peak of impulse_response is in the first half, i.e. not anti-causal
absImpulseResponse=abs(impulseResponse);
maxAbsImpulseResponse=max(absImpulseResponse(1:ceil(impulseResponseLen/2)));
ind = find(absImpulseResponse > maxAbsImpulseResponse*threshold1);
startInd = ind(1);
endInd= ind(end);
%% make the start monotonically increasing/decreasing Mohammed Abdelghany
% Monotonic start
while ((startInd>1) && (sign(impulseResponse(startInd-1))==sign(impulseResponse(startInd))) && (absImpulseResponse(startInd)>maxAbsImpulseResponse*threshold2))
    startInd=startInd-1;
end
% Monotonic end
while ((endInd<(floor(impulseResponseLen/2)+startInd)) && sign(impulseResponse(endInd+1))==sign(impulseResponse(endInd)) && (absImpulseResponse(endInd)>maxAbsImpulseResponse*threshold2))
    endInd=endInd+1;
end
%%
err=inf;
while ~all(impulseResponse==0)
    impulseResponse(1:startInd-1)=0;
    impulseResponse(endInd+1:end)=0;
    IL_modified=absImpulseResponseFreq.*exp(1j*angle(fft(impulseResponse)));
    ir_modified = ifft(IL_modified);
    delta = abs(impulseResponse-ir_modified);
    err_prev = err;
    err=max(delta)/max(impulseResponse);
    if err<0.01 || abs(err_prev-err)<0.001
        break;
    end
    impulseResponse=ir_modified;
end

causalityCorrectiondB=20*log10(norm(impulseResponse-originalImpulseResponse)/norm(impulseResponse));
