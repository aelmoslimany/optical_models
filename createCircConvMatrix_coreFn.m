function [X] = createCircConvMatrix_coreFn(x,nCol)
X=toeplitz(x,[x(1) flip(x(end-nCol+2:end).')]);
end
