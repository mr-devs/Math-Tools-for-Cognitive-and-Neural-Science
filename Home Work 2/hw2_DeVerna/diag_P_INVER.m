function [diag_transposed] = diag_P_INVER(inputDIAG)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

mat = inputDIAG;

vector_check = isvector(mat) ;

if vector_check == 1
    diag_transposed = mat' ;
    diag_transposed(1,1) = 1/mat(1,1) ;
    
else
    diag_vectorized = diag(mat,0) ;
    
    length_elements = length(diag_vectorized) ;
    diag_transposed = inputDIAG' ;
    
    % Calculate their individual inverse values
    for ii = 1:length_elements
        diag_vectorized(ii) = 1/diag_vectorized(ii) ;
        diag_transposed(ii,ii) = diag_vectorized(ii) ;
    end
end

end

