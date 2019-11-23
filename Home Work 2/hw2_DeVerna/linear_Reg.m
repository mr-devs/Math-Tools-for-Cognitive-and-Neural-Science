function beta_opt = linear_Reg(XX,y)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[U,S,V] = svd(XX)                           ;
Y_star = U'*y                               ;
b_ss = [Y_star(1:length(diag(S)))]          ;
b_star = b_ss(1:length(diag(S)))./diag(S)   ;
beta_opt = V*b_star                         ;

end

