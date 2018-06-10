function [A_bt, B_bt, C_bt, D_bt, HankSingVal] = BalancedTruncation(A, B, C, D, l)
%BalancedTruncation returns the truncated system of order l, computed using
%the balanced truncation method

S = ss(A, B, C, D); 

% Gramians
R = gram(S, 'c'); Q = gram(S, 'o');

% Decomposition Gramians
[V_ctl, D_ctl] = eig(R); 
[V_obs, D_obs] = eig(Q);

X = V_ctl*sqrt(D_ctl); % --> R = XX'
Y = V_obs*sqrt(D_obs); % --> Q = YY'
% Rem: X' is the conjugate transpose of the matrix X

% SVD of XY' and balancing transformation
[U, Sigma, V] = svd((Y')*X);
HankSingVal = svd((Y')*X);

% Careful, the last diagonal entries of Sigma are very close to zero and 
%considered to be zero when computing the rank of Sigma.
r = rank(Sigma); 

% Truncated system of order r:
% Sigma1r = Sigma(1:r, 1:r); U1r = U(:, 1:r); V1r = V(:, 1:r);
% T1r = (Sigma1r^(-1/2))*(U1r')*(Y'); S1r = X*V1r*(Sigma1r^(-1/2));
% Ar = T1r*A*S1r; Br = T1r*B; Cr = C*S1r; Dr = zeros(1, 1);
% Sr = ss(Ar, Br, Cr, Dr); % the minimal system (of order r)

% Model reduction

Sigma1 = Sigma(1:l, 1:l); U1 = U(:, 1:l); V1 = V(:, 1:l);
T1 = (Sigma1^(-1/2))*(U1')*(Y'); S1 = X*V1*(Sigma1^(-1/2));
%norm(T1*S1 - eye(l)) % in order to check

A_bt = T1*A*S1; B_bt = T1*B; C_bt = C*S1; D_bt = D;
end

