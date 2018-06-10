function [lower, error, upper] = InfinityError(A, B, C, D, Nx, l_val)
%UNTITLED14 Summary of this function goes here
%   Detailed explanation goes here
S = ss(A, B, C, D);
nb_points = length(l_val);
lower = zeros(nb_points, 1); error = zeros(nb_points, 1); upper = zeros(nb_points, 1);
for k = 1:nb_points
    [A_bt, B_bt, C_bt, D_bt, HankSingVal] = BalancedTruncation(A, B, C, D, l_val(k));
    S_bt = ss(A_bt, B_bt, C_bt, D_bt);
    E_bt = S - S_bt;
    lower(k) = HankSingVal(l_val(k)+1);
    upper(k) = 2 * sum(HankSingVal(l_val(k)+1:Nx^2));
    error(k) = norm(E_bt, Inf);
end

