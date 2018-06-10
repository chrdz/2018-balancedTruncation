function [norm_matr] = NormTransferFunction(A, B, C, D, freq)
% NormTransferFunction returns the matrix norm of the transfer function at
% frequencies frec

nb_points = length(freq);
size_A = size(A); N = size_A(1, 1);

% Error
norm_matr = zeros(1, nb_points);
for k=1:nb_points
    norm_matr(k) = norm( C * ( 1i*freq(k)*eye(N) - A )^(-1) * B + D );
end

end

