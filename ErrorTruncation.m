function [error_matr] = ErrorTruncation(A, B, C, Al, Bl, Cl, freq, Nx)
%ErrorTruncation returns the error between the original system (A, B, C, D)
%and the truncated system (Al, Bl, Cl, Dl). As in our context D=Dl, the
%matrices D and Dl are not needed. The error returned is the matricial norm
%(ie. the largest singular value) of the difference between the transfer 
%function of the two systems, evaluated at the given frequencies freq.

nb_points = length(freq);
size_Al = size(Al); l = size_Al(1, 1);

% Error system
Ae = [Al, zeros(l, Nx^2); zeros(Nx^2, l), A];
Be = [Bl; B];
Ce = [-Cl, C];
%De = zeros(Nx, Nx);

% Error
error_matr = zeros(1, nb_points);
for k=1:nb_points
    error_matr(k) = norm( Ce * ( 1i*freq(k)*eye(Nx^2+l) - Ae )^(-1) * Be );
end

end

