function [Asp,Bsp, A, B, C, D] = MatricesSemiDiscretizedPde(Nx, hx)
%MatricesSemiDiscretizedPde Returns the matrices A, B, C, D of the equation
%z' = Lapl z, in square. With 0 boundary contition execpt on one side, 
%Gamma,  where the values are given by the control. The output is the
%normal derivative of z on Gamma.

%Uses the finite difference method.

%The matrices Asp and Bsp are also returned in sparse form.

% A
Diag = sparse(1:Nx, 1:Nx, -4) + sparse(2:Nx, 1:Nx-1, 1, Nx, Nx) + sparse(1:Nx-1, 2:Nx, 1, Nx, Nx);
ExtraDiag = sparse(1:Nx, 1:Nx, 1);
Lapl = kron(eye(Nx),Diag) + kron(diag(ones(1, Nx-1),-1), ExtraDiag) + kron(diag(ones(1, Nx-1),1), ExtraDiag);
Lapl = 1/(hx*hx)*Lapl;
Asp = Lapl; A = full(Asp);

% B
Bsp = sparse(1:Nx, 1:Nx, 1/(hx*hx), Nx^2, Nx); B = full(Bsp);

% C
Csp = sparse(1:Nx, 1:Nx, -1/hx, Nx, Nx^2); C = full(Csp);

% D
D = 1/hx*eye(Nx);
end

