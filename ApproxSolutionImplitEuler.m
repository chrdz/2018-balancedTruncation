function [z] = ApproxSolutionImplitEuler(Nx, Lx, Nt, ht, Asp, Bsp, input, z0)
%ApproxSolutionImplitEuler Computes and plot the solution to the full space
%semi-discretized problem, using the implicit euler scheme.

font = 15;

% Space axis
xaxis = linspace(0, Lx, Nx+2); 
%int_mesh_xaxis=meshgrid(x1axis(1, 2:Nx+1), x1axis(1, 2:Nx+1));
[x1_axis, x2_axis] = meshgrid(xaxis, xaxis);

% Initial condition: heat the middle of the square
z = zeros(Nx^2, Nt+2); % rows<->space columns<->time
z(:,1) = reshape(z0', Nx^2, 1);

% Matrix to invert
M = sparse(1:Nx^2, 1:Nx^2, 1/ht) - Asp;
[L, U, P, Q] = lu(M);

% Implict Euler Scheme
for k=1:Nt+1
  fk = 1/ht*z(:, k) + Bsp*input(:, k+1);
  z(:, k+1) = Q*(U\(L\(P*fk))); % z(:, k+1) = M\fk;
end

%test
% zfull = zeros(Nx+2, Nx+2);
% test = z(:, 1); test(1, 1)=0.03; test(Nx, 1)=0.04; test(Nx^2-Nx+1, 1)=0.05; test(Nx^2, 1)=0.06;
% zreshape=reshape(test, Nx, Nx)';
% zfull(2:Nx+1, 2:Nx+1)=zreshape;
% surf(x1_axis, x2_axis, zfull);
% zlim([0 0.06]);
% xlabel('x1'); ylabel('x2'); zlabel('z(t, x1, x2)')
% title('Approximate solution of the PDE, at fixed time t')
% for k=2:Nt+2
%     test = z(:, k); test(1, 1)=0.03; test(Nx, 1)=0.04; test(Nx^2-Nx+1, 1)=0.05; test(Nx^2, 1)=0.06;
%     zreshape=reshape(test, Nx, Nx)';
%     zfull(2:Nx+1, 2:Nx+1)=zreshape;
%     surf(zfull);
%     zlim([0 0.06]);
%     drawnow;
%     pause(0.1);
% end
%%% To plot :
figure;
zfull = zeros(Nx+2, Nx+2);
zreshape=reshape(z(:, 1), Nx, Nx)';
zfull(2:Nx+1, 2:Nx+1) = zreshape;  
surf(x1_axis, x2_axis, zfull);
zlim([0 1]);
xlabel('x1'); ylabel('x2'); zlabel('z(t, x1, x2)')
title('Approximate solution, at time t')
liste = [5 15 25 35];
save1 = zeros(Nx+2, Nx+2); save2 = zeros(Nx+2, Nx+2); 
save3 = zeros(Nx+2, Nx+2); save4 = zeros(Nx+2, Nx+2);
for k=2:Nt+2
    zreshape=reshape(z(:, k), Nx, Nx)';
    zfull(2:Nx+1, 2:Nx+1) = zreshape;  
    
    if k==liste(1)
        save1 = zfull;
    end
    if k==liste(2)
        save2 = zfull;
    end
     if k==liste(3)
        save3 = zfull;
    end   
    if k==liste(4)
        save4 = zfull;
    end
    
    surf(x1_axis, x2_axis, zfull);
    zlim([0 1]);
    xlabel('x1'); ylabel('x2'); zlabel('z(t, x1, x2)')
    title('Approximate solution, at time t')
    set(gca,'fontsize',font)
    drawnow;
    pause(0.1);
end

figure;
surf(x1_axis, x2_axis, save1);
zlim([0 1]);
xlabel('x1'); ylabel('x2'); zlabel('z(t, x1, x2)')
title('Approximate solution, iter=5')
set(gca,'fontsize',font)

figure;
surf(x1_axis, x2_axis, save2);
zlim([0 1]);
xlabel('x1'); ylabel('x2'); zlabel('z(t, x1, x2)')
title('Approximate solution, iter=15')
set(gca,'fontsize',font)

figure;
surf(x1_axis, x2_axis, save3);
zlim([0 1]);
xlabel('x1'); ylabel('x2'); zlabel('z(t, x1, x2)')
title('Approximate solution, iter=25')
set(gca,'fontsize',font)

figure;
surf(x1_axis, x2_axis, save4);
zlim([0 1]);
xlabel('x1'); ylabel('x2'); zlabel('z(t, x1, x2)')
title('Approximate solution, iter=35')
set(gca,'fontsize',font)

end

