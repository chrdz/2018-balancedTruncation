clear; clc;

font = 18;
width = 1.5

% Space discretization
Lx = pi;
Nx = 10; % number of interior points
hx = Lx/(Nx+1);

% Space semi-discretizes PDE
[Asp, Bsp, A, B, C, D] = MatricesSemiDiscretizedPde(Nx, hx);

% Controllability and observability tests
[Abar_c,Bbar_c,Cbar_c,T_c,k_c] = ctrbf(A,B,C); 
nb_controllable_states = sum(k_c) % if sum(k_c) = Nx^2 ==> controllable 
[Abar_o,Bbar_o,Cbar_o,T_o,k_o] = obsvf(A,B,C); 
nb_observable_states = sum(k_o) % if sum(k_o) = Nx^2 ==> observable 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Approximation of the solution %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Time discretization
tmax = 1;
Nt = 200;
ht = tmax/(Nt+1);
taxis = linspace(0, tmax, Nt+2);

% Control
input = ones(Nx, Nt+2)*0;

% Initial condition
i1 = floor(2*Nx/5.0);  i2 = floor(3*Nx/5.0);
z0 = zeros(Nx, Nx); z0(i1+1:i2, i1+1:i2) = 1;

% Computation and plot of the solution before reduction:
[z] = ApproxSolutionImplitEuler(Nx, Lx, Nt, ht, Asp, Bsp, input, z0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Balanced Truncation %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S = ss(A, B, C, D); 

% Balanced truncation:
l_vec = [18 20 22];

% 1st value of l:
[A_bt_1, B_bt_1, C_bt_1, D_bt_1, HankSingVal] = BalancedTruncation(A, B, C, D, l_vec(1));
S_bt_1 = ss(A_bt_1, B_bt_1, C_bt_1, D_bt_1);

% 2nd value of l:
[A_bt_2, B_bt_2, C_bt_2, D_bt_2, temp] = BalancedTruncation(A, B, C, D, l_vec(2));
S_bt_2 = ss(A_bt_2, B_bt_2, C_bt_2, D_bt_2);

% 3rd value of l:
[A_bt_3, B_bt_3, C_bt_3, D_bt_3, temp] = BalancedTruncation(A, B, C, D, l_vec(3));
S_bt_3 = ss(A_bt_3, B_bt_3, C_bt_3, D_bt_3);

% 4th value of l:
%[A_bt_4, B_bt_4, C_bt_4, D_bt_4, temp] = BalancedTruncation(A, B, C, D, l_vec(4));
%S_bt_4 = ss(A_bt_4, B_bt_4, C_bt_4, D_bt_4);

% Hankel singular values:
figure; 
bar(HankSingVal, 'LineWidth', 1, 'BarLayout', 'stacked', 'FaceColor', 'b', 'EdgeAlpha', 0);
title('Hankel singular values');
set(gca,'fontsize',font)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Errors %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % H infinity error
% E_bt_1 = S-S_bt_1; error_approx_1 = norm(E_bt_1, Inf);
% E_bt_2 = S-S_bt_2; error_approx_2 = norm(E_bt_2, Inf);
% E_bt_3 = S-S_bt_3; error_approx_3 = norm(E_bt_3, Inf);
% E_bt_4 = S-S_bt_4; error_approx_4 = norm(E_bt_4, Inf);
% 
% % Upper bound H infinity error
% error_upper_bound_1 = 2*sum(HankSingVal(l_vec(1)+1:Nx^2));
% error_upper_bound_2 = 2*sum(HankSingVal(l_vec(2)+1:Nx^2));
% error_upper_bound_3 = 2*sum(HankSingVal(l_vec(3)+1:Nx^2));
% error_upper_bound_4 = 2*sum(HankSingVal(l_vec(4)+1:Nx^2));
% 
% % Lower bound H infinity error
% error_lower_bound_1 = HankSingVal(l_vec(1)+1);
% error_lower_bound_2 = HankSingVal(l_vec(2)+1);
% error_lower_bound_3 = HankSingVal(l_vec(3)+1);
% error_lower_bound_4 = HankSingVal(l_vec(4)+1);

% Parameters:
nb_points = 1000;
f1 = -100; f2 = 100;
freq = linspace(f1, f2, nb_points);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Computation/plot 'matrix-norm-error' across frequencies %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error_matr_1 = ErrorTruncation(A, B, C, A_bt_1, B_bt_1, C_bt_1, freq, Nx);
error_matr_2 = ErrorTruncation(A, B, C, A_bt_2, B_bt_2, C_bt_2, freq, Nx);
error_matr_3 = ErrorTruncation(A, B, C, A_bt_3, B_bt_3, C_bt_3, freq, Nx);
%error_matr_4 = ErrorTruncation(A, B, C, A_bt_4, B_bt_4, C_bt_4, freq, Nx);

% also:
% matr = freqresp(E, freq);
% error_matr = zeros(1, nb_points);
% for k=1:nb_points
%     error_matr(k) = norm(  matr(:, :, k)   );
% end
% figure;
% plot(freq, error_matr ); hold on; % plot
% line([f1, f2], [error_upper_bound, error_upper_bound]); % plot
% line([f1, f2], [error_approx, error_approx]); % plot

figure; % Plot:
plot(freq, error_matr_1, '-.', 'LineWidth', width); 
hold on; 
plot(freq, error_matr_2, '--', 'LineWidth', width);
plot(freq, error_matr_3, ':', 'LineWidth', width);
legend('ell=18', 'ell=20', 'ell=22');
%plot([f1, f2], [error_upper_bound_3, error_upper_bound_3], 'r--');
%plot([f1, f2], [error_approx_3, error_approx_3], 'b--'); 
%legend('error BT', 'upper bound error H infty', 'error H infty');
title('Matrix-norm error'); xlabel('frequency w'); ylabel('||G(iw)-G_{ell}(iw)||');
set(gca,'fontsize',font)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%% Norm Transfer functions %%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% norm_matr_S = NormTransferFunction(A, B, C, D, freq);
% norm_matr_1 = NormTransferFunction(A_bt_1, B_bt_1, C_bt_1, D_bt_1, freq);
% norm_matr_2 = NormTransferFunction(A_bt_2, B_bt_2, C_bt_2, D_bt_2, freq);
% norm_matr_3 = NormTransferFunction(A_bt_3, B_bt_3, C_bt_3, D_bt_3, freq);
% %norm_matr_4 = NormTransferFunction(A_bt_4, B_bt_4, C_bt_4, D_bt_4, freq);
% figure; 
% plot(freq, norm_matr_S, 'k-', 'LineWidth', 1.5); 
% hold on; 
% plot(freq, norm_matr_1, '-.', 'LineWidth', 1.5); 
% plot(freq, norm_matr_2, '--', 'LineWidth', 1.5);
% plot(freq, norm_matr_3, ':', 'LineWidth', 1.5);
% legend('G(iw)','G_7(iw)', 'G_8(iw)', 'G_9(1w)');
% xlabel('value of ell');
% ylabel('|| G_{ell}(iw) ||');
% title('Transfer functions of truncated systems'); xlabel('frequency w');
% set(gca, 'Fontsize', 11);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Error over ell %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l_val = 15:35;
[lower, error, upper] = InfinityError(A, B, C, D, Nx, l_val);
figure;
plot(l_val, lower, 'k--', 'LineWidth', width); 
hold on;
plot(l_val, error, 'k-', 'LineWidth', width);
plot(l_val, upper, 'k:', 'LineWidth', width);
xlabel('value of ell'); ylabel('|| G - G_{ell} ||_{Hinfinity}');
legend('lower bound Hinfinity error', 'Hinfinity error', 'upper bound Hinfinity error');
title('H infinity error');
%set(gca, 'Fontsize', 11);
set(gca,'fontsize',font)