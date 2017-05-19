% This program computes the numerical solution to 1st order PDEs, namely
% the transport equation. Two schemes are employed, UPWIND and
% LAX-WENDROFF. The results are presented side by side for comparison.


% In order to replicate the results in the written report, please run the
% program twice, one with dt = 0.033, and once with dt = 0.005.


syms x;                                 % create symbolic variable.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dx = .1;                                % spatial step size.
c = -3;                                 % wave speed.
cfl = dx/c;                             % determines max allowable time step to meet cfl cond.

disp(['By CFL condition, the maximum allowable time step = ' num2str(abs(cfl)) '.']);

dt = input(['Please enter a value for dt such that  0 < dt < ' num2str(abs(cfl)) ':']);

% dt = .033;                              % time step.


dxv = (-10:dx:10)';                     % vector of spatial steps on (-10, 10).
n = length(dxv);                        % matrix index.
sigma = c*(dt/dx);                      % CFL cond. constant.
Tsteps = 100;                           % max. number of time steps.
mytitle = {('UPWIND scheme.'); ['Wave speed c = ' num2str(c) '.'];...
    ['dx = ' num2str(dx) ', ' 'dt = ' num2str(dt) '.']};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates function values for U(t,x) @ time = 0, for all x in [-10, 10].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U_0_x = zeros(1,n)';                     % vector of u(0,x) sols.
fx = 1/(1+x^2);                          % u(0,x).

for i=1:n
    U_0_x(i) = subs(fx, x, dxv(i));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward difference (A) and UPWIND (A_hat) coefficients matrices.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = zeros(n,n);
A(1,1) = sigma+1;
A(1,2) = -sigma;
A(n,n) = sigma+1;
for i=2:n-1
    A(i,i) =sigma+1;
    A(i,i+1) = -sigma;
end

A_hat = zeros(n,n);
A_hat(1,1) = -sigma+1;
A_hat(n,n) = -sigma+1;
A_hat(n,n-1) = sigma;
for i=2:n-1
    A_hat(i,i) = -sigma+1;
    A_hat(i,i-1) = sigma;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LAX WENDROFF coefficients matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LW = zeros(n,n);
LW(1,1) = -(sigma*sigma - 1);
LW(1,2) = 0.5*sigma*(sigma-1);
LW(n,n) = -(sigma*sigma - 1);
LW(n,n-1) = 0.5*sigma*(sigma+1);

for i=2:n-1
    LW(i,i) = -(sigma*sigma - 1);
    LW(i,i+1) = 0.5*sigma*(sigma-1);
    LW(i, i-1) = 0.5*sigma*(sigma+1);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates function values for U(t,x) @ time = dt, for all x in [-10, 10]
% according to wave speed (pos. or neg.)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U_t_x = zeros(n,Tsteps);                    % solutions vector.
LWU_t_x = zeros(n,Tsteps);                  % Lax Wendroff solutions vector.
 
if c <= 0
    U_t_x(:,1) = A*(U_0_x);
    LWU_t_x(:,1) = LW*(U_0_x);
    for i=2:Tsteps
        U_t_x(:,i) = A*(U_t_x(:,i-1));
        LWU_t_x(:,i) = LW*(LWU_t_x(:,i-1));
    end
else
    U_t_x(:,1) = A_hat*(U_0_x);
    LWU_t_x(:,1) = LW*(U_0_x);
    for i=2:Tsteps
        U_t_x(:,i) = A_hat*(U_t_x(:,i-1));
        LWU_t_x(:,i) = LW*(LWU_t_x(:,i-1));
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots and animations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name', 'Numerical solution to the transport equation','Position', [100, 200, 1549, 1000]);

subplot(2,2,4)
contourf(U_t_x)
colorbar
title({('UPWIND scheme.');('Behavior along characteristic lines.')})
view(90, 90)

U_movie = zeros(size(U_t_x));
U_movie(:,1) = U_t_x(:,1);
subplot(2,2,3)
h = mesh(U_movie);
% colormap(jet);
colorbar;
xlim([0 100])
ylim([0 200])
view(95, 32)
title(mytitle)
grid on
M(1) = getframe;

for i=1:Tsteps
    subplot(2,2,1)
    plot(dxv, U_t_x(:,i),'LineWidth',1.5);
    axis([-10, 10, 1.1*min(U_t_x(:)), 1.1*max(U_t_x(:))]);
    title(mytitle)
    grid on
    
    subplot(2,2,2)
    plot(dxv, LWU_t_x(:,i),'LineWidth',1.5);
    axis([-10, 10, 1.1*min(U_t_x(:)), 1.1*max(U_t_x(:))]);
    title({('LAX WENDROFF scheme.'); ['Wave speed c = ' num2str(c) '.'];...
    ['dx = ' num2str(dx) ', ' 'dt = ' num2str(dt) '.']})
    grid on
    
    U_movie(:,i) = U_t_x(:,i);
    set(h,'ZData',U_movie);
    M(i) = getframe;
end


% pause(3.0)
% U_movie = zeros(size(U_t_x));
% U_movie(:,1) = LWU_t_x(:,1);
% subplot(1,2,2)
% h = surf(U_movie);
% colormap(jet);
% colorbar;
% xlim([0 100])
% ylim([0 200])
% zlim([0 1])
% view(95, 32)
% title({('LAX WENDROFF scheme'); ['Wave speed c = ' num2str(c) '.'];...
%     ['dx = ' num2str(dx) ', ' 'dt = ' num2str(dt) '.']})
% grid on
% M(1) = getframe;
% for i=1:Tsteps
%     U_movie(:,i) = LWU_t_x(:,i);
%     set(h,'ZData',U_movie);
%     M(i) = getframe;
% end


    
