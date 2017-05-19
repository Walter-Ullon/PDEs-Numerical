syms x
format long

l = 3;                                  % length of spatial interval.
t = 0;                                  % time.
c = 8;                                  % wave speed

dx = .1;                                % delta x -- spatial step
dt = .1;                              % delta t -- time step.
dxv = (0:dx:l)';                        % vector of spatial increments up to l
n = length(dxv);
sigma = c*(dt/dx);                      % computes cfl condition
Tsteps = 10^3;                          % # of iterations

f(x) = 1-2*abs(x-1);                    % initial cond. function

U_0_x = zeros(n,Tsteps);                % vector of initial conditions

for i=1:length(dxv)
    if dxv(i)<1/2 || dxv(i)>3/2
        U_0_x(i) = 0;
    else
        U_0_x(i) = f(dxv(i));
    end
end


b = zeros(n,1);                         % boundary conditions vector
b(1) = 0;
b(n) = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matrix of coefficients 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B = zeros(n,n);
B(1,1) = 2*(1-sigma^2);
B(1,2) = sigma^2;
B(n,n) = 2*(1-sigma^2);
B(n,n-1) = sigma^2;

for i=2:n-1
    B(i,i-1) = sigma^2;
    B(i,i) = 2*(1-sigma^2);
    B(i,i+1) = sigma^2;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution u(t,x) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U_t_x = zeros(n,Tsteps);
U_t_x(:,1) = U_0_x(:,1);
U_t_x(:,2) = 0.5*B*U_t_x(:,1);
for i=3:Tsteps
    U_t_x(:,i) = B*U_t_x(:,i-1) - U_t_x(:,i-2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots and animations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:Tsteps
    plot(U_t_x(:,i),'LineWidth',1.5);
    axis([0,n, 1.1*min(U_t_x(:)), 1.1*max(U_t_x(:))]);
    grid on
    M(i) = getframe;
end


