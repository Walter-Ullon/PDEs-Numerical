syms x
disp('Please choose a scheme:');
choice = input('Enter 1 for EXPLICIT, 2 for IMPLICIT, 3 for CRANK NICOLSON: ');


l = 1;                                  % length of spatial interval.
gamma = 1;                              % thermal diffusivity constant.
dt = 0.001;                             % delta t -- time step.
dx = 0.05;                              % delta x -- spatial step.
dxv = (0:dx:l)';                        % vector of spatial increments up to l. Determines spatial mesh size.
Tsteps = 100;                           % # of time steps to calculate.
                                        
mu = gamma*dt/(dx*dx);                  % (mesh ratio)*gamma - for finite diff. coeff. matrix.
mu1 = (1-2*mu);                         % see above.
mu2 = (1+2*mu); 
mu3 = (1+mu);                           
mu4 = (1-mu);
n = length(dxv);                        % indexes vectors.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check stability conditions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choice ==1 
    if mu>.5
        disp('The scheme is unstable. Please adjust dt, dx accordingly. Terminating now.')
        return
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates function values for U(t,x) @ time = 0, for all x in [0, l].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fx = [2*abs(x-(1/6))-(1/3), 0, 1/2-3*abs(x-(5/6))];                % f(x) interval functions.
fx = [0, 1, 0];
bounds = [0, 1/3, 2/3, 1];                                         % integral bounds.
U_dt_x = zeros(n,1);
for i=1:n
    if dxv(i)<bounds(2)
        U_dt_x(i) = subs(fx(1),x,dxv(i));
    elseif dxv(i)>=bounds(2) && dxv(i)<bounds(3)
        U_dt_x(i) = subs(fx(2),x,dxv(i));
    else
        U_dt_x(i) = subs(fx(3),x,dxv(i));
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finite diff. coefficients matrix --> U(t+1, x[0,l]) = A*U(t, x[0,l]) + B.
% Vector 'b' contains the Dirilecht boundary nodes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = zeros(n, 1);
b(1) = mu*U_dt_x(1);
b(n) = mu*U_dt_x(n);


if choice == 1                      % enter explicit scheme.
    mytitle = 'Explicit Scheme';

    A = zeros(n,n);
    A(1,1) = mu1;
    A(1,2) = mu;
    A(n, n-1) = mu;
    A(n,n) = mu1;
    for i=2:n-1
        A(i,i) = mu1;
        A(i, i-1) = mu;
        A(i, i+1) = mu;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculates the next time step by performing matrix multiplication, i.e.
    % the "explicit scheme": --> U(t+1, x[0,l]) = A*U(t, x[0,l]) + b.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Utx = zeros(n,Tsteps);
    Utx(:,1) = A*U_dt_x + b;
    for i=2:Tsteps
        Utx(:,i) = A*Utx(:,i-1) + b;
    end

elseif choice ==2                   % enter implicit scheme.
    mytitle = 'Implicit Scheme';
    
    A_hat = zeros(n,n);
    A_hat(1,1) = mu2;
    A_hat(1,2) = -mu;
    A_hat(n, n-1) = -mu;
    A_hat(n,n) = mu2;
    for i=2:n-1
      A_hat(i,i) = mu2;
      A_hat(i, i-1) = -mu;
      A_hat(i, i+1) = -mu;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculates the next time step by performing matrix multiplication, i.e.
    % the "implicit scheme": --> U(t+1, x[0,l]) = A_hat*(U(t, x[0,l]) + b).
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Utx = zeros(n,Tsteps);
    Utx(:,1) = inv(A_hat)*(U_dt_x + b);
    for i=2:Tsteps
        Utx(:,i) = inv(A_hat)*(Utx(:,i-1) + b);
    end

else                                   % enter Crank Nicolson
    mytitle = 'Crank Nicolson Scheme';
    
    B_hat = zeros(n,n);
    B_hat(1,1) = mu3;
    B_hat(1,2) = -0.5*mu;
    B_hat(n, n-1) = -0.5*mu;    
    B_hat(n,n) = mu3;
    for i=2:n-1
        B_hat(i,i) = mu3;
        B_hat(i, i-1) = -0.5*mu;
        B_hat(i, i+1) = -0.5*mu;
    end

    B = zeros(n,n);
    B(1,1) = mu4;
    B(1,2) = 0.5*mu;
    B(n, n-1) = 0.5*mu;
    B(n,n) = mu4;
    for i=2:n-1
        B(i,i) = mu4;
        B(i, i-1) = 0.5*mu;
        B(i, i+1) = 0.5*mu;
    end

    b = zeros(n, 1);
    b(1) = mu*U_dt_x(1);
    b(n) = mu*U_dt_x(n);

    Utx = zeros(n,Tsteps);
    Utx(:,1) = inv(B_hat)*(B*U_dt_x + 0.5*(b+b));
    for i=2:Tsteps
        Utx(:,i) = inv(B_hat)*(B*Utx(:,i-1) + 0.5*(b+b));
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots ans animations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:Tsteps
    plot(Utx(:,i),'LineWidth',1.5);
    axis([0, n, 1.1*min(Utx(:)), 1.1*max(Utx(:))]);
    title(mytitle)
    grid on
    M(i) = getframe;
end
