%% Initialization
clear;clc;

% Parameters
U = 100;
eta_init = 1;   % inital anomalous expectation value eta:=<b>
d = 3;          % number of dimensions
n_max = 100;    % maximal number of particles on the site
                % (matrices will have dimension (n_max+1)*(n_max+1)
t_over_U_max = 0.05;
mu_over_U_max = 4;

% Operators
bdagger = diag((1:n_max).^(1/2),-1);    % creation op.
b = transpose(bdagger);                 % annihilation op.
n = diag(0:n_max);                      % number op.


% Self-Consistency Procedure

c = (t_over_U_max+1)*(mu_over_U_max+1);
results = zeros(c,4);
k = 1;
t_intv = 0:0.0005:t_over_U_max;
mu_intv = 0:0.05:mu_over_U_max;
p_max = 0;
for t_over_U=t_intv
for mu_over_U=mu_intv
    eta = eta_init;
    eta_prev = Inf;
    p=0;
    while norm(eta_prev - eta)>1e-4
        % Calculate Hamiltonian in units of U
        H = -mu_over_U*n + 1/2*n.*(n-1)-2*d*t_over_U*eta*bdagger - 2*d*t_over_U*conj(eta)*b;
        % Diagonalize
        [V,D] = eig(H);
        E = diag(D);                        % (eigen-)energies
        [E0,idx] = min(E);                  % ground state energy
        gs = V(:,idx);                      % ground state

        % Update anomalous expectation value and see if it has changed
        % (self-consistency)
        eta_prev = eta;
        eta = gs' * b * gs;
        p=p+1;
    end;
    if p_max < p:
        p_max = p;
    end;
    results(k,1) = t_over_U;
    results(k,2) = mu_over_U;
    results(k,3) = eta;
    results(k,4) = gs' * n * gs;
    k=k+1;
end;
end;

%% Plotting
[X,Y] = meshgrid(t_intv,mu_intv);
z = results(:,3);
Z = vec2mat(z,size(X,1))';

hold on;
pcolor(X,Y,Z);
colormap(gray);
shading interp;
xlabel('t/U');
ylabel('mu/U');
title('Bose-Hubbard T=0 Phasediagram in Mean Field Approx.')
set(gca,'XTick',[0:0.01:t_over_U_max]);
set(gca,'YTick',[0:1:mu_over_U_max]);
set(gca,'box','on');
%colorbar;
print('phasediagram', '-dpng', '-r300');