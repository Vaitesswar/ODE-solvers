clear
tic
P = 4.31^2; % Thiele modulus ^ 2
n = 5; % Order of reaction
alpha = 0.9; % Relaxation factor

r = linspace(0,1,1001)';
step = r(2) - r(1);

profile = r.^2; % Initial profile

% Create const matrix A (without first and last phi rows)
a1 = 1./step;
a2 = a1.^2;
[m,p] = size(r);
A = zeros(m-2,m-2);
b = zeros(m-2,1);

for i = 2:length(A)-1
    A(i,i-1) = a2  - a1./r(i+1);
    A(i,i) = -2.*a2;
    A(i,i+1) = a2  + a1./r(i+1);
end

% Enforce the first boundary condition (dphi/dr = 0 when r = 0)
A(1,1) = -2.*a2./3 - 4.*a1/(3.*r(2));
A(1,2) = -1 .* A(1,1);

% Enforce the second boundary condition (phi = 1 when r = 1)
A(end,end-1) = a2 - a1./r(end-1);
A(end,end) = -2.*a2;

% b value for the last eq.
b(end,1) = -1*(a2 + a1./r(end-1));

for i = 1:1500
    % Store old values
    old_profile = profile;
    
    % Add non-linear part to diagonal
    A1 = A - diag(P*(old_profile(2:end-1,1).^(n-1)));
    
    % norm0 = ||A*x0 - b||
    norm0 = norm(A1*profile(2:end-1,1) - b);
    
    % new_profile = A\b
    guess_profile = A1\b;
    
    % Update profile using weighted sum of old and new values
    profile(2:end-1,1) = alpha.*old_profile(2:end-1,1) + (1-alpha).*guess_profile;
    
    % Enforce the first boundary condition (dphi/dr = 0 when r = 0)
    profile(1,1) = (4.*profile(2) - profile(3))./3;
    
    % Enforce the second boundary condition (phi = 1 when r = 1)
    profile(end,1) = 1;
    
    % norm1 = ||phi_new - phi_old||/||phi_old||
    norm1 = norm(profile - old_profile)/norm(old_profile);
    
    if min(norm0,norm1) < 1e-8
        disp(i-1); % The solution already converged in the previous iteration.
        disp(min(norm0,norm1));
        break
    end
end

if n == 1
    profile_actual = sinh(sqrt(P)*r)./(r.*sinh(sqrt(P)) + 1e-40);
    profile_actual(1,1) = (4.*profile_actual(2) - profile_actual(3))./3;
end

if n == 5
    profile_actual = sqrt(2)*(1 + r.^2 + (sqrt(1 + (4/3)*P)*(1 - r.^2))).^(-1/2);
end

hold on
if n == 1 || n == 5
    plot(r,profile_actual,'g-','LineWidth',4)
end
plot(r,profile,'b:','LineWidth',4)
legend('Actual solution','Numerical simulation','Location','northwest')
xlabel('\Upsilon')
ylabel('\phi')
title(['Plot of \phi against \Upsilon (Finite difference method for n = ', num2str(n), ')'])
toc
