%% 1-step Adams–Bashforth methods
tic
y0_possible = 0:0.0001:1;
n = 1; % Change order of reaction
A = 4.31^2; % Thiele modulus ^ 2
for i = 1:length(y0_possible)
    fun = @(t,y) dydtf(t,y,n,A); % The second last argument corresponds to order n and last argument corresponds to the constant value.
    t = linspace(0,1,1001);
    t(1) = 1e-8;
    y0 = [y0_possible(i),0];
    [t,y] = Euler(fun,t,y0);
    conc_surface = y(length(y),1);
    %disp(conc_surface)
    if 1 - conc_surface < 1e-3 
        disp('Optimal value found')
        break
    end
end

t_actual = linspace(0,1,1001);
t_actual(1) = 1e-8;

if n == 1
    y_actual = sinh(sqrt(A)*t_actual)./(t_actual.*sinh(sqrt(A)));
end

if n == 5
    y_actual = sqrt(2)*(1 + t_actual.^2 + (sqrt(1 + (4/3)*A)*(1 - t_actual.^2))).^(-1/2);
end

error = norm(y(:,1)' - y_actual);
%disp(error)

hold on
if n == 1 || n == 5
    plot(t_actual,y_actual,'g-','LineWidth',4)
end
plot(t,y(:,1),'b:','LineWidth',4)
legend('Actual solution','Numerical simulation','Location','northwest')
xlabel('\Upsilon')
ylabel('\phi')
title(['Plot of \phi against \Upsilon (Adam-Bashforth order 1 for n = ', num2str(n), ')']) 	
toc
%% 2-step Adams–Bashforth methods
tic
y0_possible = 0:0.0001:1;
n = 1; % Change order of reaction
A = 4.31^2; % Thiele modulus ^ 2
for i = 1:length(y0_possible)
    fun = @(t,y) dydtf(t,y,n,A); % The second last argument corresponds to order n and last argument corresponds to the constant value.
    t = linspace(0,1,1001);
    t(1) = 1e-8;
    y0 = [y0_possible(i),0];
    [t,y] = two_step(fun,t,y0);
    conc_surface = y(length(y),1);
    %disp(conc_surface)
    if 1 - conc_surface < 1e-3 
        disp('Optimal value found')
        break
    end
end

t_actual = linspace(0,1,1001);
t_actual(1) = 1e-8;

if n == 1
    y_actual = sinh(sqrt(A)*t_actual)./(t_actual.*sinh(sqrt(A)));
end

if n == 5
    y_actual = sqrt(2)*(1 + t_actual.^2 + (sqrt(1 + (4/3)*A)*(1 - t_actual.^2))).^(-1/2);
end

error = norm(y(:,1)' - y_actual);
%disp(error)

hold on
if n == 1 || n == 5
    plot(t_actual,y_actual,'g-','LineWidth',4)
end
plot(t,y(:,1),'b:','LineWidth',4)
legend('Actual solution','Numerical simulation','Location','northwest')
xlabel('\Upsilon')
ylabel('\phi')
title(['Plot of \phi against \Upsilon (Adam-Bashforth order 2 for n = ', num2str(n), ')']) 	
toc
%% 3-step Adams–Bashforth methods
tic
y0_possible = 0:0.0001:1;
n = 1; % Change order of reaction
A = 4.31^2; % Thiele modulus ^ 2
for i = 1:length(y0_possible)
    fun = @(t,y) dydtf(t,y,n,A); % The second last argument corresponds to order n and last argument corresponds to the constant value.
    t = linspace(0,1,1001);
    t(1) = 1e-8;
    y0 = [y0_possible(i),0];
    [t,y] = three_step(fun,t,y0);
    conc_surface = y(length(y),1);
    %disp(conc_surface)
    if 1 - conc_surface < 1e-3 
        disp('Optimal value found')
        break
    end
end

t_actual = linspace(0,1,1001);
t_actual(1) = 1e-8;

if n == 1
    y_actual = sinh(sqrt(A)*t_actual)./(t_actual.*sinh(sqrt(A)));
end

if n == 5
    y_actual = sqrt(2)*(1 + t_actual.^2 + (sqrt(1 + (4/3)*A)*(1 - t_actual.^2))).^(-1/2);
end

error = norm(y(:,1)' - y_actual);
%disp(error)

hold on
if n == 1 || n == 5
    plot(t_actual,y_actual,'g-','LineWidth',4)
end
plot(t,y(:,1),'b:','LineWidth',4)
legend('Actual solution','Numerical simulation','Location','northwest')
xlabel('\Upsilon')
ylabel('\phi')
title(['Plot of \phi against \Upsilon (Adam-Bashforth order 3 for n = ', num2str(n), ')']) 	
toc
%% 4-step Adams–Bashforth methods
tic
y0_possible = 0:0.0001:1;
n = 1; % Change order of reaction
A = 4.31^2; % Thiele modulus ^ 2
for i = 1:length(y0_possible)
    fun = @(t,y) dydtf(t,y,n,A); % The second last argument corresponds to order n and last argument corresponds to the constant value.
    t = linspace(0,1,1001);
    t(1) = 1e-8;
    y0 = [y0_possible(i),0];
    [t,y] = four_step(fun,t,y0);
    conc_surface = y(length(y),1);
    %disp(conc_surface)
    if 1 - conc_surface < 1e-3 
        disp('Optimal value found')
        break
    end
end

t_actual = linspace(0,1,1001);
t_actual(1) = 1e-8;

if n == 1
    y_actual = sinh(sqrt(A)*t_actual)./(t_actual.*sinh(sqrt(A)));
end

if n == 5
    y_actual = sqrt(2)*(1 + t_actual.^2 + (sqrt(1 + (4/3)*A)*(1 - t_actual.^2))).^(-1/2);
end

error = norm(y(:,1)' - y_actual);
%disp(error)

hold on
if n == 1 || n == 5
    plot(t_actual,y_actual,'g-','LineWidth',4)
end
plot(t,y(:,1),'b:','LineWidth',4)
legend('Actual solution','Numerical simulation','Location','northwest')
xlabel('\Upsilon')
ylabel('\phi')
title(['Plot of \phi against \Upsilon (Adam-Bashforth order 4 for n = ', num2str(n), ')']) 	
toc