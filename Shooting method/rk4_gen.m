function [t,y] = rk4_gen(dydt,t,y0) % Same as 1 ODE rk4 except that dydt and y0 are row vectors.

len = length(t);
y(1,:) = y0;
h = t(2) - t(1);

for i = 1:len-1
    % Point 1
    y1 = y(i,:);
    dydt1 = dydt(t(i),y1); % Finding the dydt for the first point.
    
    % Point 2
    t_int = (t(i) + t(i+1))/2;
    y2_1 = (dydt1.*(h/2)) + y1; % Performing Euler method (h/2 is used since 1/3 simpson's rule)
    dydt2_1 = dydt(t_int,y2_1);  
    y2_2 = (dydt2_1.*(h/2)) + y1; % Performing Euler method 
    dydt2_2 = dydt(t_int,y2_2);
    dydt2 = (dydt2_1 + dydt2_2)/2; % Finding the average dydt for second point.
    
    % Point 3
    y3 = y1 + (dydt2_2.*h); % Performing Euler method (h since this is the interval betwen 3 data points)
    dydt3 = dydt(t(i+1),y3); % Finding the dydt for the third point.
    
    % Performing 1/3 Simpson's rule
    sum = (dydt1 + dydt2*4 + dydt3) .* h./6;
    y(i+1,:) = y1 + sum;
end
    
