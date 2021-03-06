function [t,y] = three_step(dydt,t,y0) % Objective: To find y value at the prescribed t value.

len = length(t);
h = t(2) - t(1); % Equal step size

% Storing first 3 y values
y(1,:) = y0;

dydt1 = dydt(t(1),y(1,:));
y(2,:) = (dydt1 * h) + y(1,:);

dydt1 = dydt(t(2),y(2,:));
dydt2 = dydt(t(1),y(1,:));
y(3,:) = h*((dydt1 * 3/2) - (dydt2 * 1/2)) + y(2,:);

for i = 4:len
    dydt1 = dydt(t(i-1),y(i-1,:));
    dydt2 = dydt(t(i-2),y(i-2,:));
    dydt3 = dydt(t(i-3),y(i-3,:));
    y(i,:) = h*((dydt1 * 23/12) - (dydt2 * 16/12) + (dydt3 * 5/12)) + y(i-1,:);
end
