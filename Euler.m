function [t,y] = Euler(dydt,t,y0) % Objective: To find y value at the prescribed t value.

len = length(t);
h = t(2) - t(1); % Equal step size
y(1,:) = y0;

for i = 1:len-1 % No need to perform these operations for the last t point.
    dydt1 = dydt(t(i),y(i,:));
    y(i+1,:) = (dydt1 * h) + y(i,:);
end


end