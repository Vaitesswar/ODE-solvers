function [t,y] =  rk2_gen(dydt,t,y0)

len = length(t);
y(1,:) = y0;
h = t(2) - t(1);

for i = 1:len-1
    dydt1 = dydt(t(i),y(i,:));
    y2 = (dydt1.*h) + y(i,:);
    dydt2 = dydt(t(i+1),y2);
    y(i+1,:) = (dydt1 + dydt2) .* (h/2) + y(i,:);
end
