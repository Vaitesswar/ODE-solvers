function dydt = dydtf(t,y,n,A)

dydt(2) = A*(y(1)^n) - (2*y(2)/t);
dydt(1) = y(2);

end

