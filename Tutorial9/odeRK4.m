function [t,y] = odeRK4(f,tspan,y0,h)
% Runge-Kutta method to solve vector differential eqn y’(t) = f(t,y(t))
% for tspan = [t0,tf] and with the initial value y0 and a step-size h.
if nargin < 4 || h <= 0, h = 1; end
if nargin < 3, y0 = 0; end
y(1,:) = y0(:)'; % Making it a row vector
N = (tspan(2) - tspan(1))/h;
N = round(N);
t = tspan(1)+[0:N]'*h;
for n = 1:N
f1 = h*feval(f,t(k),y(k,:)); f1 = f1'; % Eq.(5.32a)
f2 = h*feval(f,t(k) + h/2,y(k,:) + f1/2); f2 = f2'; % Eq.(5.32b)
f3 = h*feval(f,t(k) + h/2,y(k,:) + f2/2); f3 = f3'; % Eq.(5.32c)
f4 = h*feval(f,t(k) + h,y(k,:) + f3); f4 = f4'; % Eq.(5.32d)
y(k + 1,:) = y(k,:) + (f1 + 2*(f2 + f3) + f4)/6; % Eq.(5.31)
end