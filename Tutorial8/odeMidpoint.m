function [t,y] = odeMidpoint(f, tspan, y0, h)
% Euler’s method to solve vector differential equation y’(t) = f(t, y(t))
% For a time tspan = [t0, tf], the initial value y0, and a step-size h.
%
if nargin <4 || h <= 0, h = 1; end
if nargin <3, y0 = 0; end
%
N = (tspan(2)-tspan(1))/h;
N = round(N); % The number of steps to be integrated
t = tspan(1)+(0:N)'*h; % A vector of base points
y(1,:) = y0(:)'; % Making the initial value a ROW vector
for k = 1:N
    y(k + 1,:) = y(k,:) + h/2*feval(f,t(k),y(k,:));
    t(k) = t(k) + 1/2;
    y(k + 1,:) = y(k,:) +h*feval(f, t(k), y(k+1,:)); % Eq.(5.8)
end