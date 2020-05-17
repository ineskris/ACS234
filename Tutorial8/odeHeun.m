function [t,y] = odeHeun(f,tspan,y0,h)
% Heun method to solve vector differential equation y’(t) = f(t,y(t))
% For a tspan = [t0,tf], with the initial value y0 and a step-size h.
if nargin <4 || h <= 0, h = 1; end
if nargin <3, y0 = 0; end
N = (tspan(2)-tspan(1))/h;
N = round(N); % The number of steps to be integrated
t = tspan(1)+[0:N]'*h; % A vector of base points
y(1,:) = y0(:)'; % Making the initial value a row vector
for k = 1:N
fk = feval(f,t(k),y(k,:));
y(k+1,:) = y(k,:)+h*fk; % Eq.(5.8)
y(k+1,:) = y(k,:) +h/2*(fk +feval(f,t(k+1),y(k+1,:))); % Eq.(5.25)
end