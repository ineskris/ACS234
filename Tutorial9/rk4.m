function [t,y] = rk4 (f, tspan, y0, h)
% this is code for RK4 according to the formulars from tutorial
if nargin <4 || h <= 0, h = 1; end
if nargin <3, y0 = 0; end
N = (tspan(2)-tspan(1))/h;
N = round(N); % The number of steps to be integrated
t = tspan(1)+(0:N)'*h; % A vector of base points
y(1,:) = y0(:)'; % Making the initial value a ROW vector
for k = 1:N
    k1 = feval(f,t(k),y(k,:));
    k2 = feval(f,t(k)+h/2,y(k,:)+h*k1/2);
    k3 = feval(f,t(k)+h/2,y(k,:)+h*k2/2);
    k4 = feval(f,t(k)+h,y(k,:)+h*k3);
    y(k + 1,:) = y(k,:) + h*(k1+2*k2+2*k3+k4)/6;
end