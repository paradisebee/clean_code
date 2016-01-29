function curvature = curv(x,y)

addpath C:\Users\Bichao\Desktop\clean_code\helper_function\Dgradient

x = x(:); y = y(:);

%% ***************** no use !! *************************
% % % manually extend the two ends of the curve to reduce error
% ax = x(2)-x(1);
% bx = x(end)-x(end-1);
% x = [x(1)-2*ax; x(1)-ax; x; x(end)+bx; x(end)+2*bx];
% ay = y(2)-y(1);
% by = y(end)-y(end-1);
% y = [y(1)-2*ay; y(1)-ay; y; y(end)+by; y(end)+2*by];
x=smooth(x,'lowess');
y=smooth(y,'lowess');


dx  = DGradient(x,1,[],'2ndOrder');
ddx = DGradient(dx,1,[],'2ndOrder');
dy  = DGradient(y,1,[],'2ndOrder');
ddy = DGradient(dy,1,[],'2ndOrder');
num   = dx .* ddy - ddx .* dy;
denom = dx .* dx + dy .* dy;
denom = sqrt(denom);
denom = denom .* denom .* denom;
curvature = abs(num) ./ denom;
curvature(denom <= 0) = NaN;


% curvature = curvature(3:end-2);