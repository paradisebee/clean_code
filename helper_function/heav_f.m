function f = heav_f(x, epsilon, h_thres)

if nargin < 2
    epsilon = 0.5;
    h_thres = 0.8;
end
f = 0.5*(1+2/pi*atan((x-h_thres)/epsilon));