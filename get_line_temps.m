function [line_temps] = get_line_temps(L, sigma)
%% Generate templates for line detectors
bmask = zeros(L);
if sigma == 0 
    bmask(ceil(L/2), :) = 1; 
else
    bmask(ceil(L/2)-sigma:ceil(L/2)+sigma, :) = 1; 
end

%% rotate by multiples of 15 degrees
line_temps{1} = bmask;
for k = 1:11
    tmask = imrotate(bmask, k*15, 'bilinear','crop');
    line_temps{k+1} = tmask;
end
    
