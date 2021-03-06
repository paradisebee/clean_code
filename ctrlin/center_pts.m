function [new_ctr, rads] = center_pts(im, ctrp, ctr_tx, ctr_ty, width, h_r)

% get cross-sectional pixels
[x, y, pos] = line2D(ctrp, [ctr_tx, ctr_ty], width);

% get local roi
r_range = min(floor(y))-1:max(ceil(y))+1;
c_range = min(floor(x))-1:max(ceil(x))+1;

y = y-min(r_range)+1;
x = x-min(c_range)+1;

roi = im(r_range,c_range);
[X,Y] = meshgrid(1:size(roi,2),1:size(roi,1));
% get intensity profile
profile = interp2(X,Y,roi,x,y);
[xmax,imax,xmin,imin] = extrema(profile);
[~, mi] = min(abs(pos-imax));
pos = imax(mi);
% find double-edge positions
rside = profile(pos:end-1)-profile(pos+1:end);
lside = profile(pos:-1:2)-profile(pos-1:-1:1);
% dx = diff(profile);
% dxx = diff(dx);
% dxx = smooth(dxx);
addpath F:\Dropbox\Code\Download\extrema
[xmax,imax,xmin,imin] = extrema(rside);
r_edge = pos+min(imax);
[xmax,imax,xmin,imin] = extrema(lside);
l_edge = pos-min(imax);

le = [y(l_edge),x(l_edge)]; re = [y(r_edge),x(r_edge)]; 
new_ctr = [(le+re)/2, (le+re)/2];
new_ctr = [new_ctr(1)+min(r_range)-1, new_ctr(2)+min(c_range)-1];
rads = sqrt((y(l_edge)-y(r_edge))^2+(x(l_edge)-x(r_edge))^2)/2;


% radius check
if h_r ~= 0 
    if rads-h_r>2
        if r_edge-pos > pos-l_edge
            new_ctr = new_ctr + (le-new_ctr).*(rads-h_r)./rads;            
        else
            new_ctr = new_ctr + (re-new_ctr).*(rads-h_r)./rads;
        end
        rads = h_r;
    end
end