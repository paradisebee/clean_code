function [x, y, pos] = line2D(ctr, vec, len)

%LINE2D   Get 2D line.
%   LINE2D(CTR, VEC, LEN) takes in a point on the line,
%   and the direction vector of the point, and return a 
%   line perpendicular to the direction vector, with total
%   length equals to (2*len+1).
%
%   Input
%       ctr : point subscript in the format of [y,x] or [r,c]
%       vec : [tx, ty], direction vector of the point
%       len : half length of the line, total length will be 2*len
%
%   Output
%       [x,y] : the output 2D line
%       pos   : the position of the ctr point in the 2D line
%
%

if (nargin ~= 3)
    ERROR('Incorrect number of arguments!\n');
else
    if length(ctr) ~= 2 || length(vec) ~= 2 || length(len) ~= 1
        ERROR('Incorrect parameter format!\n');
    end
end


if abs(vec(1)) <= 10e-4 % vertical line
    x = ctr(2)-len:ctr(2)+len;
    y = ctr(1)*ones(1,length(x));
elseif abs(vec(2)) <= 10e-4 % horizontal line
    y = ctr(1)-len:ctr(1)+len;
    x = ctr(2)*ones(1,length(y));
else
    gy = -vec(1);
    gx = vec(2);
    x1 = ctr(2)-gx*len;
    y1 = ctr(1)-gy*len;
    x2 = ctr(2)+gx*len;
    y2 = ctr(1)+gy*len;
    
    x = x1:gx:x2;
    y = y1:gy:y2;
end

% column vectors
x = x';
y = y';

[~,pos] = min((x-ctr(2)).^2+(y-ctr(1)).^2);
