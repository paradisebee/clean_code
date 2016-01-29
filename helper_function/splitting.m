function [Roi, rge] = splitting(im, sz, ovl, mask)

%% input parameters checking
% if input image is 2D
if length(size(im)) == 2
    if length(sz) == 1
        sz(2) = sz(1);
    elseif length(sz) == 2

    else
        error('Wrong dimensions of window size!');
    end

    if length(ovl) == 1
        ovl(2) = ovl(1);
    elseif length(ovl) == 2

    else
        error('Wrong dimensions of window overlap!');
    end

% if input image is 3D
elseif length(size(im)) == 3
    if length(sz) == 1
        sz(2) = sz(1);
        sz(3) = sz(2);
    elseif length(sz) == 3

    else
        error('Wrong dimensions of window size!');
    end

    if length(ovl) == 1
        ovl(2) = ovl(1);
        ovl(3) = ovl(2);
    elseif length(ovl) == 3

    else
        error('Wrong dimensions of window overlap!');
    end
    
else
    error('Wrong dimensions of input image!');
end
   
dim = size(im);
rge = cell(length(dim)+1,1);
for i = 1:length(dim)
    tmp = (sz(i)-ovl(i)+1):(sz(i)-ovl(i)):(dim(i)-sz(i)+1);
    if tmp(end) < dim(i)-sz(i)+1
        tmp(end+1) = dim(i)-sz(i)+1;
    end
    rge{i} = [1, tmp];
end

% 2D image
if length(dim) == 2
    r = rge{1};
    c = rge{2};
    Roi = cell(length(r)*length(c),1);
    n = 1;
    for j = 1:length(c)
        for i = 1:length(r)
            r_range = r(i):r(i)+sz(1)-1;
            c_range = c(j):c(j)+sz(2)-1;
            roi = im(r_range,c_range);
%              roi_mask = mask(r_range,c_range);
%                 lgc = (roi_mask(:)==1);
%                 roi(lgc) = (roi(lgc)-min(roi(lgc)))./(max(roi(lgc))-min(roi(lgc)));
%                 imwrite(roi, ['roi\roi', num2str(n),'.bmp'], 'bmp');
            Roi{n} = roi;
            n = n+1;
        end
    end
end
% 3D image
if length(dim) == 3
    r = rge{1};
    c = rge{2};
    s = rge{3};
    Roi = cell(length(r)*length(c)*length(s));
    n = 1;
    for j = 1:length(c)
        for i = 1:length(r)
            for k = 1:length(s)
                r_range = r(i):r(i)+sz(1)-1;
                c_range = c(j):c(j)+sz(2)-1;
                s_range = s(k):s(k)+sz(3)-1;
                roi = im(r_range,c_range,s_range);
                Roi{n} = roi;
                n = n+1;
            end
        end
    end
end

