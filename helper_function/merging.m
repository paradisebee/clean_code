function im = merging(Roi, sz, ovl, rge, dim)

if length(dim) == 2
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
    
elseif length(dim) == 3
    if length(sz) == 1
        sz(2) = sz(1);
        sz(3) = sz(2);
    elseif length(sz) == 3

    else
        error('Wrong dimensions of window size!');
    end
else
    error('Wrong dimensions of input image!');
end


im = zeros(dim);
% 2D image
if length(dim) == 2
    r = rge{1};
    c = rge{2};
    n = 1;
    for j = 1:length(c)
        for i = 1:length(r)
            if r(i) == 1
                r_range = 1:sz(1)-ovl(1)/2;
            elseif i == length(r)
                r_range = r(i)+ovl(1)/2:dim(1);
            else
                r_range = r(i)+ovl(1)/2:r(i)+sz(1)-ovl(1)/2-1;
            end
            if c(j) == 1
                c_range = 1:sz(2)-ovl(2)/2;
            elseif j == length(c)
                c_range = c(j)+ovl(2)/2:dim(2);
            else
                c_range = c(j)+ovl(2)/2:c(j)+sz(2)-ovl(2)/2-1;
            end

            tmp = Roi{n};
            if r(i) == 1
                rr = 1:sz(1)-ovl(1)/2;
            elseif i == length(r)
                rr = ovl(1)/2+1:size(tmp,1);
            else
                rr = ovl(1)/2+1:sz(1)-ovl(1)/2;
            end
            if c(j) == 1
                cc = 1:sz(2)-ovl(2)/2;
            elseif j == length(c)
                cc = ovl(2)/2+1:size(tmp,2);
            else
                cc = ovl(2)/2+1:sz(2)-ovl(2)/2;
            end
            
            
            %% important! different between intensity and black & white images
            disp('!!! Modified the codes for intensity or black&white image!!!');
            n = n+1;
            if isempty(tmp)
                continue;
            end
            roi = im(r_range,c_range);
%             roi(roi<tmp) = tmp(roi<tmp);
            im(r_range, c_range) = double(roi | tmp(rr,cc));
            if n > length(Roi)
                break;
            end
        end
    end
end