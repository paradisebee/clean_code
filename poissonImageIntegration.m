function output = poissonImageIntegration(LAP, DEST, mask, tx, ty, ww, hh)

% Ic = imcrop(I1,[min(c) min(r) maxW maxH]);  % crop the image in the RIO
% % crop mask - make the mask RGB (3 layers)
% Mc = zeros(size(Ic));                       % make a copy of Ic
% 
% % get laplacian of the source image
% H =fspecial('laplacian',0);
% LAP = imfilter(double(I1),H,'same');
% Id = imcrop(LAP,[min(c) min(r) maxW maxH]);
% 
% % multiple the Mask by the Image to get only the pixels in the RIO
% LAP = immultiply(Id,Mc);
% 
% 
% [H,W,C] = size(LAP);


% find DEST boundary points
se = strel('disk',1);
outer_bp = imdilate(mask,se) - mask;

fd_bp = DEST.*outer_bp;   % f* boundary points (DEST)


% only process roi with mask
roi_mask = mask(ty-1:ty+hh,tx-1:tx+ww);     % bw roi mask
I = DEST;
dest_roi = I(ty-1:ty+hh,tx-1:tx+ww);    % destination roi
f_idx = find(roi_mask==1);
[Row,Col] = ind2sub(size(roi_mask),f_idx);

roi_bp = fd_bp(ty-1:ty+hh,tx-1:tx+ww);  % destination boundary intensity
roi_outer_bp = outer_bp(ty-1:ty+hh,tx-1:tx+ww); % roi outer boundary
roi_lap = LAP(ty-1:ty+hh,tx-1:tx+ww);   % source laplacian roi

f = zeros(size(roi_bp));


dest_bp = roi_bp;
% stupid modification
% dest_bp(dest_bp~=0) = dest_bp(dest_bp~=0)-20;
%
lap = roi_lap;
indexArray1 = [];
indexArray2 = [];
valueArray = [];
b = zeros(length(f_idx),1);
% create linear system using sparse matrix
for i = 1:length(f_idx)
    n = 0;
    r = Row(i);
    c = Col(i);
    b(i) = b(i)-lap(r,c);
    if roi_mask(r+1,c)
        n = n+1;
        indexArray1(end+1) = i;
        indexArray2(end+1) = find(f_idx==(f_idx(i)+1));
        valueArray(end+1) = -1;
    else
        if roi_outer_bp(r+1,c)
            n = n+1;
            b(i) = b(i)+dest_bp(r+1,c);
        end
    end
    if roi_mask(r-1,c),
        n = n+1;
        indexArray1(end+1) = i;
        indexArray2(end+1) = find(f_idx==(f_idx(i)-1));
        valueArray(end+1) = -1;
    else
        if roi_outer_bp(r-1,c)
            n = n+1;
            b(i) = b(i)+dest_bp(r-1,c);
        end
    end
    if roi_mask(r,c+1),
        n = n+1;
        indexArray1(end+1) = i;
        indexArray2(end+1) = find(f_idx==(f_idx(i)+size(roi_mask,1)));
        valueArray(end+1) = -1;
    else
        if roi_outer_bp(r,c+1)
            n = n+1;
            b(i) = b(i)+dest_bp(r,c+1);
        end
    end
    if roi_mask(r,c-1),
        n = n+1;
        indexArray1(end+1) = i;
        indexArray2(end+1) = find(f_idx==(f_idx(i)-size(roi_mask,1)));
        valueArray(end+1) = -1;
    else
        if roi_outer_bp(r,c-1)
            n = n+1;
            b(i) = b(i)+dest_bp(r,c-1);
        end
    end
    if n
        indexArray1(end+1) = i;
        indexArray2(end+1) = i;
        valueArray(end+1) = n;
    end
end

A = sparse(indexArray1,indexArray2,valueArray,length(f_idx),length(f_idx));
x = ones(length(f_idx),1);

% jacobi method
x_old = x;
x_new = x_old;
th0 = 0;
while (1)
    sigma = A*x_old-diag(diag(A))*x_old;
    x_new = (b-sigma)./full(diag(A));
    th1 = max(abs(x_new - x_old))           % new threshold
    if max(abs(x_new - x_old)) < 0.01
        break;
    else
        th0 = max(abs(x_new - x_old));      % old threshold
        x_old = x_new;
    end
end

result = zeros(size(roi_mask));
result(f_idx) = x_new;

f = result;

tmp = dest_roi;
tmp(f_idx) = x_new;
dest_roi = tmp;


I(ty-1:ty+hh,tx-1:tx+ww) = dest_roi;    % change dest roi to calculated roi


output = I;

return
