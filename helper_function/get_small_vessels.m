function [im_obj, im_bkg] = get_small_vessels(im_in, im_mask)

%%
addpath /home/bee/Dropbox/Code/Download/FastMarching/functions

% get outer boundary points
se = strel('disk',1);
% v_boundary = v-imerode(v,se);
v_boundary = imdilate(im_in,se)-im_in;

% source points
sp = find(v_boundary);
[sp_r, sp_c] = ind2sub(size(im_in),sp);
boundaryPoints = [sp_r'; sp_c'];

% initial speed
F1 = im_in;
F1(F1==0 | ~im_in) = eps;

% fast marching for signed distance map
D1 = msfm2d(double(F1), boundaryPoints, true, true);

D1 = D1.*im_in;
D1(~im_in)=inf;

%%
D1(sp) = 0;
D1(~im_in) = 0;
d1 = D1.*im_in;
d1(d1<=3) = 0;

%%
SP = d1;
SP(SP~=0) = 1;
% figure
% imshow(SP,[])
for i = 1:10
    SP = imdilate(SP,se);
end

%%
bw = zeros(size(im_in));
cc = bwconncomp(SP);
numPixels = cellfun(@numel,cc.PixelIdxList);
[~,idx] = max(numPixels);
bw(cc.PixelIdxList{idx}) = 1;

sg = ~bw & im_in;
sg_boundary = sg-imerode(sg,se);

im_obj = sg;
im_bkg = ~(sg|bw).*im_mask;