function [img, img_raw, pmask, bigimg] = vessel_preproc(img, cmask)

%% input: normalized color fundus image (0, 1)
if size(img,3) == 3
    R = img(:,:,1);
    G = img(:,:,2);
    B = img(:,:,3);
else
    G = img;
end

img = G;
[bigimg img_raw] = getBigimg(G, cmask);
img = adapthisteq(img_raw);
bigimg = adapthisteq(bigimg);
pmask = img>=0.02;

function [bigimg smallimg] = getBigimg(img,mask)

[sizey, sizex] = size(img);

bigimg = zeros(sizey + 100, sizex + 100);
bigimg(51:(50+sizey), 51:(50+sizex)) = img;

bigmask = logical(zeros(sizey + 100, sizex + 100));
bigmask(51:(50+sizey), (51:50+sizex)) = mask;

% Creates artificial extension of image.
bigimg = fakepad(bigimg, bigmask, 5, 20);
smallimg = bigimg(51:(50+sizey), 51:(50+sizex));



