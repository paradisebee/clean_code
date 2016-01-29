function [tx, ty, gmag_norm] = get_tangent(im)

%% gradient
sh = fspecial('sobel');
sv = sh';
gy = filter2(sh, im, 'same');
gx = filter2(sv, im, 'same');

% v = [1;2;1];
% h = [1,0,-1];
% gy = conv2(conv2(im,h','same'),v','same');
% gx = conv2(conv2(im,v,'same'),h,'same');


%% normalize gradient magnitude
gmag = sqrt(gx.^2+gy.^2);
gmag_norm = (gmag-mean(gmag(:)))./sqrt(var(gmag(:)));
gmag_norm = (gmag_norm-min(gmag_norm(:)))./(max(gmag_norm(:))-min(gmag_norm(:)));

%% initial tangent map
tx = -gy;
ty = gx;

tx(gy<0 & gx>0) = -tx(gy<0 & gx>0);
ty(gy<0 & gx>0) = -ty(gy<0 & gx>0);
tx(gy<0 & gx<0) = -tx(gy<0 & gx<0);
ty(gy<0 & gx<0) = -ty(gy<0 & gx<0);

% ty((gy>0 & gx>0) | (gy<0 & gx<0)) = -ty(gy>0 & gx>0 | (gy<0 & gx<0));
% tx((gy>0 & gx<0) | (gy<0 & gx>0)) = -tx(gy>0 & gx<0 | (gy<0 & gx>0));
% tx(gx==0) = -tx(gx==0);

% normalize
tmag = sqrt(tx.^2+ty.^2);
tx = tx./tmag;
ty = ty./tmag;
tx(tmag==0) = 0;
ty(tmag==0) = 0;