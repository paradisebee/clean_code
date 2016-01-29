im = 10*fspecial('gaussian',10,5);
[r,c] = size(im);
% [circX, circY] = drawSolidCircle([11,11], 5, [r,c]);
% idx = sub2ind(size(im),circY(:),circX(:));
% im(idx) = 1;

figure
imshow(im)

tic
pairwise1 = full(setSmoothness2Dfunc(im,1,0.01,'4nbh'));
toc
tic
pairwise2 = full(setSmoothness2D_split(im,1,00.01));
toc
figure
imshow(pairwise1,[])
figure
imshow(pairwise2,[])