function pairwise = setSmoothness2D_split(im,lumbda,sigma)

indexArray1 = [];
indexArray2 = [];
valueArray = [];

[row,col] = size(im);


func_handle = @(a,b) exp(-(a-b).^2/2/sigma^2);

% 4 neighborhood
% i, j+1
cmp_im = nan(size(im));
cmp_im(:,1:end-1) = im(:,2:end);
lgc = ~isnan(cmp_im(:));
valueArray = [valueArray; func_handle(im(lgc), cmp_im(lgc))];
m = find(lgc==1);
indexArray1 = [indexArray1; m];
indexArray2 = [indexArray2; m+row];

% i, j-1
cmp_im = nan(size(im));
cmp_im(:,2:end) = im(:,1:end-1);
lgc = ~isnan(cmp_im(:));
valueArray = [valueArray; func_handle(im(lgc), cmp_im(lgc))];
m = find(lgc==1);
indexArray1 = [indexArray1; m];
indexArray2 = [indexArray2; m-row];

% i-1, j
cmp_im = nan(size(im));
cmp_im(2:end,:) = im(1:end-1,:);
lgc = ~isnan(cmp_im(:));
valueArray = [valueArray; func_handle(im(lgc), cmp_im(lgc))];
m = find(lgc==1);
indexArray1 = [indexArray1; m];
indexArray2 = [indexArray2; m-1];

% i+1, j
cmp_im = nan(size(im));
cmp_im(1:end-1,:) = im(2:end,:);
lgc = ~isnan(cmp_im(:));
valueArray = [valueArray; func_handle(im(lgc), cmp_im(lgc))];
m = find(lgc==1);
indexArray1 = [indexArray1; m];
indexArray2 = [indexArray2; m+1];

valueArray = lumbda * valueArray;
pairwise = sparse(indexArray1, indexArray2, valueArray, row*col, row*col);