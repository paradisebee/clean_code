function exu = exudates_detection(im, mask)

im = Norm(im.*mask);
h = 1/9*ones(3*3);
    
mu = conv2FFT(h, im);
mu = (mu-min(mu(:)))./(max(mu(:))-min(mu(:)));
h = 1/9*ones(3*3);
sigma = sqrt(conv2FFT(h, (im-mu).^2));
thres = 0.8;
% figure
% imshow(mu,[])
% figure
% imshow(sigma,[])

bw = ((mu > 0.7) & (sigma > 0.05)) .* mask;
if sum(bw(:)) > 0.05*sum(mask(:))
    exu = zeros(size(im));
    return;
end
cc = bwconncomp(bw);

% figure
% imshow(bw)

addpath F:\Dropbox\Code\Retinal\neat_code

% group adjacent pixels
idx = [];
index = [];
numPixels = cellfun(@numel,cc.PixelIdxList);
for i = 1:cc.NumObjects
    idx = [idx; cc.PixelIdxList{i}];
    index = [index; i*ones(numPixels(i),1)];
end
[r,c] = ind2sub(size(im),idx);
uind = unique(index);
for i = 1:length(uind)
    pr = r(index==uind(i),:);
    pc = c(index==uind(i),:);
    cmin = inf;
    gid = 0;
    for j = 1:length(uind)
        if j == i
            continue;
        end
        cr = r(index==uind(j),:);
        cc = c(index==uind(j),:);
        [cmin, pos] = min([cmin, min(min(dist([cr,cc],[pr,pc]')))]);
        if pos == 2
            gid = uind(j);
        end
    end
    if cmin < 10
        index(index==uind(i)) = gid;
    end
end

exu = zeros(size(im));
uind = unique(index);
for i = 1:length(uind)
    if sum(index==uind(i)) < 5
        continue;
    end
    tmp_mask = zeros(size(im));
    tmp_mask(idx(index==uind(i))) = 1;
    region_mask = getConvHull(tmp_mask);
    exu = exu | region_mask;    
end