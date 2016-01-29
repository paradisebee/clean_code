addpath F:\Retinal\matlabfiles
load etf_ls1.mat
lm = lak.*mask;
BW = zeros(size(im));
cc = bwconncomp(lak.*mask);
for i = 1:cc.NumObjects
    idx = cc.PixelIdxList{i};
    if length(idx)>30
        if mean(im(idx))<0.5
            BW(idx) = 1;
        end
    end
end

figure
imshow(BW.*cmask)

old_mask = zeros(size(img));
new_mask = BW.*cmask;
n = 0;
tic
%% iterative process
while sum(new_mask(:)-old_mask(:))>0.01*sum(old_mask(:))
    old_mask = new_mask;
cc = bwconncomp(old_mask);
se = strel('disk',1);
Region = cell(cc.NumObjects,1);
Seeds = cell(cc.NumObjects,1);
for i = 1:cc.NumObjects
    idx = cc.PixelIdxList{i};
    tmp = zeros(size(im));
    tmp(idx) = 1;
    tmp = imdilate(imdilate(imdilate(imdilate(tmp,se),se),se),se);
    tmp = tmp.*cmask;
    Region{i} = find(tmp==1);
    Seeds{i} = idx;
end

%%
proc_im = 1-img;
final_mask = zeros(size(img));
for i = 1:length(Region)
    pts = Region{i};
    bigmask = zeros(size(proc_im));
    bigmask(pts) = 1;
    seedP = Seeds{i};
    [pr,pc] = ind2sub(size(im),pts);
    roi = proc_im(min(pr):max(pr),min(pc):max(pc));
    roimask = bigmask(min(pr):max(pr),min(pc):max(pc));
    [sr,sc] = ind2sub(size(im),seedP);
    r_sidx = sub2ind(size(roi), sr-min(pr)+1,sc-min(pc)+1);
    % set data terms
    prob_obj = roi(:);
    prob_obj(roimask(:)==0) = 0;
    prob_obj(r_sidx) = 1;
    prob_bkg = 1-roi(:);
    prob_bkg(roimask(:)==0) = 1;
    prob_bkg(r_sidx) = 0;
    unary = zeros(2,length(roi(:)));
    unary(2,:) = -log(prob_obj(:));
    unary(1,:) = -log(prob_bkg(:));
    unary(unary==0) = 1e-5;
%     if n == 0
        unary(unary==inf) = 100;
%     else
%         unary(unary==inf) = 2;
%     end
    % set smoothness term
    lumbda = 1; sigma = 0.05; choice = '8nbh';
    pairwise = setSmoothness2Dfunc(roi,lumbda,sigma,choice,roimask);
    % initial labels and label cost
    inp = zeros(1,length(roi(:)));
    inp(r_sidx) = 1;
    labelcost = ones(2) - eye(2);

    % graph cuts
    [labels E Eafter] = GCMex(inp(:), single(unary), pairwise, single(labelcost),0);%segclass
%     gm = zeros(size(roi));
%     gm(pts) = labels;
    gm = reshape(labels, size(roi,1), size(roi,2));
%     figure
%     imshow(gm)
    final_mask(min(pr):max(pr),min(pc):max(pc)) = ...
        final_mask(min(pr):max(pr),min(pc):max(pc)) | gm;
end
    new_mask = final_mask;
    figure
    imshow(new_mask)
    n = n+1
    toc
%%
end
figure
imshow(final_mask)
