function result = patch_oper(im, im_d, pmask, cmask, win_sz, ovl, func_handle)

addpath /home/bee/Dropbox/Code/IPfunctions

[sub_img, rge] = splitting(im.*cmask, win_sz, ovl);
[sub_edge, rge] = splitting(im_d.*cmask, win_sz, ovl);
[sub_pmask, rge] = splitting(pmask, win_sz, ovl);
[sub_cmask, rge] = splitting(cmask, win_sz, ovl);

%% set parameters for line detection 
% thres = 0.65; sigma = 0; max_L = 15; 
% w = 15;

sub_result = cell(size(sub_img));
for inst = 1:length(sub_img)
    roi = sub_img{inst};
    roi_d = sub_edge{inst};
    roi_pmask = sub_pmask{inst};
    roi_cmask = sub_cmask{inst};
    if ~sum(roi_cmask(:))
        sub_result{inst} = zeros(size(roi)); 
        continue; 
    end
    tmp_result = func_handle(roi,roi_d,roi_pmask,roi_cmask);
    tmp_result = tmp_result;
    tmp_result(:,1:10) = 0;tmp_result(:,end-10:end) = 0;
    tmp_result(1:10,:) = 0;tmp_result(end-10:end,:) = 0;
    sub_result{inst} = tmp_result;
end

result = merging(sub_result, win_sz, ovl, rge, size(im));