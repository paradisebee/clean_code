
% close all; clear;
%% Assume the input file is a raw fundus image 
% [filename, data_dir] = uigetfile('.tif', 'Pick a fundus image'); 
% if filename~=0
%     im_filename = strcat(data_dir, filename);
%     im = imread(im_filename);
% else
%     disp('No image file has been selected!');
%     return;
% end
% disp('The selected fundus image has been loaded.');
% figure; imagesc(im); axis image; axis off; title('Input fundus image');
% 
% %% load the mask and groundtruth vessel mask (follow the naming convention of the DRIVE database)
% cmask_filename = strcat(data_dir(1:end-7), 'mask\', filename(1:end-4),'_mask.gif');
% cmask = imread(cmask_filename);
% cmask = cmask > 0;
% %% differentiate between test and training
% if isempty(strfind(data_dir, 'test'))
%     vmask_filename = strcat(data_dir(1:end-7), 'manual1\', filename(1:2),'_manual1.gif');
% else
%     vmask_filename = strcat(data_dir(1:end-7), '1st_manual\', filename(1:2),'_manual1.gif');
% end
% vmask = imread(vmask_filename);
% disp('All files have been imported.');
% 
% %% Pre-processing to obtain an enhance green channel image
% [img, img_raw, pmask] = vessel_preproc(double(im)./255, cmask);
% img = 1 - img; figure; imagesc(img.*pmask); colormap gray; axis image; axis off; title('Pre-processed green channel image');


%%
function result = etf_plus_fdog(img,pmask,cmask)

% addpath F:\Dropbox\Code\Retinal\new
% addpath F:\Dropbox\Code\IPfunctions
% addpath F:\Dropbox\supervisor\New\code\to_Bichao\m_code

[tx, ty, gmag_norm] = get_tangent(img);
[tx_new, ty_new] = ETF(tx,ty,gmag_norm,img);
[tx_new, ty_new] = ETF(tx_new,ty_new,gmag_norm,img);
[tx_new, ty_new] = ETF(tx_new,ty_new,gmag_norm,img);
%%
% tic
sigma_C = 1;sigma_M = 3;
Hmap = fDoG(sigma_C, sigma_M, img, tx_new, ty_new, pmask);
Hmap = Norm(Hmap.*pmask);
% figure
% imshow(Hmap.*img,[])
% toc
%% Line detection on pre-processed image
% img = imtophat(img,strel('disk',12));
w = 15; thres = 0.65; sigma = 0; max_L = 15; 
[R_combined1, R_theta1, rmask1] = detect_lines(img, pmask, cmask, w, max_L, thres, sigma);
[R_combined, R_theta, rmask] = detect_lines(Hmap, pmask, cmask, w, max_L, thres, sigma);
% [Acc, Sp, Se, MCC] = evaluation(vmask, rmask1, cmask)
%%
% cc = bwconncomp(rmask1);
% for i = 1:cc.NumObjects
%     idx = cc.PixelIdxList{i};
%     if length(idx)<=30
%         rmask1(idx) = 0;
%     end
% end
% cc = bwconncomp(rmask);
% for i = 1:cc.NumObjects
%     idx = cc.PixelIdxList{i};
%     if length(idx)<=30
%         rmask(idx) = 0;
%     end
% end
result = rmask | rmask1;
cc = bwconncomp(result);
for i = 1:cc.NumObjects
    idx = cc.PixelIdxList{i};
    if length(idx)<=30
        result(idx) = 0;
    end
end
% [Acc, Sp, Se, MCC] = evaluation(vmask, rmask1, cmask)
% [Acc, Sp, Se, MCC] = evaluation(vmask, rmask, cmask)
% [Acc, Sp, Se, MCC] = evaluation(vmask, result, cmask)

%%
% se = strel('disk',1);
% rr = imerode(result,se);
% [Acc, Sp, Se, MCC] = evaluation(vmask, rr, cmask)