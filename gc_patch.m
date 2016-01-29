function gm = gc_patch(roi,roi_d,roi_mask,roi_cmask)

addpath helper_function/GCMex

% roi_f = Norm(roi;
% % roi_f = heav_f(roi, 0.5, 0.8);
% % roi = Norm(roi);
% 
% re = roi_d;
% % [re, ~] = edge_detection(roi_f,roi_mask);
% re = Norm(re);
% re(re<0.5) = 0;
% 
% G = fspecial('gaussian', [3 3], 1);
% reb = imfilter(re,G,'same');
% hf = heav_f(reb);
% 
% % smoothness term
% epsilon = 0.3;
% h_thres = 0.2;
% fh = @(a,b) (heav_f(abs(a-b),epsilon,h_thres)>0.5)*0.9+...
%     (1-heav_f(abs(a-b),epsilon,h_thres)>0.5)*0.1;
% pairwise = setSmoothness2Dfunc(hf,1,0.1,'8nbh',roi_cmask,fh);
% % pairwise = setSmoothness2Dfunc(roi_f, 1, 0.1, '8nbh', roi_cmask);
% 
% % data term
% prob_obj = roi_f;
% prob_bkg = 1-roi_f;
% unary = zeros(2,length(roi_f(:)));
% unary(1,:) = -log(prob_obj(:));
% unary(2,:) = -log(prob_bkg(:));
% unary(unary==0) = 1e-5;
% unary(unary==inf) = 100;
% 
% % initial labels and label cost
% inp = zeros(length(roi_f(:)),1);
% inp = 1-double(roi_f>0.4);
% labelcost = ones(2) - eye(2);

[n_row, n_col] = size(roi);
num_pixels = n_row*n_col;
segclass = zeros(n_row*n_col,1);
labelcost = ones(2) - eye(2);

unary = zeros(2,num_pixels,'single');
roi_f = Norm(roi,roi_cmask);
unary(1,:) = roi_f(:);
unary(2,:) = 1-roi_f(:);

idx = find(roi_mask);
segclass(idx) = 1;
% roi_d = Norm(roi_d,roi_cmask);
% roi_d = (1-roi_d).*roi_cmask;
pairwise = setSmoothness2Dfunc(roi_f,1,0.1,'8nbh');

% graph cuts
[labels E Eafter] = GCMex(segclass, single(unary), pairwise, single(labelcost),0);%segclass
gm = reshape(labels, n_row, n_col);
