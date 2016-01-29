function gmask = graph_cut(d_term,s_term,rmask,pmask,vmask)

d_term = Norm(d_term);
%% Pre-processing to obtain an enhance green channel image
% [img, img_raw, pmask] = vessel_preproc(double(im)./255, pmask);
% img = 1 - img; figure; imagesc(img.*pmask); colormap gray; axis image; axis off; title('Pre-processed green channel image');
% 
% %% Line detection on pre-processed image
% img = imtophat(img,strel('disk',12));
% w = 15; thres = 0.65; sigma = 0; max_L = 15; 
% [R_combined, R_theta, rmask] = detect_lines(img, pmask, pmask, w, max_L, thres, sigma);

%% Evaluation of the global thresholding results
% [Acc, Sp, Se, MCC] = evaluation(vmask, rmask, pmask)

%% Implement graphcut segmentation on the line response ONLY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('helper_function/GCMex');
%% set up the graph
[n_row, n_col] = size(d_term);
segclass = zeros(n_row*n_col,1);
labelcost = ones(2) - eye(2); %% Label Cost
%% Use line response threaholding to initialize the labels
num_pixels = n_row*n_col;
vec_img = reshape(d_term.*pmask, num_pixels, 1);
vec_mask = reshape(rmask.*pmask, num_pixels, 1);
[idx] = find(vec_mask);
segclass(idx) = 1;

%% any graph
unique_pixels = 1+n_row:(n_row-1)+(n_col-2)*n_row;
row_pos = kron(ones(1,4), unique_pixels);
col_pos = [unique_pixels+1, unique_pixels+n_row, unique_pixels+n_row+1, unique_pixels-n_row+1];
%% Set smoothness cost 
% pairwise = 5*sparse(row_pos,col_pos,ones(size(row_pos)),num_pixels,num_pixels);

%% Edge detection for setting the smoothness cost 


%% Set data cost
unary = zeros(2, num_pixels,'single');
%% use heaviside function 
% epsilon = 0.5;
% h_thres = 0.8;
% f = 0.5*(1+2/pi*atan((vec_img-h_thres)/epsilon));
f = vec_img;
% figure; imagesc(reshape(f, n_row, n_col)); axis equal off tight; colormap gray; title('nonlinear transform');
unary(1,:) = f;
unary(2,:) = 1-f;

%% smoothness term
% af=reshape(f, n_row, n_col);
% [Ie, Ieb] = edge_detection(af, pmask);
% fh = @(a,b) 1-a; % if a==b, return 1, else return 0.2
fh = @(a,b) (a==b)*0.1+(a~=b)*1;
% pairwise = setSmoothness2Dfunc(d_term,1,0.5,'8nbh');
% pairwise = setSmoothness2Dfunc(s_term,1,0.5,'8nbh',fh);
pairwise = setSmoothness2Dfunc(d_term,1,0.1,'8nbh',pmask);
% pairwise = setSmoothness2Dfunc(R_theta,1,0.2,'8nbh');
% lumbda = setSmoothness2Dfunc(s_term,1,0.5,'8nbh',fh);
% pairwise = pairwise.*lumbda;%.*thet;

%% set hard constraints
% idx = find(vec_img < 0.1);
% unary(1,idx) = 0;
% unary(2,idx) = 100;
% idx = find(vec_img > 2*h_thres);
% unary(1,idx) = 100;
% unary(2,idx) = 0;
% inp = double(d_term.*pmask > 1.5);
[labels E Eafter] = GCMex(segclass, single(unary), pairwise, single(labelcost),0);%segclass
gmask = reshape(labels, n_row, n_col);
% figure; imagesc(d_term.*pmask > 1.5); axis equal off tight; colormap gray; title('input');
figure; imshow(gmask); %axis equal off tight; colormap gray; title('output');
%% Evaluation of the results
[Acc, Sp, Se, MCC] = evaluation(vmask, gmask, pmask)
[Acc, Sp, Se, MCC] = evaluation(vmask, rmask, pmask)





%% end of graphcut segmentation on the line response ONLY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Implement graphcut segmentation on the line response and the enhanced green channel image %%%%%
% u1 = reshape(unary(1,:), size(img,1),size(img,2));
% u2 = reshape(unary(2,:), size(img,1),size(img,2));
% figure
% imshow(u1,[])
% figure
% imshow(u2,[])


%%
% win_sz = 50;
% % [sub_img, rge] = splitting(R_combined.*pmask, win_sz, 0);
% [sub_img, rge] = splitting(d_term.*pmask, win_sz, 0);
% [sub_edge, rge] = splitting(s_term, win_sz, 0);
% [sub_theta, rge] = splitting(R_theta, win_sz, 0);
% se = strel('disk',1);
% gt_edge = double(imdilate(vmask,se)-vmask);
% gt_edge = imfilter(gt_edge,G,'same');
% [sub_gt, rge] = splitting(gt_edge, win_sz, 0);
% % lbls = R_combined.*pmask > 1.5;
% lbls = result;
% [sub_lbl, rge] = splitting(lbls, win_sz, 0);
% sub_gmask = cell(size(sub_img));
% % %%
% % roi = a(311:361,289:314);
% % roi_edge = I_edge_blur(311:361,289:314);
% % roi_theta = R_theta(311:361,289:314);
% % roi_lbl = lbls(311:361,289:314);
% for i = 1:length(sub_img)
%     roi = sub_img{i};
%     f = Norm(0.5*(1+2/pi.*atan((roi-h_thres)./epsilon)));
%     unary = [];
%     unary(1,:) = 10*f(:);
%     unary(2,:) = 10*(1-f(:));
% %     pairwise = setSmoothness2Dfunc(f,1,0.01,'8nbh');
% %     pairwise = setSmoothness2Dfunc(roi_theta,1,0.2,'8nbh');
% %     pairwise = setSmoothness2Dfunc(Norm(roi_edge),1,0.5,'8nbh',fh);
%     pairwise = setSmoothness2Dfunc(sub_theta{i},1,0.2,'8nbh');
%     lumbda = setSmoothness2Dfunc(Norm(sub_edge{i}),1,0.5,'8nbh',fh);
%     pairwise = pairwise.*sqrt(lumbda);
% %     pairwise = setSmoothness2Dfunc(Norm(sub_edge{i}),1,0.5,'8nbh',fh);
%     idx = find(roi(:) < 0.1);
%     unary(1,idx) = 0;
%     unary(2,idx) = 100;
%     idx = find(roi(:) > 2*h_thres);
%     unary(1,idx) = 100;
%     unary(2,idx) = 0;
%     inp = double(sub_lbl{i});
% %     inp = double(roi_lbl);
%     [labels E Eafter] = GCMex(inp(:), single(unary), pairwise, single(labelcost),0);%segclass
%     sub_gmask{i} = reshape(labels, size(roi,1), size(roi,2));
% %     figure
% %     imshow(reshape(labels, size(roi,1), size(roi,2)));
% end
% final_Re = merging(sub_gmask, win_sz, 0, rge, size(img));
% figure
% imshow(final_Re);
% [Acc, Sp, Se, MCC] = evaluation(vmask, final_Re, pmask)