%% main script
addpath 'helper_function'

[I, M, G, rgbI] = openfile('DRIVE', 'test');
% [I, M, G, rgbI] = openfile('STARE', 'labels-vk');
% for i = 1:size(I,3)
%     im = I(:,:,i);
%     imshow(im,[])
%     pause;
% end
%% set parameters for line detection 
thres = 0.65; sigma = 0; max_L = 15; 
w = 15;
n = 1;

tic
percentage = []; threshold = [];
result0 = zeros(size(I));
result1 = result0;
for i = 1:size(I,3)
    % get current instance
    im = I(:,:,i); cmask = M(:,:,i); vmask = G(:,:,i);
    %% preprocessing
%     [im, mask_out] = boundary_padding(im, cmask, 'draw');
    [img, img_raw, pmask, bigimg] = vessel_preproc(im, cmask);
    img = imtophat(1-img,strel('disk',12));
    bigimg = imtophat(1-bigimg,strel('disk',12));
    bigmask = zeros(size(bigimg));
    bigmask(51:end-50, 51:end-50) = pmask;
%     figure
%     imshow(im_out,[])
%     figure
%     imshow(mask_out)
    sz = 50; ovl = 10;
    [Roi, rge] = splitting(img, sz, ovl);
    imgmask = merging(Roi, sz, ovl, rge, size(img));
    r = rge{1}; c = rge{2};
    Mask = {};
    tic
    for j = 1:length(Roi)
        roi = Roi{j};
        [ri,ci] = ind2sub([length(r),length(c)],j);
        %% for debug only
%         if ri~=5 || ci~=4
%             continue;
%         end
        r_range = r(ri)+50-10:r(ri)+50+sz-1+10;
        c_range = c(ci)+50-10:c(ci)+50+sz-1+10;
        bigroi = bigimg(r_range,c_range);
        bigroimask = bigmask(r_range,c_range);
        if sum(bigroimask(:)==0)>0.95*length(bigroimask(:))
            continue;
        end
        bigroi = Norm(bigroi,bigroimask).*bigroimask;
        % procedure
        thres = 0.65; sigma = 0; max_L = 15; 
        w = 15;
        [R_combined, R_theta, rmask] = detect_lines(bigroi, bigroimask, bigroimask, w, max_L, thres, sigma);
        R_combined(isnan(R_combined)) = 0;
        
        [tx, ty, gmag_norm] = get_tangent(R_combined);
        [tx_new, ty_new] = ETF(tx,ty,gmag_norm,bigroi);
        [tx_new, ty_new] = ETF(tx_new,ty_new,gmag_norm,bigroi);
        [tx_new,ty_new,score1] = ETFstraight(false, 0,tx_new,ty_new,bigroi);
        % *************correct opposite direction**************
        lgc = acos([tx_new(:),ty_new(:)]*[0;-1]) < 0.35;
        ty_new(lgc) = -ty_new(lgc);
        [tx_new,ty_new,score1] = ETFstraight(false, 1,tx_new,ty_new,bigroi);
        
        [X,Y]=meshgrid(1:size(bigroi,2),1:size(bigroi,1));

        %%
        lgc = bigroimask(:)==1;
        feat = [bigroi(lgc),R_combined(lgc),tx_new(lgc),ty_new(lgc)];%,X(:),Y(:)];
        feat(:,1) = Norm(feat(:,1));
        feat(:,2) = Norm(feat(:,2));
        feat(:,3) = Norm(feat(:,3));
        feat(:,4) = Norm(feat(:,4));
        
        
        [coeff,score] = princomp(feat);
        roimask = zeros(size(bigroi));
        lgc(lgc==1) = score(:,2)>0.2;
        roimask(lgc==1) = 1;
        thresroi = bigroi>0.5;
        if sum(sum(roimask & thresroi))/sum(sum(thresroi)) < 0.3
            lgc = bigroimask(:)==1;
            roimask = zeros(size(bigroi));
            lgc(lgc==1) = score(:,1)>0.2;
            roimask(lgc==1) = 1;
        end
        Mask{j} = roimask(11:end-10,11:end-10);
        
%         figure
%         imshow(R_combined.*bigroimask,[])
%         figure
%         imshow(bigroi.*bigroimask,[])
%         figure
%         imshow(bigroi.*bigroimask,[])
%         hold on
%         quiver(tx_new,ty_new)
%         figure
%         plot(score(:,1),score(:,2),'+');
%         figure
%         imshow(roimask)
%         
        close all
    end
    toc
    imgmask = merging(Mask, sz, ovl, rge, size(img));
    figure
    imshow(imgmask)
end