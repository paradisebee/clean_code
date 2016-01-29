%% main script
addpath 'helper_function'
addpath 'ctrlin'

[I, M, G, rgbI] = openfile('DRIVE', 'test');
% [I, M, G, rgbI] = openfile('STARE', 'labels-vk');

thres = 0.65; sigma = 0; max_L = 15; 
w = 15;
for i = 1:size(I,3)
    % get current instance
    im = I(:,:,i); cmask = M(:,:,i); vmask = G(:,:,i);
%     %% preprocessing
    [img, img_raw, pmask] = vessel_preproc(im, cmask);
    img = imtophat(1-img,strel('disk',12));
    [tx, ty, gmag_norm] = get_tangent(img);
    [tx_new, ty_new] = ETF(tx,ty,gmag_norm,img);
    [tx_new, ty_new] = ETF(tx_new,ty_new,gmag_norm,img);
%     % reverse tx ty
%     [rtx_new, rty_new] = ETF(-tx,-ty,gmag_norm,img);
%     [rtx_new, rty_new] = ETF(rtx_new,rty_new,gmag_norm,img);
    [R_combined, R_theta, rmask] = detect_lines(img, pmask, cmask, w, max_L, thres, sigma);
    
    [bif, rads, vpts] = bifurcation_detection(R_combined,cmask,tx_new,ty_new,15);
    
    
    %%
%     [im, mask_bdy] = preprocess(im, cmask);
%     im = Norm(im,cmask);
%     [im_out, mu, sigma] = illumination_correction(im, mask_bdy, 100);
% 
%     roi = im_out;
%     img = Norm(roi,cmask);

    startP = [181,303;
              222,346;
              191,306;
              403,466];
    endP = [148,259;
            219,377;
            209,289;
            370,476];
      %%  
      figure
      imshow(img.*cmask,[])
      hold on
      plot(startP(1:3,2),startP(1:3,1),'ro','MarkerSize',10)
      plot(endP(1:3,2),endP(1:3,1),'yo','MarkerSize',10)
      for sli = 1:3
        if sli == 2
             flag = -1;
        else
            flag = 1;
        end
        [s, radius] = ctrlin_trace(img, tx_new, ty_new, startP(sli,:), endP(sli,:), 3, 10, flag);
        plot(s(:,2),s(:,1),'LineWidth',3,'MarkerSize',5);
        for ff = 1:size(s,1)
            ss = round(s(ff,:));
            [x, y, pos] = line2D(ss, [tx_new(ss(1),ss(2)), ty_new(ss(1),ss(2))],...
                radius(ff));
%             plot(x,y,'g-','lineWidth',3);
        end
      end
    
    %%
    
    exu = exudates_detection(im, cmask);
    
    figure
    subplot(1,2,1);
    imshow(im,[])
    subplot(1,2,2);
    imshow(exu)
    
end