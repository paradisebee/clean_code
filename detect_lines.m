function [R_combined, R_theta, rmask] = detect_lines(img, pmask, cmask, w, max_L, thres, sigma)
%% detect vessels by detecting lines at different scales
%% set default parameters
if nargin < 4
    thres = 0.65; sigma = 0; max_L = 15; 
    w = 15;
end
%% obtain average image
cim = img.*pmask;
m_filter = ones(w);
I_ave = conv2(cim, m_filter,'same')/sum(m_filter(:));
% I_denorm = max(conv2(double(pmask), m_filter, 'same'), 10);
% I_ave = conv2(cim, m_filter,'same')./I_denorm;

%% obtain image dimension
n_row = size(cim, 1);
n_col = size(cim, 2);

%% Obtian a mask for handling the boundary
erosionsize = 5;
tmask = imerode(cmask, strel('disk', erosionsize, 0));
bmask = pmask - tmask; %% boundary mask
% [px, py] = find(bmask);
% cx = mean(px);
% cy = mean(py);
% theta = angle((py-cy) + sqrt(-1)*(px-cx));
% sidx = find(abs(theta)<pi/4 | abs(theta)>3*pi/4);
% tmask1 = bmask;
% tmask1((py(sidx)-1)*n_row + px(sidx)) = 0;
% tmask2 = bmask - tmask1;
% figure; imagesc(tmask1); colormap gray; axis image;
% figure; imagesc(tmask2); colormap gray; axis image;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Obtain responses at differnt scales L = 3, 5, 7, ..., max_L
%% init R_all
R_all = zeros(12, n_row, n_col);
count = 0;
for L = 3:2:max_L
    [line_temps] = get_line_temps(L, sigma);
    for k = 1:12
        tmask = line_temps{k};
        R = conv2(cim, tmask, 'same')./sum(tmask(:));
        R_all(k,:,:) = R;        
    end
    [Rmax(:,:), ctheta(:,:)] = max(R_all, [], 1);
    count = count+1;
    R = Rmax - I_ave;
    Rim{count} = (R - mean(R(:)))/std(R(:)); % normalization    
    theta{count} = ctheta;
end

%% Combination
R_combined = cim;
theta_combined = -ones(n_row, n_col).*(1-pmask);
for k = 1:count
    R_combined = R_combined + Rim{k}.*(1-bmask);
    theta_combined = theta_combined + theta{k}; 
end
R_combined = R_combined/(count+1);
theta_combined = theta_combined/count;
R_theta = theta_combined;
% figure; subplot(2,2,1); imagesc(cim); axis image; colormap gray; axis off; title('Preprocessed green channel');
% subplot(2,2,2); imagesc(R_combined.*cmask); axis image; colormap gray; axis off; title('Combined response');
% subplot(2,2,3); imagesc(R_combined.*cmask > thres); axis image; colormap gray; axis off; title('Result mask');
% subplot(2,2,4); imagesc(theta_combined.*cmask); axis image; colormap gray; axis off; title('Angle');

rmask = (R_combined.*cmask) > thres;
disp('Line detection is done.');


