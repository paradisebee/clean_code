function [R_combined, R_theta, rmask] = detect_lines_ETF(tx,ty, pmask, cmask, w, max_L, thres, sigma)
%% detect vessels by detecting lines at different scales
%% set default parameters
if nargin < 4
    thres = 0.65; sigma = 0; max_L = 15; 
    w = 15;
end
%% obtain average image
cx = tx.*pmask;
cy = ty.*pmask;
m_filter = ones(w);
I_ave = conv2(sqrt(cx.^2+cy.^2), m_filter,'same')/sum(m_filter(:));
%% obtain image dimension
n_row = size(cx, 1);
n_col = size(cx, 2);

%% Obtian a mask for handling the boundary
erosionsize = 5;
tmask = imerode(cmask, strel('disk', erosionsize, 0));
bmask = pmask - tmask; %% boundary mask

%% Obtain responses at differnt scales L = 3, 5, 7, ..., max_L
%% init R_all
R_all = zeros(12, n_row, n_col);
count = 0;
for L = 3:2:max_L
    [line_temps] = get_line_temps(L, sigma);
    for k = 1:12
        tmask = line_temps{k};
        Rx = conv2(cx, tmask, 'same')./sum(tmask(:));
        Ry = conv2(cy, tmask, 'same')./sum(tmask(:));
        R_all(k,:,:) = sqrt(Rx.^2+Ry.^2);        
    end
    [Rmax(:,:), ctheta(:,:)] = max(R_all, [], 1);
    count = count+1;
    R = Rmax;% - I_ave;
    Rim{count} = (R - mean(R(:)))/std(R(:)); % normalization    
    theta{count} = ctheta;
end

%% Combination
R_combined = zeros(size(cx));
theta_combined = -ones(n_row, n_col).*(1-pmask);
for k = 1:count
    R_combined = R_combined + Rim{k}.*(1-bmask);
    theta_combined = theta_combined + theta{k}; 
end
R_combined = R_combined/(count+1);
theta_combined = theta_combined/count;
R_theta = theta_combined;

rmask = (R_combined.*cmask) > thres;
disp('Line detection is done.');