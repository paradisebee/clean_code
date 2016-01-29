function [I_edge, I_edge_blur] = edge_detection(im, mask)

[n_row, n_col] = size(im);

I_edge = zeros(n_row, n_col);
thres_levels = 0.1:0.05:0.5;
n_thres = length(thres_levels);
for k = 1:n_thres
    [I_e] = edge(im.*mask, 'canny', thres_levels(k));
    I_edge = I_edge + I_e;
end
I_edge = I_edge/n_thres;

%% Gaussian smoothing on the edge image
% Create Gaussian Filter
G = fspecial('gaussian', [3 3], 1);
% Blur Image
I_edge_blur = imfilter(I_edge,G,'same');