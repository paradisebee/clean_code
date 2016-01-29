function [im_out, mu, sigma] = illumination_correction(im_in, im_mask, N)
% illumination correction

if nargin < 3
    N =250;
end

%% boundary padding
dim = size(im_in);

[rm, cm] = find(im_mask);
idm = find(rm==cm, 1);
sub_im = im_in(rm(idm):rm(idm)+300,cm(idm):cm(idm)+300);

im_temp = repmat(sub_im, ceil((dim(1)+N*2)/300), ceil((dim(2)+N*2)/300));
im_temp = im_temp(1:dim(1)+N*2,1:dim(2)+N*2);

im_large = zeros(dim+N*2);
im_large(N+1:end-N,N+1:end-N) = im_mask;
lgc = (im_large==0);
im_large(N+1:end-N,N+1:end-N) = im_in;
im_large(lgc) = im_temp(lgc);


%% simple correction
im_bp = im_large;

h = 1/N*ones(N);

mu = conv2FFT(h, im_bp);
mu = (mu-min(mu(:)))./(max(mu(:))-min(mu(:)));
h = 1/26*ones(26);
sigma = sqrt(conv2FFT(h, (im_bp-mu).^2));

ni = (im_bp-mu)./sigma;
ni = ni(N+1:end-N,N+1:end-N);%.*im_mask;
% figure
% imshow(ni,[]);

im_out = (ni-min(ni(:)))./(max(ni(:))-min(ni(:))).*im_mask;
sigma = (sigma-min(sigma(:)))./(max(sigma(:))-min(sigma(:)));
% figure
% imshow(mu,[]);
% figure
% imshow(sigma,[]);

%%


im_out = 1./(1+exp(-20.*(im_out-mean(im_out(im_mask==1))))).*im_mask;
mu = mu(N+1:end-N,N+1:end-N);
sigma = sigma(N+1:end-N,N+1:end-N);
% figure
% imshow(im_h.*im_mask,[])