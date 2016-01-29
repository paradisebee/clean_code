function [outIm,whatScale,Direction,Dx,Dy] = try_vesselness_ETF(tx,ty,tmag,I,sigmas,BlackWhite)

tx = tx.*tmag;
ty = ty.*tmag;
gx = ty;
gy = -tx;

% gx(ty<0 & tx>0) = -gx(ty<0 & tx>0);
% gy(ty<0 & tx>0) = -gy(ty<0 & tx>0);
% gx(ty<0 & tx<0) = -gx(ty<0 & tx<0);
% gy(ty<0 & tx<0) = -gy(ty<0 & tx<0);

% Make matrices to store all filterd images
ALLfiltered=zeros([size(I) length(sigmas)]);
ALLangles=zeros([size(I) length(sigmas)]);
ALLIx = ALLangles;
ALLIy = ALLangles;

beta  = 2*0.5^2;
c     = 2*15^2;

%% comparison
[Dxx,Dxy,Dyy] = Hessian2D(I,1);


addpath F:\Dropbox\Code\Download\frangi_filter_version2a
for i = 1:length(sigmas)
    [Dxx,Dxy,Dyy] = Hessian2D(I,sigmas(i));

    [Lambda1,Lambda2,Ix,Iy]=eig2image(Dxx,Dxy,Dyy,gx,gy);
    
    % Compute the direction of the minor eigenvector
    angles = atan2(Ix,Iy);

    % Compute some similarity measures
    Lambda1(Lambda1==0) = eps;
    Rb = (Lambda2./Lambda1).^2;
    S2 = Lambda1.^2 + Lambda2.^2;
   
    % Compute the output image
    Ifiltered = exp(-Rb/beta) .*(ones(size(I))-exp(-S2/c));
    
    if(BlackWhite)
        Ifiltered(Lambda1<0)=0;
    else
        Ifiltered(Lambda1>0)=0;
    end
    
    % store the results in 3D matrices
    ALLfiltered(:,:,i) = Ifiltered;
    ALLangles(:,:,i) = angles;
    ALLIx(:,:,i) = Ix;
    ALLIy(:,:,i) = Iy;
end

% Return for every pixel the value of the scale(sigma) with the maximum 
% output pixel value
if length(sigmas) > 1,
    [outIm,whatScale] = max(ALLfiltered,[],3);
    outIm = reshape(outIm,size(I));
    if(nargout>1)
        whatScale = reshape(whatScale,size(I));
    end
    if(nargout>2)
        Direction = reshape(ALLangles((1:numel(I))'+(whatScale(:)-1)*numel(I)),size(I));
        Dx = reshape(ALLIx((1:numel(I))'+(whatScale(:)-1)*numel(I)),size(I));
        Dy = reshape(ALLIy((1:numel(I))'+(whatScale(:)-1)*numel(I)),size(I));
    end
else
    outIm = reshape(ALLfiltered,size(I));
    if(nargout>1)
        whatScale = ones(size(I));
    end
    if(nargout>2)
        Direction = reshape(ALLangles,size(I));
        Dx = reshape(ALLIx,size(I));
        Dy = reshape(ALLIy,size(I));
    end
end