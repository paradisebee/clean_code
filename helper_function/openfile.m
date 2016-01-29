function [I, M, G, rgbI] = openfile(db, ds)

% addpath D:/Data/Retinal/DRIVE/test/images
% addpath F:/Retinal/DRIVE/training/images

%% load images DRIVE
if strcmpi('DRIVE', db)
    str1 = 'training'; str2 = 'test';
    str = ds;
    root = ['/home/bee/Research/Retinal/DRIVE/',str,'/'];
    % root = 'F:/Retinal/DRIVE/test/';
    I_all = dir([root,'images/*_',str,'.tif']);
    M_all = dir([root,'mask/*_',str,'_mask.gif']);
    if strcmpi(str,str2)
        G_all = dir([root,'1st_manual/*_manual1.gif']);
    else
        G_all = dir([root,'manual1/*_manual1.gif']);
    end
    for i = 1:length(I_all)
        im_rgb = imread([root, 'images/', I_all(i).name]);
        im = double(im_rgb);
        im = im./max(im(:));
        im = im(:,:,2);
        if ~exist('I','var')
            I = zeros(size(im,1),size(im,2),length(I_all));
            M = I;
            G = I;
            rgbI = cell(1,length(I_all));
        end
        I(:,:,i) = im;
        rgbI{i} = im_rgb;

        m = double(imread([root, 'mask/', M_all(i).name]));
        m = m(:,:,1);
        m(m~=0) = 1;
        M(:,:,i) = m;
        if strcmpi(str,str2)
            g = double(imread([root, '1st_manual/', G_all(i).name]));
        else
            g = double(imread([root, 'manual1/', G_all(i).name]));
        end
        g = g(:,:,1);
        g(g~=0) = 1;
        G(:,:,i) = g;
    end
end

%% load images STARE
if strcmpi('STARE', db)
    str1 = 'labels-vk'; str2 = 'labels-ah';
    str = 'stare-images';
    root = ['/home/bee/Research/Retinal/STARE/'];
    % root = 'F:/Retinal/DRIVE/test/';
    I_all = dir([root,str,'/*.ppm']);
%     M_all = dir([root,'mask/*_',str,'_mask.gif']);
    if strcmpi(str1,ds)
        G_all = dir([root,str1,'/*.ppm']);
    else
        G_all = dir([root,str2,'/*.ppm']);
    end
    for i = 1:length(I_all)
        im_rgb = imread([root, str, '/', I_all(i).name]);
        im = double(im_rgb);
        m = double(sum(im,3)>100);
        im = im./max(im(:));
        im = im(:,:,2);
        if ~exist('I','var')
            I = zeros(size(im,1),size(im,2),length(I_all));
            M = I;
            G = I;
            rgbI = cell(1,length(I_all));
        end
        I(:,:,i) = im;
        rgbI{i} = im_rgb;
        
%         m = double(imread([root, 'mask/', M_all(i).name]));
%         m = m(:,:,1);
%         m(m~=0) = 1;
        M(:,:,i) = m;
        if strcmpi(str1,ds)
            g = double(imread([root, str1, '/', G_all(i).name]));
        else
            g = double(imread([root, str2, '/', G_all(i).name]));
        end
        g = g(:,:,1);
        g(g~=0) = 1;
        G(:,:,i) = g;
    end
end