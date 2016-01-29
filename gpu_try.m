
% a = rand(605,700);
% b = rand(605,700);
% e = zeros(605,700);
I = double(imread('pout.tif'))';
a = I;
b = I;
e = zeros(size(a));
% d = a+b;

a = gpuArray(a);
b = gpuArray(b);
e = gpuArray(e);

cudaFilename = '/home/bee/TRY/add_kernel.cu';
ptxFilename = '/home/bee/TRY/add_kernel.ptx';
kernel = parallel.gpu.CUDAKernel(ptxFilename, cudaFilename);

[c,r] = size(a);

gridDim_x = ceil(c/32);
gridDim_y = ceil(r/32);

kernel.ThreadBlockSize = [32,32];
kernel.GridSize = [gridDim_x,gridDim_y];

[e1,e2] = feval(kernel, e, a, b, r,c);

e1 = gather(e1);
e2 = gather(e2);

figure
e1(e1<0) = 0;
imshow(e1,[])
figure
imshow(e2,[])