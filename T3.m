%
% T3 Demo code
%
I1 = imread('moon.jpg');            % SOURCE IMAGE
I2 = imread('mountain.jpg');        % DESTINATION IMAGE
% I1 = double(repmat(I1,[1,1,3]));
% I2 = double(repmat(I2,[1,1,3]));

PIE_Gui(I1,I2);