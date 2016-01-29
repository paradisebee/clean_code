%%
% Author: Michael S. Brown, National University of Singapore
%
% Simple Matlab GUI for selection a Region of Interest from Image1, and
% pasting this into Image2.
%
% It calls a function called: PossionImageIntegration() 
%       Currently this function is just a dummy function adds the ROI and
%       Dest together.  This should generate the PIE image.
%
% Main functions: PIE_Gui
%   - Call back functions: myButtonPressDown, myButtonPressUp,
%   myMouseMotion, myKeypress
%
function PIE_Gui(I1, I2)

%
% -------------------  SELECT SOURCE ROI ------------------
%
disp('USAGE: select a polygon region by using left mouse clicks to draw the vertices');
disp('       right nouse click to finish');

h = figure('MenuBar', 'none', 'Toolbar', 'none');  % open window
[BW, xi, yi] = roipoly(I1);                         % this returns a binary image with white (1) in the mask

% extract mask (crop image)
[r,c] = find(BW == 1);                      % find the max values
maxH = max(r) - min(r);                     % extract the height
maxW = max(c) - min(c);                     % extract the width
Ic = imcrop(I1,[min(c) min(r) maxW maxH]);  % crop the image in the RIO

% crop mask - make the mask RGB (3 layers)
Mc = zeros(size(Ic));                       % make a copy of Ic 
Mc(:,:,1) = imcrop(BW,[min(c) min(r) maxW maxH]);
Mc(:,:,2) = imcrop(BW,[min(c) min(r) maxW maxH]);
Mc(:,:,3) = imcrop(BW,[min(c) min(r) maxW maxH]);

% get laplacian of the source image
H =fspecial('laplacian',0);
LAP = imfilter(double(I1),H,'same');
Id = imcrop(LAP,[min(c) min(r) maxW maxH]);

% multiple the Mask by the Image to get only the pixels in the RIO
LAP = immultiply(Id,Mc);
lap = immultiply(double(Ic),Mc);

%
% NOW SELECT PLACE TO PASTE
%
imshow(I2(:,:,1),[]);
title('Click and drag region to desired location. Press any key to integrate, press q to quit');
lh = line(xi, yi, 'Marker','.','LineStyle','-', 'Color', 'r', 'LineWidth',2);

% Set up units and callback functions
set(h, 'Units', 'pixels');
set(h,'WindowButtonDownFcn',@myButtonPressDown);
set(h,'WindowButtonUpFcn',@myButtonPressUp);
set(h, 'WindowButtonMotionFcn', @myMouseMotion);
set(h, 'KeyPressFcn', @myKeyPress);

myData.xi = xi-min(xi);
myData.yi = yi-min(yi);
myData.lap = lap;
myData.LAP = LAP;
myData.DEST = I2;
myData.MASK = Mc;
myData.pressDown = 0;
myData.line = lh;
myData.curX = -1;
myData.curY = -1;

set(h, 'UserData', myData); 

return


%% 
% When button is pressed, call this function
%
function myButtonPressDown(obj,event_obj)

    myData = get(obj, 'UserData');      % get the user data (variable name does not have to be "myData"
    myData.pressDown = 1;               % set mouse press = true
    p = get(gca,'CurrentPoint');        % get current position of mouse on the image
    curX = p(1,1);                      % extract the X position (it's a floating point value)
    curY = p(1,2);                      % extract the Y positions 
    myData.curX = curX;
    myData.curY = curY;
    set(myData.line,'XData', myData.xi+curX, 'YData', myData.yi+curY);

    % Save the myData variable back to the object
    set(obj, 'UserData', myData);
return

%% 
% When button is released, call this function
%
function myButtonPressUp(obj,event_obj)

    myData = get(obj, 'UserData');  % get the user data
    myData.pressDown = 0;           % set mouse press to be false
    set(obj, 'UserData', myData);   % set the uer data (i.e. record mouse is not longer being pressed)
    
return

%% 
% Called anytime the mouse is moved
% 
function myMouseMotion(obj,event_obj)

    myData = get(obj, 'UserData');  % get the user data
  
    if (myData.pressDown == 1)              % we are only interested if the mouse is down
        p = get(gca,'CurrentPoint');        % get the current point from the image
        curX = p(1,1);                      % extract the point from the strange matlab datastructure return by previous line of code
        curY = p(1,2);
        set(myData.line,'XData', myData.xi+curX, 'YData', myData.yi+curY);                                       
        myData.curX = curX;
        myData.curY = curY;
        set(obj, 'UserData', myData);
    end
return


%% 
% Call when key any pressed any key 
%
function myKeyPress(obj, event_obj)

    if (event_obj.Key == 'q')
        close(obj);
        return;
    end
    
     % Update the userdata in the object
    myData = get(obj, 'UserData');
    if (myData.pressDown == 0)          % if mouse is not pressed

        if (myData.curX == -1)
            disp('Select a location');
            return;
        end
 
        %
        % Get the source and destination image
        % Compute a new image (SImage) where the source is translated to
        % the correct position based on the last mouse position.
        % 
        %
        DEST = myData.DEST;
        lap = myData.lap;
        LAP = myData.LAP;
        Mc = myData.MASK;
        tx = round(myData.curX);
        ty = round(myData.curY);
        
        [hh ww depth] = size(lap);
        SImage = zeros(size(DEST));
        OImage = zeros(size(DEST));
        Mask = zeros(size(DEST));
        Target = DEST;
        
        
        SImage( ty:(ty+hh-1), tx:(tx+ww-1), 1 ) =  LAP(:,:,1);
        SImage( ty:(ty+hh-1), tx:(tx+ww-1), 2 ) =  LAP(:,:,2);
        SImage( ty:(ty+hh-1), tx:(tx+ww-1), 3 ) =  LAP(:,:,3);
        OImage( ty:(ty+hh-1), tx:(tx+ww-1), 1 ) =  lap(:,:,1);
        OImage( ty:(ty+hh-1), tx:(tx+ww-1), 2 ) =  lap(:,:,2);
        OImage( ty:(ty+hh-1), tx:(tx+ww-1), 3 ) =  lap(:,:,3);
        Mask( ty:(ty+hh-1), tx:(tx+ww-1), 1 ) =  Mc(:,:,1);
        Mask( ty:(ty+hh-1), tx:(tx+ww-1), 1 ) =  Mc(:,:,2);
        Mask( ty:(ty+hh-1), tx:(tx+ww-1), 1 ) =  Mc(:,:,3);     
        
        idx = find(Mask(:,:,1)==1);
        for i = 1:3
            tmp = Target(:,:,i);
            lap = OImage(:,:,i);
            tmp(idx) = lap(idx);
            Target(:,:,i) = tmp;
        end
        imwrite(uint8(Target), 'source.tif');
        
        % Call the PIE function.  It will returned the integrated image
        newI = PossionImageIntegration(SImage, DEST, Mask, tx, ty, ww, hh);
        
        figure;
        image(newI);
        
    end
      
return

%
%    YOUR FUNCTION HERE ---
%
function I = PossionImageIntegration(LAP, DEST, MASK, tx, ty, ww, hh)
    
    [H,W,C] = size(LAP);
    LAP = double(LAP);
    DEST = double(DEST);
    MASK = double(MASK);
        
        % find DEST boundary points
    mask = MASK(:,:,1);
    se = strel('disk',1);
    outer_bp = imdilate(mask,se) - mask;

    fd_bp = DEST.*repmat(outer_bp,[1,1,3]);   % f* boundary points (DEST)
    
    
    % only process roi with mask
    roi_mask = mask(ty-1:ty+hh,tx-1:tx+ww);     % bw roi mask
    I = DEST;
    dest_roi = I(ty-1:ty+hh,tx-1:tx+ww,1:3);    % destination roi
    f_idx = find(roi_mask==1);
    [Row,Col] = ind2sub(size(roi_mask),f_idx);
    
    roi_bp = fd_bp(ty-1:ty+hh,tx-1:tx+ww,1:3);  % destination boundary intensity
    roi_outer_bp = outer_bp(ty-1:ty+hh,tx-1:tx+ww); % roi outer boundary
    roi_lap = LAP(ty-1:ty+hh,tx-1:tx+ww,1:3);   % source laplacian roi

    f = zeros(size(roi_bp));
    % deal with RGB channels separately
    for chn = 1:3
        dest_bp = roi_bp(:,:,chn);
        lap = roi_lap(:,:,chn);
        indexArray1 = [];
        indexArray2 = [];
        valueArray = [];
        b = zeros(length(f_idx),1);
        % create linear system using sparse matrix
        for i = 1:length(f_idx)
            n = 0;
            r = Row(i);
            c = Col(i);
            b(i) = b(i)-lap(r,c);
            if roi_mask(r+1,c)
                n = n+1;
                indexArray1(end+1) = i;
                indexArray2(end+1) = find(f_idx==(f_idx(i)+1));
                valueArray(end+1) = -1;
            else
                if roi_outer_bp(r+1,c)
                    n = n+1;
                    b(i) = b(i)+dest_bp(r+1,c);
                end
            end
            if roi_mask(r-1,c),
                n = n+1;
                indexArray1(end+1) = i;
                indexArray2(end+1) = find(f_idx==(f_idx(i)-1));
                valueArray(end+1) = -1;
            else
                if roi_outer_bp(r-1,c)
                    n = n+1;
                    b(i) = b(i)+dest_bp(r-1,c);
                end
            end
            if roi_mask(r,c+1),
                n = n+1;
                indexArray1(end+1) = i;
                indexArray2(end+1) = find(f_idx==(f_idx(i)+size(roi_mask,1)));
                valueArray(end+1) = -1;
            else
                if roi_outer_bp(r,c+1)
                    n = n+1;
                    b(i) = b(i)+dest_bp(r,c+1);
                end
            end
            if roi_mask(r,c-1),
                n = n+1;
                indexArray1(end+1) = i;
                indexArray2(end+1) = find(f_idx==(f_idx(i)-size(roi_mask,1)));
                valueArray(end+1) = -1;
            else
                if roi_outer_bp(r,c-1)
                    n = n+1;
                    b(i) = b(i)+dest_bp(r,c-1);
                end
            end
            if n
                indexArray1(end+1) = i;
                indexArray2(end+1) = i;
                valueArray(end+1) = n;
            end
        end

        A = sparse(indexArray1,indexArray2,valueArray,length(f_idx),length(f_idx));
        x = ones(length(f_idx),1);
    
        % jacobi method
        x_old = x;
        x_new = x_old;
        th0 = 0;
        while (1)
            sigma = A*x_old-diag(diag(A))*x_old;
            x_new = (b-sigma)./full(diag(A));
            th1 = max(abs(x_new - x_old))           % new threshold
            if max(abs(x_new - x_old)) < 0.01
                break;
            else
                th0 = max(abs(x_new - x_old));      % old threshold          
                x_old = x_new;
            end
        end
        
        result = zeros(size(roi_mask));
        result(f_idx) = x_new;

        f(:,:,chn) = result;
        
        tmp = dest_roi(:,:,chn);
        tmp(f_idx) = x_new;
        dest_roi(:,:,chn) = tmp;

    end
    
    I(ty-1:ty+hh,tx-1:tx+ww,1:3) = dest_roi;    % chane dest roi to calculated roi

    
    figure
    image(uint8(I));
    DEST = uint8(DEST);
    I = uint8(I);
    imwrite(DEST,'dest.tif');
    imwrite(I, 'result.tif');

return
