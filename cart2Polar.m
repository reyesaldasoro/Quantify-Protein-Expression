function  [polarData,paddedCartData]=cart2Polar(cartData,centreRotation)
%function  [polarData,paddedCartData]=cart2Polar(cartData,centreRotation)

% transforms a 2D data set (an image typically) from cartesian coordinates to polar coordinates
% centred around a specific point of rotation, the default is the centre of the data
% the data is padded to leave the centreRotation at the centre of the padded image.
% The transformation is obtained by rotating the image and taking the "positive x-axis" values

[rows,cols,levs]=size(cartData);

if (levs~=1)||(rows==1)||(cols==1)
    %return as the data is not 2D
    disp('data is not 2D, please try again')
    polarData =[];
    return;
end

if ~exist('centreRotation','var')
    centreRotation(1)     = floor(rows/2)+1;
    centreRotation(2)     = floor(cols/2)+1;
end

% first, determine centre of rotation and pad so that the centre
% of rotation is the centre of a square

rows_above_centre   = centreRotation(1)-1;
rows_below_centre   = rows-centreRotation(1);
cols_left_centre    = centreRotation(2)-1;
cols_right_centre   = cols-centreRotation(2);
% the padded data will be a *square* of even dimensions
% in which the original data is accommodated
newDims = max(1+2*[rows_above_centre rows_below_centre cols_left_centre cols_right_centre]);
if ((newDims/2)==floor(newDims/2))
    newDims         = newDims+1;
end
centre_newDims      = 1+(newDims-1)/2;
paddedCartData      = zeros(newDims);
position_in_r   = centre_newDims-rows_above_centre:centre_newDims+rows_below_centre;
position_in_c   = centre_newDims-cols_left_centre:centre_newDims+cols_right_centre;


paddedCartData(position_in_r,position_in_c)=cartData;


%% fourth method, rotate the image around its centre to obtain the projection rays

stepFi                          = 6;
angleRange                      = 0:stepFi:359;
angleIndex                      = 1:length(angleRange);
polarData                       = zeros(centre_newDims,length(angleRange));
%polarData2                      = zeros(centre_newDims,360);
%polarData                       = zeros(centre_newDims,360);
for countAngle=angleIndex
    %rotated data
    dataRotated                 = imrotate(paddedCartData,-angleRange(countAngle),'crop');
    polarData(:,countAngle)     = dataRotated(centre_newDims,centre_newDims:end);
    %polarData(:,countAngle)     = dataRotated(centre_newDims,centre_newDims:-1:1);
end
if stepFi>1
[xx,yy]                         = meshgrid(linspace(1,length(angleRange),360),(1:centre_newDims));
polarData                       = interp2(polarData,xx,yy,'nearest');
end
%%    