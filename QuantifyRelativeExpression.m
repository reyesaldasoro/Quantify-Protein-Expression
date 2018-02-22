function [totExpression,locExpression,numFinalNuclei,finalData,ElastinVsCollagen,NucleiVsElastinBoundary] = QuantifyRelativeExpression(dataInName)

%presentDir =pwd;
if nargin ==0
    %----- no data received, Open question dialog and pass to next section to analyse
    button                                  = questdlg('Please specify the location of FOUR CHANNELS','Select Input','Multiple Files in a Folder','Cancel','Cancel');
    if strcmp(button(1),'M')
        %----- a folder, can be a folder with matlab files, or a folder with folders and tiff images
        [pathname]                          =  uigetdir('*.*','Please select folder where FOUR folders are located');
        if pathname~=  0
            % pass the pathname to same function to process
            dir0                            = pathname;
        else
            %disp('Folder not found');
            totExpression=[];locExpression=[];numFinalNuclei=[];finalData=[];
            return;
        end
    else
        totExpression=[];locExpression=[];numFinalNuclei=[];finalData=[];
        return;
    end
else
    dir0                                       = dataInName;
end
%%
dir00                                           = dir(strcat(dir0,'/channel*'));
numSubDir00                                      = size(dir00,1);

if (numSubDir00~=4)
    disp('The folders are not arranged in the correct way.') 
    disp('Please save them with the following names:') 
    disp('     channel0 - YELLOW = Phalloidin staining of ACTIN fibers') 
    disp('     channel1 - RED    = Alexa staining of ELASTIN') 
    disp('     channel2 - BLUE   = DAPI nuclei') 
    disp('     channel3 - GREEN  = Second harmonic imaging of COLLAGEN') 
    totExpression =[];
    return;
end
%%
%determine the directories of the channels
dir_ch00                                         = ((strcat(dir0,'/',dir00(1).name)));
dir_ch10                                         = ((strcat(dir0,'/',dir00(2).name)));
dir_ch20                                         = ((strcat(dir0,'/',dir00(3).name)));
dir_ch30                                         = ((strcat(dir0,'/',dir00(4).name)));

%%
% Read the folders 
dir_ch0                                         = dir((strcat(dir_ch00,'/*.tif')));
dir_ch1                                         = dir((strcat(dir_ch10,'/*.tif')));
dir_ch2                                         = dir((strcat(dir_ch20,'/*.tif')));
dir_ch3                                         = dir((strcat(dir_ch30,'/*.tif')));
%%
numSubDir                                       = size(dir_ch0,1);

firstSample                                     = imread(strcat(dir0,'/',dir00(1).name,'/',(dir_ch0(1).name)));
[rows,cols,levs]                                = size(firstSample);

%% Read red channel of the ELASTIN to get the boundaries of the channel
channel1_Elastin                                = zeros(rows/2,cols/2,numSubDir);
for counterDir=1:numSubDir
    %read the images
    nameImage_ch1                                       = strcat(dir0,'/',dir00(2).name,'/',dir_ch1(counterDir).name);
    currImage_ch1                                       = (imread((nameImage_ch1)));
    %subsample once to reduce the size of the data 

    channel1_Elastin(:,:,counterDir)                    = reduceu(currImage_ch1);

end

% The red channel is now a volumetric matrix, it is necessary to unwrap the vessel to calculate the
% extent of the surfaces described by the ELASTIN, loop over the rows subsampling
%
%subsample along the columns and rows
channel1_ElastinR                                       = 0.5*(channel1_Elastin(:,1:2:end,:) +channel1_Elastin(:,2:2:end,:));
channel1_ElastinC                                       = 0.5*(channel1_ElastinR(1:2:end,:,:) +channel1_ElastinR(2:2:end,:,:));

% Find the centre of rotation for the vessel

vessProjection                                          = (squeeze(mean(mean(channel1_ElastinC(:,:,:),1),3)));
averageVesselValue                                      = mean(channel1_ElastinC(:));
initValue                                               = find(vessProjection>averageVesselValue,1,'first');
finalValue                                              = find(vessProjection>averageVesselValue,1,'last');
centreRotation                                          = [size(channel1_ElastinC,3) round((initValue+finalValue)/2)];

% Loop over the axial projections of the vessel to unwrap and find the boundaries of the ELASTIN
% To find the two edges of the ELASTIN expression, it is necessary to read all slices, create
% Volume to then unwrap, find surfaces, find maxima, and return.
for counterRows=2:2:size(channel1_ElastinC,1)
    %disp(counterRows)
    transposedCh1               =  (squeeze(channel1_ElastinC(counterRows,:,:)))';
    % this converts the cartesian data into a polar one rotating around the defined origin
    [polarData]                 = cart2Polar(transposedCh1,[centreRotation(1) centreRotation(2)]);
    %smooth the data to avoid small peaks
    polarDataSm                 = imfilter(polarData,gaussF(5,5,1),'replicate');

    % now find the peaks along each angle of the rotation, only 0-180 as the other half is empty
    peakPosition                = zeros(180,2);
    for k=1:180
        [peaksVess,indPeaks]    = findpeaks(polarDataSm(:,k),'npeaks',6,'minpeakheight',3,'sortstr','descend','minpeakdistance',2);
        if numel(peaksVess)==1
            peakPosition(k,:)       = [indPeaks indPeaks+2];
        elseif numel(peaksVess)>1
            peakPosition(k,:)       = sort(indPeaks(1:2));
        end
    end
    %filter the positions to smooth the lines
    %peakPosition2 = round((peakPosition + [peakPosition(1,:); peakPosition(1:end-1,:)]+ [peakPosition(2:end,:); peakPosition(end,:)])/3);

    
    cartData2=zeros(size(transposedCh1));
    % Return the peaks to the cartesian Data
    for k=1:180
        tempAngleSin        = sin(pi*k/180);
        tempAngleCos        = cos(pi*k/180);
        % The maximum places it above zero and the minimum places it below rows/4,
        % cols/4
        rr                  = min(max(1,round(-peakPosition(k,1)*tempAngleSin+centreRotation(1))),rows/4);
        cc                  = min(max(1,round( peakPosition(k,1)*tempAngleCos+centreRotation(2))),cols/4);
        cartData2(rr,cc)    = 1;
        rr                  = min(max(1,round(-peakPosition(k,2)*tempAngleSin+centreRotation(1))),rows/4);
        cc                  = min(max(1,round( peakPosition(k,2)*tempAngleCos+centreRotation(2))),cols/4);
        

        %try
        cartData2(rr,cc)    = 1;
        %catch
        %    q=1;
        %end
    end

    % close the regions between the two surfaces
    cartData3 =(imclose(cartData2,strel('disk',16)));

    try
    channel1_Elastin_Seg(counterRows/2,:,:) = cartData3'; %#ok<AGROW>
    catch
        q=1;
    end
end
%%
% Loop over the sagital projection of the boundaries to smooth it
channel1_Elastin_Seg2 = zeros(size(channel1_Elastin_Seg));
for counterCols=1:size(channel1_Elastin_Seg,2)
    channel1_Elastin_Seg2(:,counterCols,:) = imerode(imclose(imopen(squeeze(channel1_Elastin_Seg(:,counterCols,:)),ones(3)),ones(3)),ones(3));
end
%% Process the DAPI blue channel separately to obtain the nuclei

[numNuclei,channel2_Nuclei] = QuantifyNuclei(dir_ch20);




%% Quantification of the Expression
% Now with the boundaries of the Elastin, loop over the rest of the channels 
% and quantify their expression within the boundaries and outside the boundaries

%% Quantify expresion of four channels
% This Loop  will read the image of the four folders, one at a time from each, quantify the expression
% of each channel and then return the values of each, in the process, two layers of Elastin will be
% determined to quantify the expression inside that region as well
[rows_Elastin,cols_Elastin,levs_Elastin]                = size(channel1_Elastin_Seg2);

%%
[xx,yy]             = meshgrid(linspace(1,cols_Elastin,cols),linspace(1,rows_Elastin,rows));
totExpression       = zeros(numSubDir,4);
locExpression       = zeros(numSubDir,4);
totExpression_ch0_3D = zeros(rows_Elastin,cols_Elastin/2,numSubDir);
totExpression_ch1_3D = zeros(rows_Elastin,cols_Elastin/2,numSubDir);
totExpression_ch2_3D = zeros(rows_Elastin,cols_Elastin/2,numSubDir);
totExpression_ch3_3D = zeros(rows_Elastin,cols_Elastin/2,numSubDir);

for counterDir=1:numSubDir
    %counterDir=1;
    %recover the current elastin boundaries
    %try
    currentElastinBoundaries                            = interp2(channel1_Elastin_Seg2(:,:,counterDir),xx,yy,'nearest');
    %catch
    %    q=1;
    %end
    %read the images
    nameImage_ch0                                       = strcat(dir0,'/',dir00(1).name,'/',dir_ch0(counterDir).name);
    nameImage_ch1                                       = strcat(dir0,'/',dir00(2).name,'/',dir_ch1(counterDir).name);
    nameImage_ch2                                       = strcat(dir0,'/',dir00(3).name,'/',dir_ch2(counterDir).name);
    nameImage_ch3                                       = strcat(dir0,'/',dir00(4).name,'/',dir_ch3(counterDir).name);
    currImage_ch0                                       = (imread((nameImage_ch0)));
    currImage_ch1                                       = (imread((nameImage_ch1)));
    currImage_ch2                                       = (imread((nameImage_ch2)));
    currImage_ch3                                       = (imread((nameImage_ch3)));
    % find the otsu threshold
    thresLevel_ch0                                      = 255*graythresh(currImage_ch0(:));
    thresLevel_ch1                                      = 255*graythresh(currImage_ch1(:));
    thresLevel_ch2                                      = 255*graythresh(currImage_ch2(:));
    thresLevel_ch3                                      = 255*graythresh(currImage_ch3(:));
    % threshold the image and apply a majority vote to quantify the expression
    totExpression_ch0                                      = bwmorph(currImage_ch0>thresLevel_ch0,'majority');
    totExpression_ch1                                      = bwmorph(currImage_ch1>thresLevel_ch1,'majority');
    totExpression_ch2                                      = bwmorph(currImage_ch2>thresLevel_ch2,'majority');
    totExpression_ch3                                      = bwmorph(currImage_ch3>thresLevel_ch3,'majority');
    % calculate the localised expression within the boundaries of the Elastin
    locExpression_ch0                                       = (totExpression_ch0.*currentElastinBoundaries);
    locExpression_ch1                                       = (totExpression_ch1.*currentElastinBoundaries);
    locExpression_ch2                                       = (totExpression_ch2.*currentElastinBoundaries);
    locExpression_ch3                                       = (totExpression_ch3.*currentElastinBoundaries);
   
    %Quantify expression per slice
    totExpression(counterDir,:)   = [sum(totExpression_ch0(:)) sum(totExpression_ch1(:)) sum(totExpression_ch2(:)) sum(totExpression_ch3(:))]; 
    locExpression(counterDir,:)   = [sum(locExpression_ch0(:)) sum(locExpression_ch1(:)) sum(locExpression_ch2(:)) sum(locExpression_ch3(:))];

    % save subsampled version for a final reconstruction
    
    totExpression_ch0_3D (:,:,counterDir)                   = totExpression_ch0(1:8:end,1:8:end);
    totExpression_ch1_3D (:,:,counterDir)                   = totExpression_ch1(1:8:end,1:8:end);
    totExpression_ch2_3D (:,:,counterDir)                   = totExpression_ch2(1:8:end,1:8:end);
    totExpression_ch3_3D (:,:,counterDir)                   = totExpression_ch3(1:8:end,1:8:end);
    
end

totExpression                       = sum(totExpression);
locExpression                       = sum(locExpression);
%%


channel2_Nuclei_Lab                 = bwlabeln(channel2_Nuclei(1:2:end,:,:));

%nucleiInsideElastin                 = channel2_Nuclei_Lab.*channel1_Elastin_Seg2;
%locNumNuclei                        = numel(unique(nucleiInsideElastin(:)))-1;

nucleiInsideElastin_3D              = channel2_Nuclei_Lab.*channel1_Elastin_Seg2;
nucleiOutsideElastin_3D             = channel2_Nuclei_Lab.*(1-channel1_Elastin_Seg2);
nucleiInsideElastin                 = unique(nucleiInsideElastin_3D(:));
nucleiOutsideElastin                = unique(nucleiOutsideElastin_3D(:));

nucleiInBoundaryElastin             = intersect(nucleiInsideElastin,nucleiOutsideElastin);

nucleiInsideElastin                 = setdiff(nucleiInsideElastin,nucleiInBoundaryElastin);
nucleiOutsideElastin                = setdiff(nucleiOutsideElastin,nucleiInBoundaryElastin);


%%
finalNuclei                         = double(ismember(channel2_Nuclei_Lab,nucleiInsideElastin))+...
                                      double(ismember(channel2_Nuclei_Lab,nucleiOutsideElastin))*2+...
                                      double(ismember(channel2_Nuclei_Lab,nucleiInBoundaryElastin(2:end)))*3;

                                  
numFinalNuclei                      = [numel(nucleiInsideElastin)  numel(nucleiOutsideElastin)  numel(nucleiInBoundaryElastin(2:end))];            
channel1_Elastin_Seg3               = imerode(channel1_Elastin_Seg2,ones(4,4,3));
channel1_Elastin_Seg4               = channel1_Elastin_Seg2& (1-channel1_Elastin_Seg3); 

%% Return a composite subsampled version of the channels
%To avoid overlapf of classes give priority:
%   1,2,3 - nuclei              Channel 2
%   4     - ELASTIN boundary    Channel 1
%   5     - Elastin             Channel 1
%   6     - Actin               Channel 0
%   7     - Collagen            Channel 3

finalData = finalNuclei(:,1:2:end,:) + 4*(channel1_Elastin_Seg4(:,1:2:end,:)).*(1-(finalNuclei(:,1:2:end,:)>0));

finalData = finalData + 5*(totExpression_ch1_3D).*(1-(finalData>0));
finalData = finalData + 6*(totExpression_ch0_3D).*(1-(finalData>0));
finalData = finalData + 7*(totExpression_ch3_3D).*(1-(finalData>0));
%% Create a composite image of Elastin vs Collagen

ElastinVsCollagen           =  (totExpression_ch1_3D)+2*(totExpression_ch3_3D);
NucleiVsElastinBoundary     = finalNuclei(:,1:2:end,:) + 4*(channel1_Elastin_Seg4(:,1:2:end,:)).*(1-(finalNuclei(:,1:2:end,:)>0));

%%

% kk=60;
% 
% chan0_SL =  squeeze(sum(totExpression_ch0_3D,1));
% chan1_SL =  squeeze(sum(totExpression_ch1_3D,1));
% chan2_SL =  squeeze(sum(finalNuclei,1));
% chan3_SL =  squeeze(sum(totExpression_ch3_3D,1));
% chan1_Bound_SL =  squeeze(sum(channel1_Elastin_Seg4(:,1:end,:),1));
% 
% clear combinedView1
% combinedView1(:,:,1) =chan1_SL'+0.49674*chan0_SL'+chan1_Bound_SL';
% combinedView1(:,:,2) =chan3_SL'+0.49674*chan0_SL'+chan1_Bound_SL';
% combinedView1(:,:,3) =6*chan2_SL'+chan1_Bound_SL';
% 
% 
% figure(1)
% combinedView1 = combinedView1/max(combinedView1(:));
% imagesc(combinedView1)

%%
% % call the other routines per channel
% %[numHorizonalNuclei,numVerticalNuclei,numIntenseNuclei,dataOut] = QuantifyNuclei('images/channel2');
% [AbsExpression0,RelExpression0,dataOut0] = QuantifyExpression('images/channel0');
% [AbsExpression1,RelExpression1,dataOut1] = QuantifyExpression('images/channel1');
% [AbsExpression2,RelExpression2,dataOut2] = QuantifyExpression('images/channel2');
% [AbsExpression3,RelExpression3,dataOut3] = QuantifyExpression('images/channel3');
% 
% %%
% data1 = dataOut(:,:,110);
% data2 = imopen(data1,strel('disk',3));
% data3 = imclose(data2,strel('line',20,90));
% %data4 = imopen(data3,strel('disk',1));
% %data4 = imclose(data2,strel('disk',32));
% 
% imagesc(data3+2*data2);
% %%
% k=420;
% imagesc(squeeze(dataOut0(k,:,:))'+2*squeeze(dataOut1(k,:,:))')
% %%
% 
% 
% 
