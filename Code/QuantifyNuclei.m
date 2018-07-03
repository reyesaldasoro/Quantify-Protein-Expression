function [numNuclei,dataOut] = QuantifyNuclei(dataInName)

%function [numHorizonalNuclei,numVerticalNuclei,numIntenseNuclei,dataOut] = QuantifyNuclei(dataInName)

    %presentDir =pwd;
if nargin ==0
        %----- no data received, Open question dialog and pass to next section to analyse
        button                                  = questdlg('Please specify the location of Image data sets (one channel at a time)','Select Input','Multiple Files in a Folder','Cancel','Cancel');
        if strcmp(button(1),'M')
            %----- a folder, can be a folder with matlab files, or a folder with folders and tiff images
            [pathname]                          =  uigetdir('*.*','Please select folder where the images/data are located');
            if pathname~=  0
                % pass the pathname to same function to process
                dir0                            = pathname;
            else
                %disp('Folder not found');
                numHorizonalNuclei=[];
                return;
            end
        else
            numHorizonalNuclei=[];
            return;
        end
else
     dir0                                       = dataInName;    
end
%cd(dir0)

 
%clear data* thres*
dir1                                            = dir(strcat(dir0,'/*.tif'));
%dir1                                            = dir('*.tif');
numSubDir                                       = size(dir1,1);

firstSample                                     = imread(strcat(dir0,'/',(dir1(1).name)));
[rows,cols,levs]                                = size(firstSample);

disp(strcat(num2str(numSubDir),' images in folder, with dimensions: ',num2str(rows),' x  ',num2str(cols),' x  ',num2str(levs)));
 dataIn2(rows/4,cols/4,numSubDir)                   =0;
 dataIn3(rows/4,rows/4,numSubDir)                   =0;
 thresLevel2(numSubDir)                         =0;
% %counterDir=31;
for counterDir=1:numSubDir
     tempDir                                             = strcat(dir0,'/',dir1(counterDir).name);
     dataIn                                              = imread((tempDir));
     dataIn2(:,:,counterDir)                             = double(reduceu(dataIn,2));
     % Obtain thresholding levels and threshold
     thresLevel                                          = graythresh(dataIn(:));
     thresLevel2(counterDir)                             = thresLevel*255;
end
%% %
%threshold at the mean level of all the sets and find majority
%
%% Segment nuclei by intensity in two rounds
% there are some very intense nuclei and some dimmer ones, a single pass would merge
% some of the intense and discard some dim ones, do one loop for the very intense ones
% and then repeat for the dim ones, removing the bright intensities

thresLevelLow_Br                               =-2+ 0.75*mean(thresLevel2)+0.25*max(thresLevel2);
thresLevelHigh_Br                              = 2+ 0.25*mean(thresLevel2)+0.75*max(thresLevel2);
thresLevelLow_dim                              =-2+ 0.25*mean(thresLevel2)+0.75*min(thresLevel2);
thresLevelHigh_dim                             = 2+ 0.50*mean(thresLevel2)+0.50*min(thresLevel2);
%%
for counterDir=1:numSubDir
    %First pass for the bright ones
    [NucleiLowT,numNucleiLowT]              = bwlabeln(dataIn2(:,:,counterDir)>thresLevelLow_Br);
    AreaNucleiLowT                          = regionprops(NucleiLowT,'Area');
    LargeAreaNuclei                         = find([AreaNucleiLowT.Area]>3);
    
    % Only process if there are any nuclei there
    if numNucleiLowT>0

        NucleiHighT                         = (dataIn2(:,:,counterDir)>thresLevelHigh_Br);
        prodHighLowT                        = NucleiLowT.*NucleiHighT;
        NucleiLowT_2D                       = (sum(NucleiLowT,3))>0;
        clear rr cc
        [rr,cc]                             = find(NucleiLowT_2D);
        rr                                  = unique(rr);
        cc                                  = unique(cc);
        NucleiToKeep                        = unique(prodHighLowT(rr,cc,:));        
        dataIn3(:,:,counterDir)             = double(ismember(NucleiLowT,intersect(LargeAreaNuclei,  NucleiToKeep(2:end))));
    else

        dataIn3(1,1,counterDir)       = 0;
    end
%    %second pass for the dim ones
    dilatedBrightNuclei                     = 1-(imdilate(dataIn3(:,:,counterDir),ones(3)));
    try
    [NucleiLowT,numNucleiLowT]              = bwlabeln((dataIn2(:,:,counterDir).*(dilatedBrightNuclei))>thresLevelLow_dim);
    catch
        qqq=1;
    end
    AreaNucleiLowT                          = regionprops(NucleiLowT,'Area');
    LargeAreaNuclei                         = find([AreaNucleiLowT.Area]>3);
    % Only process if there are any nuclei there
    if numNucleiLowT>0
        NucleiHighT                         = ((dataIn2(:,:,counterDir).*(dilatedBrightNuclei))>thresLevelHigh_dim);
        prodHighLowT                        = NucleiLowT.*NucleiHighT;
        NucleiLowT_2D                       = (sum(NucleiLowT,3))>0;
        clear rr cc
        [rr,cc]                             = find(NucleiLowT_2D);
        rr                                  = unique(rr);
        cc                                  = unique(cc);
        NucleiToKeep                        = unique(prodHighLowT(rr,cc,:));        
        dataIn3(:,:,counterDir)             = dataIn3(:,:,counterDir)+double(ismember(NucleiLowT,intersect(LargeAreaNuclei,  NucleiToKeep(2:end))));
    else
        dataIn3(1,1,counterDir)       = 0;
    end
    
    %Do a small close to join possibly separated objects
    dataIn3(:,:,counterDir)                 = (imclose(dataIn3(:,:,counterDir),ones(2)));
    dataIn3(:,:,counterDir)                 = bwmorph(dataIn3(:,:,counterDir),'majority');
end
%%


% %% To get more uniform results find majority
% for counterDir=1:numSubDir
%     dataIn3A(:,:,counterDir)                              = bwmorph(dataIn3(:,:,counterDir),'majority');
% end
% %% To get a wider area with the low level nuclei do a closing
% closeElement = strel('disk',3);
% for counterDir=1:numSubDir
%     dataIn3B(:,:,counterDir)                              = imclose(dataIn3A(:,:,counterDir),closeElement);
% end
% dataIn3                                    = (imclose(dataIn3,strel('ball',3,3)));
%dataIn3                                     = dataIn3>0.5;
%% find the edges of the vessel


% for counterDir=1:numSubDir
%     tempEdges                               = max(dataIn3(:,:,counterDir));
% 
%     edgeVesselLeft(counterDir)              = find(tempEdges,1,'first');
%     edgeVesselRight(counterDir)             = find(tempEdges,1,'last');
% 
% end



%% 
% Finally remove all those 3D structures that are small, 15 voxels in total
[AllNuclei,numLargeNuclei]                  = bwlabeln(dataIn3(:,:,:));
AreaNuclei                                  = regionprops(AllNuclei,'Area','BoundingBox');
LargeAreaNuclei                             = find([AreaNuclei.Area]>15);
[LargeNuclei,numNuclei]                     = bwlabeln(ismember(AllNuclei,LargeAreaNuclei));
%AreaNuclei                                  = regionprops(LargeNuclei,'Area','BoundingBox');

dataOut = (LargeNuclei>0);

%% The lines below classify according to shape, but may not be necessary with the ELASTIN boundaries,
%% remove for the time being
% % Classify the nuclei into three groups: 
% %   1 Bright at the edges of the vessel
% %   2 Faint below the others perpendicular to the vessel
% %   3 faint at the centre parallel to the vessel
% 
% %%
% %propsLargeNuclei                            = regionprops(LargeNuclei(:,:,60:70),'BoundingBox');
% %%
% % Find average intensity of the cells use this in conjunction to position and orientation to discriminate
% 
% %redLargeNuclei                              = LargeNuclei(1:2:end,1:2:end,:);
% avIntensityNuclei(numLargeNuclei)           = 0;
% avOrientation(numLargeNuclei)               = 0;
% avEccentricity(numLargeNuclei)              = 0; 
% %%
% %thresholdArea                               = 25;
% for counterNuclei = 1:numLargeNuclei
%     
%     if (AreaNuclei(counterNuclei).Area)>thresholdArea
%         tempRows                                = max(1,floor(AreaNuclei(counterNuclei).BoundingBox(2))):min(rows/4,ceil(AreaNuclei(counterNuclei).BoundingBox(2)+AreaNuclei(counterNuclei).BoundingBox(5)));
%         tempCols                                = max(1,floor(AreaNuclei(counterNuclei).BoundingBox(1))):min(cols/4,ceil(AreaNuclei(counterNuclei).BoundingBox(1)+AreaNuclei(counterNuclei).BoundingBox(4)));
%         tempLevs                                = max(1,floor(AreaNuclei(counterNuclei).BoundingBox(3))):min(numSubDir,ceil(AreaNuclei(counterNuclei).BoundingBox(3)+AreaNuclei(counterNuclei).BoundingBox(6)));
%         try
%         %disp(counterNuclei)
%         tempNuclei                              = AllNuclei(tempRows,tempCols,tempLevs)==counterNuclei;
%         tempData                                = dataIn2(tempRows,tempCols,tempLevs);
%         tempIntensity                           = tempData(tempNuclei);
%         catch
%             q=1;
%         end
%         avIntensityNuclei(counterNuclei)        = mean(tempIntensity);
%         tempProps                               = regionprops(max(tempNuclei,[],3),'Orientation','Eccentricity');
%         avOrientation(counterNuclei)            = tempProps.Orientation;
%         avEccentricity(counterNuclei)           = tempProps.Eccentricity;
%         
%         avDistanceFromEdge(counterNuclei)       = min(abs(tempCols(1)-edgeVesselLeft(tempLevs(1))),abs(tempCols(end)-edgeVesselRight(tempLevs(1))));
%     else
%         avIntensityNuclei(counterNuclei)        = 0;
%         avOrientation(counterNuclei)            = 0;
%         avDistanceFromEdge(counterNuclei)       = 0;
%     end
% end
% 
% 
% %
% 
% %Distances From Edges
% distanceClose                               =  50/4; 
% distanceMed                                 =  80/4; 
% distanceFar                                 = 110/4; 
% distanceVFar                                = 120/4; 
% 
% 
% avAreaNuclei                                = ([AreaNuclei.Area])>thresholdArea;
% 
% %threshold to distinguish the very intense from the not so
% intensityThreshold                          = graythresh(uint8(avIntensityNuclei(avAreaNuclei)))*255;
% distanceFromEdge                            = 100;
% EccentricityThreshold                       = 0.9;
% 
% %First rule, *very-very* intense,  any shape any position
% VeryBright                                  = (avIntensityNuclei>(0.5*intensityThreshold+0.5*255));
% 
% VeryBright2                                  = (avDistanceFromEdge<distanceMed).*(avIntensityNuclei>(0.75*intensityThreshold+0.25*255));
% %Second rule, very intense,  closer to edges any shape but not horizontal
% BrightCloserToEdge                          = (avDistanceFromEdge<distanceClose).*(avIntensityNuclei>intensityThreshold).*((avOrientation<-45)+(avOrientation>45));
% 
% %Second rule,  intense,  close to edges not very elongated 
% BrightCloseToEdgeRoundish                   = (avDistanceFromEdge<distanceFar).*(avIntensityNuclei>intensityThreshold).*(avEccentricity<0.9);
% %third rule very elongated and very directional
% VeryElongatedHorizontal                     = (avOrientation>-15).*(avOrientation<15).*(avEccentricity>0.95);
% VeryElongatedVertical                       = ((avOrientation<-75)+(avOrientation>75)).*(avEccentricity>0.95);
% %fourth rule not so elongated, very directional and far from edges
% ElongatedHorizontal                         = (avOrientation>-15).*(avOrientation<15).*(avEccentricity>0.9).*(avDistanceFromEdge>70);
% ElongatedVertical                           = ((avOrientation<-75)+(avOrientation>75)).*(avEccentricity>0.9).*(avDistanceFromEdge>70);
% RoundishVerticalFarFromEdge                 = (avDistanceFromEdge>distanceVFar).*(avIntensityNuclei<intensityThreshold).*((avOrientation<-25)+(avOrientation>25));
% 
% % intenseAndCloseToEdge                       = (-avDistanceFromEdge+avIntensityNuclei)>50;
% % FaintFarFromEdge                            = (-avDistanceFromEdge+avIntensityNuclei)<50;
% % HorizontalOrientation                       = (avOrientation>-45).*(avOrientation<45);
% % VerticalOrientation                         = (avOrientation<-45)+(avOrientation>45);
% % 
% % veryIntenseNuclei                           = ismember (AllNuclei,find(intenseAndCloseToEdge.*avAreaNuclei));
% % %faintNuclei                                 = ismember (AllNuclei,find((-avDistanceFromEdge+avIntensityNuclei)<50));
% 
% %
% PeriferalNucleiID                       = (VeryBright+VeryBright2+BrightCloseToEdgeRoundish+BrightCloserToEdge)>0;
% HorizontalNucleiID                      = (VeryElongatedHorizontal+ElongatedHorizontal)>0;
% VerticalNucleiID                        = (1-PeriferalNucleiID).*(ElongatedVertical+VeryElongatedVertical+RoundishVerticalFarFromEdge)>0;
% 
% PeriferalNuclei                         = ismember (AllNuclei,find(PeriferalNucleiID.*avAreaNuclei));
% HorizontalNuclei                        = ismember (AllNuclei,find(HorizontalNucleiID.*avAreaNuclei));
% VerticalNuclei                          = ismember (AllNuclei,find(VerticalNucleiID.*avAreaNuclei));
% 
% 
% 
% 
% [HorizontalLab,numHorizonalNuclei ]         = bwlabeln( HorizontalNuclei);
% [VerticalLab,numVerticalNuclei ]            = bwlabeln( VerticalNuclei);
% [IntenseLab,numIntenseNuclei ]              = bwlabeln( PeriferalNuclei);
% 
% dataOut                                     =   PeriferalNuclei +2*HorizontalNuclei +3*VerticalNuclei;
% 

%%
% % find the nuclei of the cells by using a separate threshold for very high levels for the external nuclei, use
% % the hysteresis to segment the cells and then find the average intensity of the cells and separate into
% % those that are very uniformly intense and those that are lower intensity, also use the edge of the
% % vessel (sum(dataIn)) to determine those that are on the edges then use bwlabeln to find the individual
% % cells and regionprops to discard those that are too small, hopefully there will not be too large. Once
% % that has been established, get orientation to detect the ones that are parallel and perpendicular.
% 
% % Contar el numero de numero de nucleos del C2, distinguir entre los nucleos mas intensos y regulares y los mas alargados y menos intensos
% % 
% % Paralelos ? endotelio
% % Perpendiculars ? musculo liso
% % Circulares exteriores - fibroblastos
% % 
% % Sobre el canal rojo, contar la actina entre esos mantos
% % 
% % Azul nucleos dna
% % Rojo elastina
% 
%cd(presentDir);