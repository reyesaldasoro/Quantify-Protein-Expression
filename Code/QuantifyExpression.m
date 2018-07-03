function [AbsExpression,RelExpression,dataIn3,dataIn2] = QuantifyExpression(dataInName)

presentDir =pwd;
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
            dataIn=[];handles=[];
            return;
        end
    else
        dataIn=[];handles=[];
        return;
    end
else
    dir0                                       = dataInName;
end
%cd(dir0)


clear data* thres*
dir1                                            = dir(strcat(dir0,'/*.tif'));
%dir1                                            = dir('*.tif');
numSubDir                                       = size(dir1,1);

firstSample                                     = imread(strcat(dir0,'/',(dir1(1).name)));
[rows,cols,levs]                                = size(firstSample);

disp(strcat(num2str(numSubDir),' images in folder, with dimensions: ',num2str(rows),' x  ',num2str(cols),' x  ',num2str(levs)));
dataIn2(rows,cols,numSubDir)                   =0;
dataIn3(1024,1024,numSubDir)                   =0;
thresLevel2(numSubDir)                         =0;
% %counterDir=31;
for counterDir=1:numSubDir
    tempDir                                             = strcat(dir0,'/',dir1(counterDir).name);
    dataIn                                              = imread((tempDir));
    dataIn2(:,:,counterDir)                             = double(dataIn);
    % Obtain thresholding levels and threshold
    thresLevel                                          = graythresh(dataIn(:));
    thresLevel2(counterDir)                             = thresLevel*255;
end
% %%
%threshold at the mean level of all the sets and find majority
thresLevel3                         =mean(thresLevel2);
%
for counterDir=1:numSubDir
    dataIn3(:,:,counterDir)                             = bwmorph(dataIn2(:,:,counterDir)>thresLevel3,'majority');
end
%
% To get more uniform results find majority again%
% for counterDir=1:numSubDir
%      dataIn3(:,:,counterDir)                              = bwmorph(dataIn3(:,:,counterDir),'majority');
% end
%
%
AbsExpression                           = sum(dataIn(:));
RelExpression                           = sum(dataIn(:))/rows/cols/numSubDir;
%cd(presentDir);
%
%
%
% %% Processing of the DAPI channel to distinguish the nuclei of cells
% load('/Users/ccr22/Academic/work/microscopicCells/LuisMartinezLemus/channel0.mat')
%
% %%
%
% imagesc(dataIn(:,:,45))
% %%
% NucleiCells     = bwlabeln(dataIn);
%
%
% %%
%
% surfdat(NucleiCells(1:4:end,1:4:end,50:80),'all')
%
%
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
%
% %mkdir(dir2)
%
% %vesselsMetrics(numSubDir,1)        = 0;
% %histsLengths2(numSubDir,1)          = 0;
% %vesselsMetrics(numSubDir,23)        = 0;
% %
% %resultsFluoresence(numSubDir,11)    =0;
% %jet5=jet;
%
%
