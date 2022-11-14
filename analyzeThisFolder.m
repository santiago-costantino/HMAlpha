function [alphaVolumeBM, alphaVolumeRet]=analyzeThisFolder(varargin)

if nargin>0
    folderPath=varargin{1};
else
    folderPath='/Users/santiago/Dropbox (Biophotonics)/Projects/Remy/Hypotony/Data/hypotony_images/Marchand HM 2019-3-14';
end

imagePaths=dir(fullfile(folderPath, '*.tif'));
analyzeBscans(imagePaths);
[alphaVolumeBM, alphaVolumeRet]=analyzeVolume(imagePaths);



