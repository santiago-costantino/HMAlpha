clear
% runs the analysis with no GUI for a series of folders. It assumes ROI and
% flagged files are there

baseFolder='/Users/santiago/Dropbox (Biophotonics)/Projects/Remy/Hypotony/Data/hypotony_images';
volumeFolders=getFolders(baseFolder);

for p=1:size(volumeFolders, 1)
    [alphaVolumeBM(p), alphaVolumeRet(p)]=analyzeThisFolder(fullfile(volumeFolders(p).folder, volumeFolders(p).name));
end

