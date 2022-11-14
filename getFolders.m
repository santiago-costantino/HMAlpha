function outFolders=getFolders(inFolder)

outFolders=dir(inFolder);

isub = [outFolders(:).isdir];
outFolders=outFolders(isub);
outFolders(1:2)=[];