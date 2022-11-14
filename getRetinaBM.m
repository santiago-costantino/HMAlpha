function [yRet, yBM] = getRetinaBM(im)

% Copyright (C) 2022, Rémy Dumas, Santiago Costantino 
% Hopital Maisonneuve-Rosemont, 
% Centre de Recherche
% www.biophotonics.ca
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

    averagingSizeX = 5;
    gradientLarge = 8;
    gradientSmall = 2;
    params.colSpacing = 2;
    params.weightThreshold = 0.1;
    
    % Preprocess the image
    imF = imfilter(im,fspecial('gaussian', [1 2*averagingSizeX], averagingSizeX));
    imStep = imfilter(imF, myHeaviside(-gradientLarge:gradientLarge)' - 0.5);
    imStepInv = - imStep;
    imStepInv(imStepInv < 0) = 0;
    imStepInv([1:2*gradientLarge, ...
               end-2*gradientLarge:end],:) = 0;
    
    % Reduce smoothing to get better precision
    imStepPrecise = - imfilter(imF, myHeaviside(-gradientSmall:gradientSmall)' - 0.5);
    imStepPrecise(imStepPrecise < 0) = 0;
    imStepPrecise([1:2*gradientLarge,...
                   end-2*gradientLarge:end],:) = 0;
    
    % Initialization
    yFirst  = NaN(1,size(imStep,2));
    wFirst  = NaN(1,size(imStep,2));
    ySecond = NaN(1,size(imStep,2));
    wSecond = NaN(1,size(imStep,2));
    yThird  = NaN(1,size(imStep,2));
    
    % Process each column
    for k = 1:size(imStep,2)
        % Find peaks on imStep image
        [pks,locs] = findpeaks(imStep(:,k),'SortStr','descend','MinPeakDistance',11);
        if numel(locs) < 2, continue, end
        locs = locs(1:2);
        pks = pks(1:2);
        
        % Reordering
        [locs, ix] = sort(locs);
        pks        = pks(ix);
        
        % Take the two best peaks
        yFirst(k)  = locs(1); % Retina
        wFirst(k)  = pks(1);
        ySecond(k) = locs(2); % top
        wSecond(k) = pks(2);
        
        % Find peaks on imStepInv image
        [pks,locs] = findpeaks(imStepInv(:,k),'SortStr','descend','MinPeakDistance',11);
        ix = find(locs > ySecond(k),1,'first'); % Choose the best one above ySecond
        if isempty(ix), continue, end
        yThird(k) = locs(ix); % bottom
    end
    
    % Remove outliers
    [yFirst, ~] = cleanTrace(yFirst,size(im,1));
    
    % Refine RPE
    RPEthickness       = yThird - ySecond;
    medianRPEthickness = nanmedian(RPEthickness);
    stdRPEthickness    = nanstd(RPEthickness);
    rpeMskWrong = RPEthickness > (medianRPEthickness + 3 * stdRPEthickness) |...
    RPEthickness < (medianRPEthickness - 3 * stdRPEthickness); % Find outliers
    
    ySecond(rpeMskWrong) = NaN;
    yThird(rpeMskWrong)  = NaN;
    [ySecond, ~] = cleanTrace(ySecond,size(im,1));
    [yThird, ~]  = cleanTrace(yThird,size(im,1));

    [xBM,yBM] = findRPEbottom(imF, imStepPrecise, ySecond, yThird, params);

%     Commented Javier's code end replaced by next line
%     CHcurve = floor(smooth(yBM(:),0.3,'rloess'));
%     CHcurve(end+1)= CHcurve(end);

    % Interpolate with NN the retina to have the same number of elements as BM
    yRet = interp1((1:numel(yFirst))',yFirst,(1:numel(yBM))','nearest','extrap');
    yRet(yRet<1)=1;
    yBM(yBM<1)=1;
end