function [wholeTable_bm, wholeTable_ret]=analyzeBscans(varargin)



% quantify Hypotony Maculopathy on sigle B-Scans
% image paths is an array pointing a single OCT images. The output consists
% of two tables containing the info of every fold found, both for the
% retina and the BM. Results are saved as mat and xls files

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


imagePaths=varargin{1};
if nargin==2
    scaleX=varargin{2}(1)
    scaleY=varargin{2}(2);
else
    scaleX = 5.7;
    scaleY = 3.9;
end


% paths point are the structures out di the function DIR, with fallged
% images excluded

%% settings
visualDir=strcat(imagePaths(1).folder, filesep,'Visualization')
if ~isfolder(visualDir)
    mkdir(visualDir);
end

imageOne = imread(fullfile(imagePaths(1).folder, imagePaths(1).name));
if isfile(fullfile(imagePaths(1).folder, 'ROI.mat'))
    load(fullfile(imagePaths(1).folder, 'ROI.mat'));
else
    ROI=[1, 1, size(imageOne, 2)-1, size(imageOne, 1)-1];
end



% Read crop information
xmin = ROI(1);
xmax = ROI(1) + ROI(3);
ymin = ROI(2);
ymax = ROI(2) + ROI(4);

    % Initialization
    wholeTable_ret = [];
    wholeTable_bm = [];
    for k=1:size(imagePaths, 1)
%         try
            % Load image
            imgPath = fullfile(imagePaths(k).folder, imagePaths(k).name);
            im = imread(imgPath);
            im = im(ymin:ymax,xmin:xmax,1); % Cropping

            % Find and analyze ILM folds
            [bumpsProps_ret_tmp,yILM,yPolyILM] = getBumps(im, 'Boundary', 'Retina');
            N_folds = height(bumpsProps_ret_tmp);
            bumpsProps_ret.id_scan = repmat(k, N_folds, 1);
            bumpsProps_ret.file = repmat({imgPath}, N_folds, 1);
            bumpsProps_ret.imWidth = repmat(size(im,2), N_folds, 1)*scaleX;
            bumpsProps_ret.width = sqrt((scaleX*(bumpsProps_ret_tmp.xw2 - bumpsProps_ret_tmp.xw1)).^2 +...
                                       (scaleY*(bumpsProps_ret_tmp.yw2 - bumpsProps_ret_tmp.yw1)).^2);
            bumpsProps_ret.height = sqrt((scaleX*(bumpsProps_ret_tmp.xh2 - bumpsProps_ret_tmp.xh1)).^2 +...
                                        (scaleY*(bumpsProps_ret_tmp.yh2 - bumpsProps_ret_tmp.yh1)).^2);
            bumpsProps_ret.sharpness = bumpsProps_ret.height./bumpsProps_ret.width;
            bumpsProps_ret.area = bumpsProps_ret_tmp.area*scaleX*scaleY;
%             bumpsProps_ret.area2 = bumpsProps_ret_tmp.area2*scaleY*sqrt(scaleX);
            
            % Find and analyze BM folds
            [bumpsProps_bm_tmp,yBM,yPolyBM] = getBumps(im, 'Boundary', 'BM');
            N_folds = height(bumpsProps_bm_tmp);
            bumpsProps_bm.id_scan = repmat(k, N_folds, 1);
            bumpsProps_bm.file = repmat({imgPath}, N_folds, 1);
            bumpsProps_bm.imWidth = repmat(size(im,2), N_folds, 1)*scaleX;
            bumpsProps_bm.width = sqrt((scaleX*(bumpsProps_bm_tmp.xw2 - bumpsProps_bm_tmp.xw1)).^2 +...
                                       (scaleY*(bumpsProps_bm_tmp.yw2 - bumpsProps_bm_tmp.yw1)).^2);
            bumpsProps_bm.height = sqrt((scaleX*(bumpsProps_bm_tmp.xh2 - bumpsProps_bm_tmp.xh1)).^2 +...
                                        (scaleY*(bumpsProps_bm_tmp.yh2 - bumpsProps_bm_tmp.yh1)).^2);
            bumpsProps_bm.sharpness = bumpsProps_bm.height./bumpsProps_bm.width;
            bumpsProps_bm.area = bumpsProps_bm_tmp.area*scaleX*scaleY;
%             bumpsProps_bm.area2 = bumpsProps_bm_tmp.area2*scaleY*sqrt(scaleX);
            
            % Concatenation
            wholeTable_ret = [wholeTable_ret;struct2table(bumpsProps_ret)];
            wholeTable_bm = [wholeTable_bm;struct2table(bumpsProps_bm)];
            
            % Visualization
            imshow(im); hold on; % Display original image (cropped)
            hILM = plot(yILM, 'LineWidth',2); % Plot boundaries and polynomial approximations
            hPolyILM = plot(yPolyILM, 'LineWidth',2);
            hBM = plot(yBM, 'LineWidth',2);
            hPolyBM = plot(yPolyBM, 'LineWidth',2);
            X_w = [bumpsProps_ret_tmp.xw1' bumpsProps_bm_tmp.xw1';...
                   bumpsProps_ret_tmp.xw2' bumpsProps_bm_tmp.xw2'];
            Y_w = [bumpsProps_ret_tmp.yw1' bumpsProps_bm_tmp.yw1';...
                   bumpsProps_ret_tmp.yw2' bumpsProps_bm_tmp.yw2'];
            X_h = [bumpsProps_ret_tmp.xh1' bumpsProps_bm_tmp.xh1';...
                   bumpsProps_ret_tmp.xh2' bumpsProps_bm_tmp.xh2'];
            Y_h = [bumpsProps_ret_tmp.yh1' bumpsProps_bm_tmp.yh1';...
                   bumpsProps_ret_tmp.yh2' bumpsProps_bm_tmp.yh2'];
            hH = plot(X_h, Y_h, 'm', 'LineWidth',2); % Display the detected fold heights
            hold off;
            exportgraphics(gca, fullfile(visualDir, sprintf('bscan_%03d.png', k))); % Save png
%             savefig(f, fullfile(visualDir, sprintf('bscan_%03d.fig', k))); % Save fig
%         catch
%             fprintf("Error on volume %d - bscan %d\n", k);
%         end
    end

    % Save
    save(fullfile(visualDir, 'regionProps_Ret.mat'), 'wholeTable_ret');
    save(fullfile(visualDir, 'regionProps_BM.mat'), 'wholeTable_bm');
    
    writetable(wholeTable_ret, fullfile(visualDir, 'regionProps_Ret.xls'));
    writetable(wholeTable_ret, fullfile(visualDir, 'regionProps_BM.xls'));


