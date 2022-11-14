function [bumpsProps, y, yPoly] = getBumps(im, varargin)

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

    % Read parameters
    p = inputParser;
    p.addParameter('Boundary', 'BM');
    parse(p, varargin{:});
    [yILM, yBM] = getRetinaBM(double(im)); % Segment ILM and BM
    
    if isequal(p.Results.Boundary, 'BM')
        y = yBM;
        order = 4;
    elseif isequal(p.Results.Boundary, 'Retina')
        y = yILM;
        order = 12;
    end
    
    y = fillmissing(y, 'nearest');
    X = (1:numel(y))';
    [p,S,mu] = polyfit(X, y, order); % Fit polynomial
    yPoly = polyval(p,X,S,mu); % Get polynomial trace trace

    yPoly = fillmissing(yPoly, 'nearest');
    yPoly(yPoly<1) = 1; % Polynomial shouldn't be outside the image area

    % Find intersections between the trace and its associated polynomial
    X_inter = find(diff(sign(y - yPoly)));
    
    % Get segments (one segment, one fold)
    couplesX = table(X_inter(1:end-1), X_inter(2:end),...
                    'VariableNames', {'X1','X2'});
                
    % Get geometrical properties for each fold
    bumpsProps = rowfun(@(x1,x2) getBumpProps(x1,x2,y,yPoly), couplesX).Var1;
end

