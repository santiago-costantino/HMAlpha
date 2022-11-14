function [bumpProps, resImg] = getBumpProps(x1,x2,y,yPoly)

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

    bumpProps = [];
    
    % Compute bump area
    area = sum(abs(y(x1:x2) - yPoly(x1:x2)));
    y1 = y(x1); y2 = y(x2);

    % Change coordinate basis
    alpha = atan((y2-y1)/(x2-x1)); % Calculate rotation angle
    M = [cos(alpha) sin(alpha) 0;...
        -sin(alpha) cos(alpha) 0;...
                 0          0  1]; % Define tranformation matrix
    T = maketform('affine', M); % Create tform object

    % Apply tranformation to the trace and its associated polynomial
    [x_new, y_new] = tforminv(T,(x1:x2)',y(x1:x2));
    [x_new, I] = unique(x_new);
    y_new = y_new(I);
    [xPoly_new, yPoly_new] = tforminv(T,(x1:x2)',yPoly(x1:x2));
    [xPoly_new, I] = unique(xPoly_new);
    yPoly_new = yPoly_new(I);
    
    % Make the X reference (t) the same for both traces
    t = linspace(x_new(1), x_new(end), 100);
    v = interp1(x_new, y_new, t);
    vPoly = interp1(xPoly_new, yPoly_new, t);
    % v and vPoly are now comparable
    [~, imax] = max(abs(v - vPoly));

    % Assign properties
    bumpProps.area = area;
%     bumpProps.area2 = sqrt(nansum((v - vPoly).^2));
    [bumpProps.xw1, bumpProps.yw1] = tformfwd(T,t(1),v(1)); % Left point
    [bumpProps.xw2, bumpProps.yw2] = tformfwd(T,t(end),v(end)); % Right point
    [bumpProps.xh1, bumpProps.yh1] = tformfwd(T,t(imax),v(imax)); % Down point
    [bumpProps.xh2, bumpProps.yh2] = tformfwd(T,t(imax),vPoly(imax)); % Up point
    bumpProps = struct2table(bumpProps);
end
