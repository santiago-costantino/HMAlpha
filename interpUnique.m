function outY = interpUnique(inX,inY,outX,varargin)
% Interpolates by removing repetitions in the X coordinates, as it is
% forbidden by interp1.
% NOTE: This may not be the best approach. The unique values are chosen 
% arbitrary and they should be chosen acording to the neighbouring values.

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


     [auxX, auxY] = getUniqueX(inX,inY);
     outY = interp1(auxX,auxY,outX,varargin{:});

end