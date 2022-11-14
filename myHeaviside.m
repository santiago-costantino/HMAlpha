function y=myHeaviside(x)

% a simple version of the heaviside function, as it cannot be used by the Matlab
% compiler because it belongs to the symboic toolbox. Nice explanation.
% This simpler version assumes is only made to solve the problem of using
% it for getRetinaBM


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

% This function analyzes Hypotony Maculopathy from OCT images. Input is the 
% folder where images (tiff or png, only) are located and the x and y scales of the OCT machine used. 

% Detailed explanations of the algorithm and use are available in an article 
% published in the Journal of Glaucoma, 
% "Quantification of hypotony maculopathy using spectral-domain optical coherence tomography" by R. Dumas et al.

y(x>0)=1;
y(x==0)=0.5;
y(x<0)=0;

