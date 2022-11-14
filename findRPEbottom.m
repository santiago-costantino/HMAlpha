function [xOut,traceOut] = findRPEbottom(imOri,invGrad,top, bottom, params)
% Compute bottom of RPE layer
% Parameters :
%     - imOri : original image (with gaussian filter)
%     - invGrad : image with gradient filter and inversed values
%     - top & bottom : layers found on imStep and imStepInv images
%     - params : general parameters of the graph-search method

sz = size(invGrad);
msk   = false(sz);
nodes = false(sz);

% Find the start and end of all the traces together
start = find(~isnan(top) & ~isnan(bottom),1,'first');
fin   = find(~isnan(top) & ~isnan(bottom),1,'last');
outX    = (start:fin)';
top(1:start) = top(start);
top(fin:end) = top(fin);
bottom(1:start) = bottom(start);
bottom(fin:end) = bottom(fin);

top(top<1) = 1;
bottom(bottom<1) = 1;
top(top>sz(1)) = sz(1);
bottom(bottom>sz(1)) = sz(1);
% Computes RPE msk
for k = 1:numel(top)
      msk(top(k):bottom(k),k) = true;
end

wInt = imOri;
wInt = (wInt - min(wInt(msk))) / (max(wInt(msk)) - min(wInt(msk))); 
wInt(~msk) = false;

% Find nodes
for x = 1:params.colSpacing:numel(top)
      vals = wInt(top(x):bottom(x),x);
      if isempty(vals)
          continue
      end
      [~,ix] = sort(vals,'descend');
      ix = top(x) + ix(1:min(3,numel(ix))) - 1;
      nodes(ix,x) = true;
end


[I,J] = find(nodes);
X = unique(J);
Y = nan(size(X));
for i=1:numel(X)
x = X(i);
Y(i) = mean(I(J==x));
end
yRPE = interpUnique(X,Y,outX,'linear','extrap');
[yRPE,~] = cleanTrace(yRPE,sz(1));

% Estimates the bottom edge of the RPE
rpeThickness = nanmedian(bottom - top);

botLim = round(min(sz(1), yRPE + rpeThickness));
msk(:) = false;
nodes(:) = false;

for k = 1:numel(yRPE)
      msk(yRPE(k):botLim(k),k) = true;
end

wGrad = invGrad;
wGrad = (wGrad - min(wGrad(msk))) / (max(wGrad(msk)) - min(wGrad(msk))); 
wGrad(~msk) = false;
% Find nodes on invGrad image
for x = 1:params.colSpacing:numel(yRPE)
      vals = wGrad(yRPE(x):botLim(x),x);
      [~,ix] = sort(vals,'descend');
      ix = yRPE(x) + ix(1:3) - 1;
      nodes(ix,x) = true;
end
% botGraph = GraphIm(nodes, wGrad); % Create GraphIm
% botGraph.getMinPath; % Find path through nodes
% 
% xBot = outX;
% yBot = interpUnique(botGraph.minPath.x,botGraph.minPath.y,xBot,'linear'); % Interpolation
% [yBot, xBot]  = cleanTrace(yBot,sz(1)); % Smoothing
% delete(botGraph);

[I,J] = find(nodes);
X = unique(J);
Y = nan(size(X));
for i=1:numel(X)
x = X(i);
Y(i) = mean(I(J==x));
end
yBot = interpUnique(X,Y,outX,'linear','extrap');
[yBot,xBot] = cleanTrace(yBot,sz(1));

% Estimated distance to BM
estDistBM = nanmedian(yBot - yRPE(xBot));
xOut = 1:numel(yRPE);
% traceOut = yRPE + estDistBM;
traceOut = yBot;

end
