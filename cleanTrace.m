function [traceOut, outX] = cleanTrace(traceIn,nRows)

    % Make the trace a column vector
    traceIn = traceIn(:);

    start = find(~isnan(traceIn),1,'first');
    fin   = find(~isnan(traceIn),1,'last');

    % outX    = (start:fin)'; changed by remy
    % traceIn = traceIn(start:fin);
    outX = (1:numel(traceIn))';
    nonanMask = ~isnan(traceIn);

    % Remove nans for analysis
    nonanX  = outX(nonanMask);
    traceIn = traceIn(nonanMask);

    % Get the medianfilter version of the trace as a reference to find jumps
    interpTrace = interpUnique(nonanX,traceIn,outX,'linear');
    yMed        = medfilt1(interpTrace,15);

    % Remove large deviations (jumps)
    dev       = abs(traceIn(:) - yMed(nonanMask));
    threshold = max(5 * median(dev),1);
    validMask = abs(traceIn(:) - yMed(nonanMask)) < threshold;

    nonanX  = nonanX(validMask);
    traceIn = traceIn(validMask);

    % Compute values that were NaN using linear interpolation
    traceOut = interpUnique(nonanX,traceIn,outX,'linear', 'extrap');

    % Coerce to image limits
    traceOut = max(1,min(nRows,traceOut));

    % Smooth the surface using robust-weighted 2nd-degree local-regression
    traceOut = smooth(traceOut,15,'rloess');
    traceOut = round(traceOut);

end
