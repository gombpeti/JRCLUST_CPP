function colors = selectionColors(clusterIds, nClusters)
    %SELECTIONCOLORS Each cluster's own colour, as the waveform view draws it.
    %
    %   Returns one RGB row per requested cluster id. The colours are the
    %   ones plotMeanWaveforms and plotFigTraces already use - row iCluster
    %   of lines(nClusters) - so a unit looks the same in every window.
    %
    %   This matters once more than two units can be selected. The views
    %   that compare a pair used fixed colours instead: black for the first
    %   selection and red for the second, from hCfg.colorMap. With four or
    %   five units that convention runs out, and worse, the colours did not
    %   agree with the waveform view, so identifying which curve belonged to
    %   which unit meant re-reading the legend in every figure.
    %
    %   Selection is still distinguishable in the waveform view because the
    %   overlay is drawn thicker than the background traces, not because it
    %   is a different hue.
    if nargin < 2 || isempty(nClusters)
        nClusters = max(clusterIds(:));
    end
    nClusters = max(nClusters, 1);

    cmap = lines(nClusters);

    % lines() returns exactly nClusters rows, but guard anyway: a stale
    % cluster id outside the current clustering should give a colour rather
    % than an index error while a figure is mid-refresh.
    idx = mod(clusterIds(:) - 1, size(cmap, 1)) + 1;
    colors = cmap(idx, :);
end
