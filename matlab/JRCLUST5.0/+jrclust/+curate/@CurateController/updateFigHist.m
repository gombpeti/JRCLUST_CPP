function updateFigHist(obj)
    %UPDATEFIGHIST Plot ISI histogram
    if isempty(obj.selected) || ~obj.hasFig('FigHist')
        return;
    end

    plotFigHist(obj.hFigs('FigHist'), obj.hClust, obj.hCfg, obj.selected);
end

%% LOCAL FUNCTIONS
function hFigHist = plotFigHist(hFigHist, hClust, hCfg, selected)
    %PLOTFIGHIST Plot ISI histogram
    % One histogram per selected unit, each in that unit's own cluster
    % colour - the same colour the waveform view uses - so a curve can be
    % matched to its unit without a legend. Previously only two were drawn,
    % in fixed black and red that corresponded to nothing else on screen.
    nSel = numel(selected);
    colors = jrclust.views.selectionColors(selected, hClust.nClusters);

    nBinsHist = 50; %TODO: put this in param file (maybe)
    XData = logspace(-1, 4, nBinsHist); % logarithmic scale from 0.1ms to 10s

    % draw the plot
    if ~hFigHist.hasAxes('default')
        hFigHist.addAxes('default');
        hFigHist.axApply('default', @set, 'XLim', [min(XData) max(XData)], 'XScale', 'log'); % ms
        hFigHist.axApply('default', @grid, 'on');
        hFigHist.axApply('default', @xlabel, 'ISI (ms)');
        hFigHist.axApply('default', @ylabel, 'Prob. Density');
    end

    % drop curves left by a larger previous selection
    k = nSel + 1;
    while hFigHist.hasPlot(sprintf('hPlot%d', k))
        hFigHist.rmPlot(sprintf('hPlot%d', k));
        k = k + 1;
    end

    for iSel = 1:nSel
        key = sprintf('hPlot%d', iSel);
        YData = getISIHistogram(selected(iSel), XData, hClust, hCfg);
        if ~hFigHist.hasPlot(key)
            hFigHist.addPlot(key, @stairs, nan, nan, 'Color', colors(iSel, :), 'LineWidth', 1.2);
        else
            hFigHist.plotApply(key, @set, 'Color', colors(iSel, :));
        end
        hFigHist.updatePlot(key, XData, YData);
    end

    if nSel == 1
        hFigHist.axApply('default', @title, sprintf('Unit %d', selected(1)));
    else
        hFigHist.axApply('default', @title, ...
            sprintf('Units %s (coloured as in the waveform view)', ...
                    strjoin(arrayfun(@num2str, selected, 'UniformOutput', 0), ', ')));
    end
end

function YData = getISIHistogram(iCluster, XData, hClust, hCfg)
    %GETISIHISTOGRAM Get a histogram of ISI values
    %   TODO: replace hist call with histogram
    clusterTimes = double(hClust.spikeTimes(hClust.spikesByCluster{iCluster}))/hCfg.sampleRate;
    YData = hist(diff(clusterTimes)*1000, XData);
    YData(end) = 0;
    YData = YData./sum(YData); % normalize
end
