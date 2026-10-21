function updateFigISI(obj)
    %UPDATEFIGISI Plot return map
    if isempty(obj.selected) || ~obj.hasFig('FigISI')
        return;
    end

    nSpikes = cellfun(@(c) numel(c), obj.hClust.spikesByCluster(obj.selected));
    if all(nSpikes > 1)
        plotFigISI(obj.hFigs('FigISI'), obj.hClust, obj.hCfg, obj.selected);
    else
        clf(obj.hFigs('FigISI'));
    end
end

%% LOCAL FUNCTIONS
function hFigISI = plotFigISI(hFigISI, hClust, hCfg, selected)
    %PLOTFIGISI Plot return map (ISI vs. previous ISI) for every selection
    %
    %   One series per selected unit, each in that unit's own cluster colour
    %   - the same colour the waveform view uses - so the series can be
    %   matched to a unit without a legend. Previously only two units were
    %   drawn, in fixed black and red that did not correspond to anything
    %   else on screen.
    nSel = numel(selected);
    colors = jrclust.views.selectionColors(selected, hClust.nClusters);

    if ~hFigISI.hasAxes('default')
        hFigISI.addAxes('default');
        hFigISI.axApply('default', @set, 'XScale','log', 'YScale','log');
        hFigISI.axApply('default', @xlabel, 'ISI_{k} (ms)');
        hFigISI.axApply('default', @ylabel, 'ISI_{k+1} (ms)');
        hFigISI.axApply('default', @axis, [1 10000 1 10000]);
        hFigISI.axApply('default', @grid, 'on');

        % refractory period reference lines
        hFigISI.addPlot('refracLine1', @line, hFigISI.axApply('default', @get, 'XLim'), hCfg.refracInt*[1 1], 'Color', [1 0 0]);
        hFigISI.addPlot('refracLine2', @line, hCfg.refracInt*[1 1], hFigISI.axApply('default', @get, 'YLim'), 'Color', [1 0 0]);
    end

    % Drop series left behind by a larger previous selection.
    k = nSel + 1;
    while hFigISI.hasPlot(sprintf('foreground%d', k))
        hFigISI.rmPlot(sprintf('foreground%d', k));
        k = k + 1;
    end
    % the original two keys, no longer used
    if hFigISI.hasPlot('foreground'),  hFigISI.rmPlot('foreground');  end
    if hFigISI.hasPlot('foreground2'), hFigISI.rmPlot('foreground2'); end

    for iSel = 1:nSel
        key = sprintf('foreground%d', iSel);
        [isiK, isiK1] = getReturnMap(selected(iSel), hClust, hCfg);
        if ~hFigISI.hasPlot(key)
            hFigISI.addPlot(key, nan, nan, 'Color', colors(iSel, :), ...
                            'Marker', 'o', 'LineStyle', 'none');
        else
            hFigISI.plotApply(key, @set, 'Color', colors(iSel, :));
        end
        hFigISI.updatePlot(key, isiK, isiK1);
    end

    if nSel == 1
        hFigISI.axApply('default', @title, sprintf('Unit %d', selected(1)));
    else
        hFigISI.axApply('default', @title, ...
            sprintf('Units %s (coloured as in the waveform view)', ...
                    strjoin(arrayfun(@num2str, selected, 'UniformOutput', 0), ', ')));
    end
end

function [isiK, isiK1] = getReturnMap(iCluster, hClust, hCfg)
    %GETRETURNMAP subset
    clusterTimes = double(hClust.spikeTimes(hClust.spikesByCluster{iCluster}))/hCfg.sampleRate;
    if numel(clusterTimes) < 3
        isiK = 0;
        isiK1 = 0;
        return;
    end
    clusterISI = diff(clusterTimes * 1000); % in msec

    subset = randperm(numel(clusterISI) - 1, min(hCfg.nSpikesFigISI, numel(clusterISI) - 1));
    isiK = clusterISI(subset);
    isiK1 = clusterISI(subset+1);
end
