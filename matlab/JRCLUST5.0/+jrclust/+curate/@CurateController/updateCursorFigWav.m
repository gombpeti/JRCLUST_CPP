function updateCursorFigWav(obj)
    %UPDATECURSORFIGWAV Plot mean waveforms on top of selected cluster(s)
    if isempty(obj.selected) || ~obj.hasFig('FigWav')
        return;
    end
    hFigWav = obj.hFigs('FigWav');

    % Overlay every selected unit, not just the first two.
    %
    % Each is drawn in its OWN cluster colour - row iCluster of
    % lines(nClusters), the same colour plotMeanWaveforms gives it - so a
    % unit looks the same here, in the ISI plot, the amplitude histogram and
    % the feature views. Previously the overlay was black for the first
    % selection and red for the second, which ran out at two units and did
    % not agree with the colours underneath.
    %
    % The overlay still stands out because it is drawn at twice the line
    % width, not because of its hue.
    nSel = numel(obj.selected);
    colors = jrclust.views.selectionColors(obj.selected, obj.hClust.nClusters);

    % Remove overlays left by a larger previous selection.
    k = nSel + 1;
    while hFigWav.hasPlot(sprintf('selected%d', k))
        hFigWav.rmPlot(sprintf('selected%d', k));
        k = k + 1;
    end

    for iSel = 1:nSel
        plotSelectedMeans(hFigWav, obj.hClust, obj.selected(iSel), ...
                          sprintf('selected%d', iSel), obj.maxAmp, obj.hCfg, ...
                          obj.spatial_idx, colors(iSel, :));
    end
end

%% LOCAL FUNCTIONS
function iCluster = plotSelectedMeans(hFigWav, hClust, iCluster, plotKey, maxAmp, hCfg, spatial_idx, color)
    %PLOTSELECTEDMEANS Plot an overlay on selected cluster's mean traces
    showSubset = hFigWav.figData.showSubset;

    if ~hFigWav.hasPlot(plotKey)
        hFigWav.addPlot(plotKey, nan, nan, 'Color', color, 'LineWidth', 2);
    else
        % the same slot may now hold a different unit, so recolour it
        hFigWav.plotApply(plotKey, @set, 'Color', color);
    end

    if hCfg.showRaw
        meanWf = hClust.meanWfGlobalRaw(:, :, iCluster);
    else
        meanWf = hClust.meanWfGlobal(:, :, iCluster);
    end

    hFigWav.multiplot(plotKey, maxAmp, jrclust.views.getXRange(find(iCluster == showSubset), size(meanWf, 1), hCfg), meanWf(:,spatial_idx));
end
