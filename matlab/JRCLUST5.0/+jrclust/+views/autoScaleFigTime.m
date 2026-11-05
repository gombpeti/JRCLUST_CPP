function autoScaleFigTime(hFigTime, hClust, selected, iSite, channel_idx)
    %AUTOSCALEFIGTIME Automatically scale features in FigTime
    hCfg = hClust.hCfg;
    autoscalePct = hCfg.getOr('autoscalePct', 99.5)/100;

    % compute the autoscalePct'th quantile of every selected unit's display
    % features on iSite and set ylimits to the largest, so no unit is
    % clipped. Only called right after plotFigTime, so the features are read
    % back from the series it just drew rather than computed a second time.
    nSel = numel(selected);
    timeData = cell(1, nSel);
    for iSel = 1:nSel
        if iSel == 1
            key = 'foreground';
        else
            key = sprintf('foreground%d', iSel);
        end
        yData = hFigTime.plotApply(key, @get, 'YData');
        timeData{iSel} = yData(~isnan(yData));
    end

    dispFeaturesTime = timeData{1};
    if numel(dispFeaturesTime)<2 || ~any(dispFeaturesTime)
        timeData{1} = getFigTimeFeatures(hClust, iSite,[],channel_idx);
    end

    if all(cellfun(@isempty, timeData))
        return;
    end

    timeScale = max(cellfun(@(x) quantile(x(:), autoscalePct), timeData));
    if isnan(timeScale)
        timeScale = 1;
    end
    hFigTime.axApply('default', @set, 'YLim', [0, 1]*timeScale);
    hFigTime.axApply('histogram',@set,'YLim',[0,1]*timeScale);
    imrectSetPosition(hFigTime, 'hRect', [], [0, timeScale]);
end