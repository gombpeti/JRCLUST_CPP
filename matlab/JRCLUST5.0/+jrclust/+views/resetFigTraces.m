function resetFigTraces(hFigTraces, rawTraces, hCfg)
    %RESETFIGTRACES
    if hCfg.nSegmentsTraces > 1
        tlim1 = ([0, size(rawTraces, 2)] + hFigTraces.figData.windowBounds(1) - 1) / hCfg.sampleRate;
        tlim1 = round(tlim1*1000)/1000;
    else
        tlim1 = hFigTraces.figData.windowBounds/hCfg.sampleRate;
    end

    shankOrder = hFigTraces.figData.shankOrder;
    shankVisible = hFigTraces.figData.shankVisible;
    shankAxes = hFigTraces.figData.shankAxes;
    visibleFullIdx = find(shankVisible);

    for v = 1:numel(visibleFullIdx)
        shankID = shankOrder(visibleFullIdx(v));
        nSitesShank = sum(hCfg.shankMap == shankID);
        axis(shankAxes(v), [tlim1, 0, nSitesShank + 1]);
    end
end
