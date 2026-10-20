function beginSession(obj)
    %BEGINSESSION Start curating clusters
    if ~isempty(obj.hFigs) % session already running
        return;
    end

    closeStaleSplitFigures();

    obj.isEnding = 0;
    obj.cRes = struct('hClust', obj.hClust);
    obj.maxAmp = obj.hCfg.maxAmp;

    % select first cluster
    obj.selected = 1;
    obj.currentSite = obj.hClust.clusterSites(1);

    obj.plotAllFigures();
end

function closeStaleSplitFigures()
    try
        delete(findall(groot, 'Type', 'figure', 'Tag', 'FigSplit'));
    catch
    end
    try
        delete(findall(groot, 'Type', 'figure', 'Name', 'Figures'));
    catch
    end
    try
        hFigs = findall(groot, 'Type', 'figure');
        figNames = get(hFigs, 'Name');
        if ischar(figNames)
            figNames = {figNames};
        end
        delete(hFigs(startsWith(figNames, 'Split cluster')));
    catch
    end
end
