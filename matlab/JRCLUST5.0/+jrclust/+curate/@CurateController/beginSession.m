function beginSession(obj)
    %BEGINSESSION Start curating clusters
    if ~isempty(obj.hFigs) % session already running
        return;
    end

    closeStaleSplitFigures();

    obj.isEnding = 0;
    obj.isBusy = 0;
    obj.cRes = struct('hClust', obj.hClust);
    fprintf(['JRCLUST: File > Save saves right away; File > Save and exit (Ctrl+Q) or closing the\n' ...
             '         waveform window saves and closes. If the windows ever stop responding:\n' ...
             '         press Ctrl+C here, then type  jrc_save_curation  (or  jrc_save_curation close).\n']);
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
