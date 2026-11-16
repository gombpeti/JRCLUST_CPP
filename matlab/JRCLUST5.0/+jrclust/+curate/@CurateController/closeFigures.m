function closeFigures(obj)
    %CLOSEFIGURES Close all open figures
    %   Deletes the windows directly instead of sending them close requests,
    %   so no close handler (and no question) can run while closing.
    figKeys = keys(obj.hFigs);
    for i = 1:numel(figKeys)
        try
            hFig = obj.hFigs(figKeys{i});
            if hFig.isReady
                hFig.figApply(@set, 'CloseRequestFcn', 'closereq');
                hFig.figApply(@delete);
            end
        catch ME
            warning('Could not close %s: %s', figKeys{i}, ME.message);
        end
    end
    obj.hFigs = containers.Map();

    % leftovers of an interrupted split
    try
        delete(findall(groot, 'Type', 'figure', 'Tag', 'FigSplit'));
        hFigs = findall(groot, 'Type', 'figure');
        figNames = get(hFigs, 'Name');
        if ischar(figNames)
            figNames = {figNames};
        end
        delete(hFigs(startsWith(figNames, 'Split cluster')));
    catch
    end
end
