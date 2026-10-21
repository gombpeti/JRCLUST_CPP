function mergeSelected(obj)
    %MERGESELECTED Merge a pair of clusters
    if obj.isWorking
        jrclust.utils.qMsgBox('An operation is in progress.');
        return;
    end

    if numel(obj.selected) < 2
        return;
    end

    % Merge EVERY selected unit, not just the first and last.
    %
    % This used to take min(selected) and max(selected), which was right
    % while the selection was capped at two. Now that any number can be
    % selected it would silently merge the outermost pair and leave the
    % middle ones behind - the worst kind of failure, since the GUI would
    % report a successful merge and the remaining units look merged until
    % you count them.
    %
    % Merged into the lowest id, highest first, so that the renumbering
    % after each step cannot invalidate the ids still to be merged.
    targets = sort(unique(obj.selected(:)'), 'descend');
    iCluster = targets(end);          % the survivor
    toMerge = targets(1:end-1);

    obj.isWorking = 1;

    showSubset = obj.showSubset;
    spikeClusters = obj.hClust.spikeClusters;
    metadata = [];
    merged = [];

    for jCluster = toMerge
        res = obj.hClust.mergeUnits(spikeClusters, iCluster, jCluster);
        if isempty(res.metadata)
            warning('failed to merge units %d and %d', iCluster, jCluster);
            continue;
        end
        spikeClusters = res.spikeClusters;
        metadata = res.metadata;
        merged(end+1) = jCluster; %#ok<AGROW>

        showSubset(showSubset == jCluster) = [];
        mask = showSubset > jCluster;
        showSubset(mask) = showSubset(mask) - 1;
    end

    if isempty(metadata)
        obj.isWorking = 0;
        return;
    end

    msg = sprintf('merge %s into %d', ...
                  strjoin(arrayfun(@num2str, sort(merged), 'UniformOutput', 0), ','), ...
                  iCluster);
    try
        obj.hClust.commit(spikeClusters, metadata, msg);
        obj.showSubset = showSubset;
    catch ME
        warning('Failed to merge: %s', ME.message);
        jrclust.utils.qMsgBox('Operation failed.');
    end

    obj.isWorking = 0; % in case updateSelect needs to zoom

    obj.selected = iCluster; % fix OOB error when merging last cluster

    % replot
    obj.replot();
    obj.updateSelect(iCluster);
end