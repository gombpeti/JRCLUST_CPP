function deleteClusters(obj, deleteMe)
    %DELETECLUSTERS Delete clusters either specified or selected
    if obj.isWorking
        jrclust.utils.qMsgBox('An operation is in progress.');
        return;
    end

    if nargin < 2 && numel(obj.selected) > 1
        return;
    elseif nargin < 2
        deleteMe = obj.selected(1);
    end

    obj.isWorking = 1;
    
    % speculatively delete clusters
    %
    % Only the LATEST snapshot is ever needed: as the input to the next
    % iteration, and as the argument to commit() below. The original code
    % appended every snapshot to a growing struct array and so held one
    % complete copy of the clustering PER DELETED UNIT, all at the same time.
    %
    % Each snapshot carries the whole spike table plus every per-cluster
    % field - spikesByCluster, four mean-waveform arrays, the nClusters x
    % nClusters similarity matrix. Measured on a 512-site, 1143-unit
    % recording that is 514 MB each, so Edit > Delete annotated on 120 units
    % asked for 60 GB and MATLAB stalled partway through. The cost grew with
    % the number of units being removed, which is why this worked on smaller
    % probes and stopped working here rather than breaking outright.
    %
    % Keeping a single snapshot makes the memory cost independent of how many
    % units are deleted. Behaviour is otherwise unchanged: the same snapshots
    % are produced in the same order, only the earlier ones are released.
    res = [];
    deleteMe = sort(deleteMe, 'desc');
    showSubset = obj.showSubset;

    for iCluster = 1:numel(deleteMe) % go backwards to avoid deleting the wrong units after a reorder
        if isempty(res)
            args = {obj.hClust.spikeClusters, deleteMe(iCluster), struct()};
        else
            args = {res.spikeClusters, deleteMe(iCluster), res.metadata};
        end

        res_ = obj.hClust.deleteUnit(args{:});

        showSubset(showSubset == deleteMe(iCluster)) = [];
        mask = showSubset > deleteMe(iCluster);
        showSubset(mask) = showSubset(mask) - 1;

        % operation found to be inconsistent
        if isempty(res_.metadata)
            warning('failed to delete unit %d', deleteMe(iCluster));
            continue;
        end

        res = res_;
    end

    if ~isempty(res)
        msg = ['delete ' strjoin(arrayfun(@num2str, deleteMe, 'UniformOutput', 0), ',')];
        try
            obj.hClust.commit(res.spikeClusters, res.metadata, msg);

            obj.showSubset = showSubset;
        catch ME
            warning('Failed to delete: %s', ME.message);
            jrclust.utils.qMsgBox('Operation failed.');
        end
        
        obj.isWorking = 0; % in case updateSelect needs to zoom

        obj.selected = min([max(obj.showSubset), obj.selected]); % fix OOB error when deleting last cluster

        % replot
        obj.updateFigWav();
        obj.updateFigRD(); % centers changed, need replotting
        obj.updateFigSim();
        obj.updateSelect(obj.selected, 1);
    end

    obj.isWorking = 0;
end