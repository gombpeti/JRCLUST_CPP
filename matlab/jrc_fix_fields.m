function jrc_fix_fields(prmFile, varargin)
%JRC_FIX_FIELDS Make every per-cluster field the right length again.
%
%   Repairs the failure behind
%
%       Warning: failed to delete unit 1199
%       Index exceeds the number of array elements. Index must not exceed 1159.
%           in jrclust.interfaces.Clustering/deleteUnit
%
%   when Edit > Delete annotated is used after manual curation.
%
%   JRCLUST keeps a set of per-cluster fields - notes, sites, counts, quality
%   metrics, mean waveforms, the similarity matrix - and deleteUnit subsets
%   ALL of them at once with
%
%       keepSubset = setdiff(goodUnits, unitID)
%       metadata.(fn) = vals(keepSubset)
%
%   `keepSubset` holds cluster IDs, so it indexes up to max(spikeClusters).
%   If any one of those fields is shorter than that, the indexing throws, and
%   because the fields are subsetted in a loop the deletion dies partway
%   through, leaving the clustering half-modified.
%
%   Splitting a cluster in the GUI grows the cluster count, and every
%   per-cluster field has to grow with it. A field that was written into the
%   _res.mat from outside JRCLUST - by the C++ import, or by any script that
%   added one - is not in JRCLUST's own list, so it is never extended, and
%   the lengths drift apart silently until a deletion exposes them.
%
%   This function brings every field back to exactly max(spikeClusters)
%   entries: recomputing the ones derivable from the spike table, padding the
%   rest. Quality metrics for clusters created by splitting cannot be
%   recovered and are set to NaN, which the GUI shows as blank rather than
%   as a wrong number; recompute them with a fresh import if you need them.
%
%   jrc_fix_fields(prmFile)
%   jrc_fix_fields(prmFile, 'outputDir', 'F:\...\resort', 'dryRun', true)
%
%   Options
%     'outputDir'  where the _res.mat lives
%     'dryRun'     report the mismatches without writing (default false)

    opts = struct('outputDir', '', 'dryRun', false);
    for k = 1:2:numel(varargin)
        name = varargin{k};
        if ~isfield(opts, name)
            error('jrc_fix_fields:badOption', 'unknown option ''%s''', name);
        end
        opts.(name) = varargin{k+1};
    end

    resFile = jrc_res_path(prmFile, opts.outputDir);
    fprintf('Loading %s\n', resFile);
    res = load(resFile);

    clu = double(res.spikeClusters(:));
    nClusters = max(clu);
    fprintf('  %d spikes, max(spikeClusters) = %d\n\n', numel(clu), nClusters);

    % JRCLUST's own per-cluster fields, from json/Clustering.json.
    vectorFields = {'clusterNotes','clusterSites','spikesByCluster','unitCount', ...
                    'unitPeaks','unitPeaksRaw','unitPeakSites','unitVpp', ...
                    'unitVppRaw','unitISIRatio','unitIsoDist','unitLRatio'};
    wfFields  = {'meanWfGlobal','meanWfGlobalRaw','meanWfLocal','meanWfLocalRaw', ...
                 'meanWfRawLow','meanWfRawHigh'};

    % ---------------------------------------------------------- report
    fprintf('%-20s %10s %10s   %s\n', 'field', 'found', 'expected', 'action');
    nFixed = 0;

    % Fields that can be rebuilt exactly from the spike table come first, so
    % that padding never invents a count or a site that the data disagrees
    % with.
    counts = accumarray(clu(clu > 0), 1, [nClusters 1]);
    byCluster = arrayfun(@(iC) find(clu == iC), (1:nClusters)', 'UniformOutput', 0);

    [res, nFixed] = setField(res, 'unitCount', counts, nClusters, nFixed, 'rebuilt');
    [res, nFixed] = setField(res, 'spikesByCluster', byCluster, nClusters, nFixed, 'rebuilt');

    if isfield(res, 'spikeSites')
        sites = zeros(nClusters, 1);
        ss = double(res.spikeSites(:));
        for iC = 1:nClusters
            if ~isempty(byCluster{iC}), sites(iC) = mode(ss(byCluster{iC})); end
        end
        [res, nFixed] = setField(res, 'clusterSites', sites, nClusters, nFixed, 'rebuilt');
        [res, nFixed] = setField(res, 'unitPeakSites', sites, nClusters, nFixed, 'rebuilt');
    end

    % The rest are padded or trimmed.
    for i = 1:numel(vectorFields)
        fn = vectorFields{i};
        if ~isfield(res, fn) || isempty(res.(fn)), continue; end
        if any(strcmp(fn, {'unitCount','spikesByCluster','clusterSites','unitPeakSites'}))
            continue;   % already rebuilt above
        end
        v = res.(fn);
        n = numel(v);
        if n == nClusters
            fprintf('%-20s %10d %10d   ok\n', fn, n, nClusters);
            continue;
        end
        res.(fn) = padVector(v, nClusters);
        fprintf('%-20s %10d %10d   %s\n', fn, n, nClusters, ...
                ternary(n < nClusters, 'padded', 'trimmed'));
        nFixed = nFixed + 1;
    end

    % Mean waveforms: clusters are the 3rd dimension.
    for i = 1:numel(wfFields)
        fn = wfFields{i};
        if ~isfield(res, fn) || isempty(res.(fn)), continue; end
        v = res.(fn);
        n = size(v, 3);
        if n == nClusters
            fprintf('%-20s %10d %10d   ok\n', fn, n, nClusters);
            continue;
        end
        if n < nClusters
            filler = zeros(size(v,1), size(v,2), nClusters - n, 'like', v);
            res.(fn) = cat(3, v, filler);
        else
            res.(fn) = v(:,:,1:nClusters);
        end
        fprintf('%-20s %10d %10d   %s\n', fn, n, nClusters, ...
                ternary(n < nClusters, 'padded', 'trimmed'));
        nFixed = nFixed + 1;
    end

    % clusterCentroids is nClusters x 2.
    if isfield(res, 'clusterCentroids') && ~isempty(res.clusterCentroids)
        v = res.clusterCentroids;
        n = size(v, 1);
        if n ~= nClusters
            if n < nClusters
                res.clusterCentroids = [v; nan(nClusters-n, size(v,2))];
            else
                res.clusterCentroids = v(1:nClusters, :);
            end
            fprintf('%-20s %10d %10d   %s\n', 'clusterCentroids', n, nClusters, ...
                    ternary(n < nClusters, 'padded', 'trimmed'));
            nFixed = nFixed + 1;
        else
            fprintf('%-20s %10d %10d   ok\n', 'clusterCentroids', n, nClusters);
        end
    end

    % waveformSim is square, nClusters x nClusters.
    if isfield(res, 'waveformSim') && ~isempty(res.waveformSim)
        v = res.waveformSim;
        n = size(v, 1);
        if n ~= nClusters
            W = nan(nClusters, nClusters, 'like', v);
            m = min(n, nClusters);
            W(1:m, 1:m) = v(1:m, 1:m);
            res.waveformSim = W;
            fprintf('%-20s %10d %10d   %s\n', 'waveformSim', n, nClusters, ...
                    ternary(n < nClusters, 'padded', 'trimmed'));
            nFixed = nFixed + 1;
        else
            fprintf('%-20s %10d %10d   ok\n', 'waveformSim', n, nClusters);
        end
    end

    % Cluster centres must point at a real spike, or the GUI throws
    % "Array indices must be positive integers" the moment it opens.
    if isfield(res, 'clusterCenters')
        res = rebuildCenters(res, byCluster, nClusters);
        fprintf('%-20s %10s %10d   rebuilt\n', 'clusterCenters', '-', nClusters);
    end

    fprintf('\n%d field(s) needed repair\n', nFixed);

    if opts.dryRun
        fprintf('dry run: nothing written\n');
        return;
    end
    if nFixed == 0
        fprintf('nothing to write\n');
        return;
    end

    fprintf('\nWriting %s\n', resFile);
    jrclust.utils.saveStruct(res, resFile);
    fprintf('Done. Reopen with jrc manual and try Delete annotated again.\n');
end

% ---------------------------------------------------------------------------
function [res, nFixed] = setField(res, fn, value, nClusters, nFixed, action)
    if ~isfield(res, fn), return; end
    old = numel(res.(fn));
    if old == nClusters && isequal(size(res.(fn)(:)), size(value(:)))
        % Even at the right length the values may be stale, so rebuild anyway
        % but report it as unchanged in size.
        res.(fn) = reshapeLike(res.(fn), value);
        fprintf('%-20s %10d %10d   ok\n', fn, old, nClusters);
        return;
    end
    res.(fn) = reshapeLike(res.(fn), value);
    fprintf('%-20s %10d %10d   %s\n', fn, old, nClusters, action);
    nFixed = nFixed + 1;
end

function out = reshapeLike(oldVal, newVal)
    %RESHAPELIKE Keep the original orientation, which JRCLUST relies on.
    %   spikesByCluster in particular must stay a column cell array, and
    %   spikesBySite a row one; concatenation in the GUI fails otherwise.
    out = newVal;
    if isrow(oldVal) && ~isrow(out), out = out(:)'; end
    if iscolumn(oldVal) && ~iscolumn(out), out = out(:); end
end

function v = padVector(v, n)
    cur = numel(v);
    if cur == n, return; end
    if cur > n
        v = v(1:n);
        return;
    end
    if iscell(v)
        filler = repmat({''}, n - cur, 1);
        if isrow(v), filler = filler'; end
        v = [v(:); filler(:)];
    else
        filler = nan(n - cur, 1, 'like', double(v));
        v = [double(v(:)); filler];
    end
end

function res = rebuildCenters(res, byCluster, nClusters)
    %REBUILDCENTERS Point each cluster at one of its own surviving spikes.
    %
    %   A centre of 0, or one pointing at a spike that no longer exists, makes
    %   the GUI fail with "Array indices must be positive integers or logical
    %   values" as soon as it loads.
    centers = zeros(nClusters, 1, 'int32');
    useRho = isfield(res, 'spikeRho') && ~isempty(res.spikeRho);
    if useRho, rho = double(res.spikeRho(:)); end
    for iC = 1:nClusters
        idx = byCluster{iC};
        if isempty(idx), centers(iC) = 1; continue; end   % never 0
        if useRho
            [~, best] = max(rho(idx));
            centers(iC) = int32(idx(best));
        else
            centers(iC) = int32(idx(1));
        end
    end
    res.clusterCenters = reshapeLike(res.clusterCenters, centers);
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
