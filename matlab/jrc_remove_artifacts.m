function jrc_remove_artifacts(prmFile, varargin)
%JRC_REMOVE_ARTIFACTS Strip spikes from artifact windows, without re-sorting.
%
%   jrc_remove_artifacts('mysession.prm', 'dryRun', true)   % look first
%   jrc_remove_artifacts('mysession.prm')                   % apply
%
%   A real population of neurons fires at a roughly stable aggregate rate.
%   An electrical artifact - an amplifier being stopped, a connector touched,
%   the animal hitting something - crosses threshold on many channels at once
%   and produces a burst of detections orders of magnitude above that rate,
%   confined to a short window.
%
%   Those detections are not spikes, but density-peak clustering has no way to
%   reject them: every spike is assigned by following a chain to its nearest
%   denser neighbour, and that chain always terminates in some cluster. So the
%   artifact gets adopted, and no amount of merging or re-clustering removes
%   it. The only fix is to take the spikes out.
%
%   This finds those windows from the aggregate rate alone and sets the spikes
%   in them to cluster 0. Run jrc_prune afterwards to delete them for good.
%
%   Options:
%     'binSeconds'  bin width for the rate estimate (default 1)
%     'madThresh'   flag bins above median + this many MADs (default 10)
%     'padSeconds'  widen each flagged window by this much (default 0.5),
%                   since an artifact rarely aligns with a bin boundary
%     'windows'     explicit [start end; ...] in seconds, instead of detecting
%     'minRateRatio' never flag a bin below this multiple of the median
%                   (default 3), so a quiet recording cannot flag itself
%     'outputDir'   where the _res.mat lives, if the .prm's outputDir is stale
%     'dryRun'      report only, write nothing (default false)
%
%   This is deliberately blunt: it removes EVERYTHING in the flagged windows,
%   including any genuine spikes. At a fraction of a percent of recording time
%   that is a good trade, but it is a real loss and the report states it.

    opts = struct('binSeconds', 1, 'madThresh', 10, 'padSeconds', 0.5, ...
                  'windows', [], 'minRateRatio', 3, 'outputDir', '', ...
                  'dryRun', false);
    for k = 1:2:numel(varargin)
        name = varargin{k};
        if ~isfield(opts, name)
            error('jrc_remove_artifacts:badOption', 'unknown option ''%s''', name);
        end
        opts.(name) = varargin{k+1};
    end

    resFile = jrc_res_path(prmFile, opts.outputDir);
    fprintf('Loading %s\n', resFile);
    res = load(resFile);

    hCfg = jrclust.Config(prmFile);
    sampleRate = hCfg.sampleRate;

    t = double(res.spikeTimes(:)) / sampleRate;
    clu = double(res.spikeClusters(:));
    nSpikes = numel(t);
    duration = max(t);
    nClustersBefore = max(clu);

    fprintf('  %d spikes, %d clusters, %.0f s (%.2f h)\n', ...
            nSpikes, nClustersBefore, duration, duration/3600);

    % ------------------------------------------------- find artifact windows
    if isempty(opts.windows)
        edges = 0:opts.binSeconds:(duration + opts.binSeconds);
        counts = histcounts(t, edges);

        med = median(counts);
        mad = median(abs(counts - med)) * 1.4826;
        if mad <= 0, mad = max(sqrt(max(med,1)), 1); end

        cutoff = max(med + opts.madThresh*mad, opts.minRateRatio*med);

        fprintf('\n  population rate per %g s bin:\n', opts.binSeconds);
        fprintf('    median %.0f   MAD %.0f   max %.0f  (%.1fx median)\n', ...
                med, mad, max(counts), max(counts)/max(med,1));
        fprintf('    cutoff %.0f  (median + %g MAD, floor %gx median)\n', ...
                cutoff, opts.madThresh, opts.minRateRatio);

        badBins = find(counts > cutoff);
        if isempty(badBins)
            fprintf('\n  no bin exceeds the cutoff; nothing to remove\n');
            return;
        end

        % Merge adjacent bins into windows, then pad.
        starts = edges(badBins);
        stops  = edges(badBins) + opts.binSeconds;
        windows = mergeWindows([starts(:) stops(:)], opts.binSeconds);
        windows(:,1) = max(windows(:,1) - opts.padSeconds, 0);
        windows(:,2) = windows(:,2) + opts.padSeconds;
    else
        windows = opts.windows;
        fprintf('\n  using %d explicitly supplied window(s)\n', size(windows,1));
    end

    % ------------------------------------------------------------- report
    drop = false(nSpikes, 1);
    fprintf('\n  artifact windows:\n');
    for i = 1:size(windows,1)
        inWin = t >= windows(i,1) & t <= windows(i,2);
        drop = drop | inWin;
        fprintf('    %8.1f - %8.1f s  (%5.1f s)  %8d spikes\n', ...
                windows(i,1), windows(i,2), diff(windows(i,:)), sum(inWin));
    end

    droppedTime = sum(diff(windows, 1, 2));
    fprintf('\n  removing %d spikes (%.2f%% of all)\n', sum(drop), ...
            100*sum(drop)/nSpikes);
    fprintf('  discarding %.1f s of recording (%.3f%% of %.0f s)\n', ...
            droppedTime, 100*droppedTime/duration, duration);

    % Which units does this gut? A unit that loses nearly everything was
    % mostly artifact to begin with and will disappear when pruned.
    affected = clu > 0 & drop;
    lost = accumarray(clu(affected), 1, [nClustersBefore 1]);
    total = accumarray(clu(clu>0), 1, [nClustersBefore 1]);
    frac = zeros(nClustersBefore,1);
    nz = total > 0;
    frac(nz) = lost(nz) ./ total(nz);

    fprintf('\n  units affected          : %d / %d\n', sum(lost>0), nClustersBefore);
    fprintf('  units losing >90%% spikes : %d  (these were essentially all artifact)\n', ...
            sum(frac > 0.9));
    fprintf('  units losing >50%%        : %d\n', sum(frac > 0.5));
    fprintf('  units losing >20%%        : %d\n', sum(frac > 0.2));

    worst = find(frac > 0.9);
    if ~isempty(worst)
        [~, ord] = sort(total(worst), 'descend');
        show = worst(ord(1:min(10,numel(ord))));
        fprintf('\n  largest units that were pure artifact:\n');
        fprintf('    unit   spikes   lost\n');
        for i = 1:numel(show)
            fprintf('    %5d  %7d  %5.1f%%\n', show(i), total(show(i)), ...
                    100*frac(show(i)));
        end
    end

    if opts.dryRun
        fprintf('\ndry run: nothing written\n');
        return;
    end

    % --------------------------------------------------------------- apply
    clu(drop) = 0;
    res.spikeClusters = int32(clu);

    % Keep the by-cluster table consistent with what we just changed.
    nClusters = max(clu);
    if isfield(res, 'spikesByCluster')
        res.spikesByCluster = arrayfun(@(iC) find(clu == iC), 1:nClusters, ...
                                       'UniformOutput', 0)';
    end
    if isfield(res, 'unitCount')
        res.unitCount = cellfun(@numel, res.spikesByCluster);
    end

    fprintf('\nWriting %s\n', resFile);
    jrclust.utils.saveStruct(res, resFile);

    fprintf('\nDone. The spikes are now at cluster 0.\n');
    fprintf('To delete them for good:\n    jrc_prune(''%s'')\n', prmFile);
end

%% ------------------------------------------------------------------ helpers

function merged = mergeWindows(w, gap)
    %MERGEWINDOWS Join windows separated by less than `gap`.
    w = sortrows(w);
    merged = w(1,:);
    for i = 2:size(w,1)
        if w(i,1) - merged(end,2) <= gap
            merged(end,2) = max(merged(end,2), w(i,2));
        else
            merged(end+1,:) = w(i,:); %#ok<AGROW>
        end
    end
end
