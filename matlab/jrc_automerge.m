function jrc_automerge(prmFile, varargin)
%JRC_AUTOMERGE Merge over-split units, with the cross-correlogram as referee.
%
%   JRCLUST's own Edit > Merge auto thresholds on mean-waveform correlation
%   alone. That cannot distinguish one neuron split in two from two neurons
%   sitting at the same distance from the same sites - their mean waveforms
%   are equally alike in both cases - so the threshold has to be set
%   conservatively, and over-split units survive.
%
%   The spike trains settle it. For a candidate pair, count spike pairs at
%   short lag against pairs at 10-25 ms, scaled by the window widths:
%
%     ratio near 0   the two units AVOID each other at short lag. Only one
%                    refractory period can do that, so they are one cell that
%                    has been split - MERGE.
%
%     ratio near 1   the units are independent, which is what two different
%                    neurons look like however similar their waveforms -
%                    LEAVE THEM ALONE.
%
%     ratio far      they fire at the same instant: the same event detected
%     above 1        twice rather than two cells - a duplicate, which merging
%                    would bury rather than fix.
%
%   So waveform similarity proposes and the cross-correlogram disposes. Only
%   pairs that pass both are merged.
%
%   jrc_automerge(prmFile, 'dryRun', true)
%   jrc_automerge(prmFile, 'maxUnitSim', 0.95)
%
%   Options
%     'maxUnitSim'     waveform correlation above which a pair is a CANDIDATE
%                      (default 0.95). Measured on a 918-unit recording, a
%                      threshold of 0.98 - JRCLUST's default - proposed zero
%                      pairs, which is why auto-merge appeared to do nothing.
%     'maxCrossRatio'  cross-correlogram ratio below which a candidate pair is
%                      accepted as one cell (default 0.5). A clean refractory
%                      dip gives well under this; independent cells give ~1.
%     'minSpikes'      skip units smaller than this; their correlograms are
%                      too sparse to judge (default 50)
%     'outputDir'      where the _res.mat lives
%     'dryRun'         report only, write nothing (default false)
%
%   Merging makes the stored quality metrics and mean waveforms stale. Those
%   are set to NaN rather than left wrong; recompute them in the GUI with
%   Edit > Recompute all mean waveforms.

    opts = struct('maxUnitSim', 0.95, 'maxCrossRatio', 0.5, 'minSpikes', 50, ...
                  'outputDir', '', 'dryRun', false);
    for k = 1:2:numel(varargin)
        name = varargin{k};
        if ~isfield(opts, name)
            error('jrc_automerge:badOption', 'unknown option ''%s''', name);
        end
        opts.(name) = varargin{k+1};
    end

    resFile = jrc_res_path(prmFile, opts.outputDir);
    fprintf('Loading %s\n', resFile);
    res = load(resFile);

    if ~isfield(res, 'waveformSim') || isempty(res.waveformSim)
        error('jrc_automerge:noSim', ...
              'this _res.mat has no waveformSim; re-import with jrc_import_cpp');
    end

    fs  = sampleRateOf(prmFile);
    t   = double(res.spikeTimes(:)) / fs;
    clu = double(res.spikeClusters(:));
    nClusters = max(clu);
    counts = accumarray(clu(clu > 0), 1, [nClusters 1]);
    fprintf('  %d spikes, %d clusters\n\n', numel(t), nClusters);

    % ------------------------------------------------------- candidates
    W = double(res.waveformSim);
    W = W(1:nClusters, 1:nClusters);
    W(1:nClusters+1:end) = -Inf;            % never propose a unit against itself
    W(~isfinite(W)) = -Inf;
    [ia, ib] = find(triu(W >= opts.maxUnitSim, 1));

    fprintf('candidate pairs at waveformSim >= %.2f : %d\n', opts.maxUnitSim, numel(ia));
    if isempty(ia)
        fprintf(['\nNothing proposed. The threshold is above every pair in ' ...
                 'this\nrecording - lower it and run again.\n']);
        showSpectrum(W, nClusters);
        return;
    end

    % Group the spikes once; masking the full vector per pair would dominate.
    byCluster = cell(nClusters, 1);
    keep = clu > 0;
    [sc, ix] = sort(clu(keep));
    tt = t(keep); tt = tt(ix);
    ed = [0; find(diff(sc)); numel(sc)];
    for k = 1:numel(ed)-1
        byCluster{sc(ed(k)+1)} = tt(ed(k)+1 : ed(k+1));
    end

    % ------------------------------------------------------- referee
    fprintf('\n%12s %12s %12s  %s\n', 'pair', 'waveformSim', 'crossRatio', 'verdict');
    accept = false(numel(ia), 1);
    ratios = nan(numel(ia), 1);
    for k = 1:numel(ia)
        a = ia(k); b = ib(k);
        if counts(a) < opts.minSpikes || counts(b) < opts.minSpikes
            verdict = 'too few spikes to judge';
        else
            ratios(k) = crossRatio(byCluster{a}, byCluster{b});
            if ~isfinite(ratios(k))
                verdict = 'no overlap to judge';
            elseif ratios(k) < opts.maxCrossRatio
                verdict = 'ONE CELL - merge';
                accept(k) = true;
            elseif ratios(k) > 5
                verdict = 'duplicate detections - not a merge';
            else
                verdict = 'independent cells - keep separate';
            end
        end
        if k <= 40 || accept(k)
            fprintf('%6d-%-5d %12.4f %12.3f  %s\n', a, b, W(a,b), ratios(k), verdict);
        end
    end
    if numel(ia) > 40
        fprintf('  ... %d further candidates not listed\n', numel(ia) - 40);
    end

    fprintf('\n%-38s %6d\n', 'candidates', numel(ia));
    fprintf('%-38s %6d\n', 'accepted as one cell', sum(accept));
    fprintf('%-38s %6d\n', 'rejected (independent or duplicate)', sum(~accept));

    if ~any(accept)
        fprintf('\nnothing to merge\n');
        return;
    end

    % ------------------------------------------------------- merge
    % Union-find, so that a-b and b-c collapse into one unit rather than
    % merging twice and leaving an orphan.
    parent = (1:nClusters)';
    for k = find(accept)'
        ra = findRoot(parent, ia(k));
        rb = findRoot(parent, ib(k));
        if ra ~= rb, parent(max(ra,rb)) = min(ra,rb); end
    end
    for i = 1:nClusters, parent(i) = findRoot(parent, i); end

    newClu = clu;
    pos = clu > 0;
    newClu(pos) = parent(clu(pos));

    survivors = unique(newClu(newClu > 0));
    fprintf('\n%d units -> %d units (%d merged away)\n', ...
            nClusters, numel(survivors), nClusters - numel(survivors));

    if opts.dryRun
        fprintf('\ndry run: nothing written\n');
        return;
    end

    res.spikeClusters = int32(newClu);
    res = compactAfterMerge(res, nClusters);

    fprintf('\nWriting %s\n', resFile);
    jrclust.utils.saveStruct(res, resFile);

    % Merging makes the mean waveforms, similarity matrix, centroids and
    % quality scores describe a clustering that no longer exists, and the GUI
    % reads all of them from the file rather than regenerating them on load.
    % waveformSim matters most: it is what the space key uses to jump to the
    % most similar unit, so a stale one silently breaks navigation and shows
    % a similarity matrix of the wrong clustering. Telling the user to run
    % Edit > Recompute all mean waveforms was not enough - that rebuilds the
    % waveforms but is easy to forget, and nothing warns when it is skipped.
    jrc_recompute(prmFile, 'outputDir', opts.outputDir);
end

% ---------------------------------------------------------------------------
function showSpectrum(W, ~)
    %SHOWSPECTRUM How many pairs each similarity threshold would propose.
    %   Printed when nothing is proposed at the chosen cut, which is the
    %   useful moment: maxUnitSim 0.98 proposed zero pairs on a real
    %   recording, so Merge auto appeared to do nothing at all.
    U = triu(W, 1);
    U = U(isfinite(U));
    fprintf('\n%10s %12s\n', 'threshold', 'pairs');
    for th = [0.98 0.97 0.96 0.95 0.93 0.90 0.85]
        fprintf('%10.2f %12d\n', th, sum(U >= th));
    end
end

function r = findRoot(parent, i)
    while parent(i) ~= i, i = parent(i); end
    r = i;
end

function res = compactAfterMerge(res, nOld)
    %COMPACTAFTERMERGE Renumber survivors to 1..N and subset per-unit fields.
    good = unique(res.spikeClusters(res.spikeClusters > 0));
    if isempty(good), return; end
    keepUnits = false(nOld, 1);
    keepUnits(good) = true;
    oldToNew = zeros(nOld, 1, 'int32');
    oldToNew(good) = int32(1:numel(good));

    pos = res.spikeClusters > 0;
    res.spikeClusters(pos) = oldToNew(res.spikeClusters(pos));
    n = numel(good);

    % Rebuilt exactly from the spike table, so they cannot disagree with it.
    clu = double(res.spikeClusters(:));
    res.unitCount = accumarray(clu(clu > 0), 1, [n 1]);
    res.spikesByCluster = arrayfun(@(iC) find(clu == iC), (1:n)', 'UniformOutput', 0);
    if isfield(res, 'spikeSites')
        ss = double(res.spikeSites(:));
        sites = zeros(n, 1);
        for iC = 1:n
            idx = res.spikesByCluster{iC};
            if ~isempty(idx), sites(iC) = mode(ss(idx)); end
        end
        res.clusterSites = sites;
        if isfield(res, 'unitPeakSites'), res.unitPeakSites = sites; end
    end

    % Subset what survives unchanged.
    for fn = {'clusterNotes','clusterCenters','nSitesOverThresh'}
        f = fn{1};
        if isfield(res, f) && numel(res.(f)) >= nOld
            res.(f) = res.(f)(keepUnits);
        end
    end

    % A merged unit's shape and quality no longer describe it. NaN says
    % "unknown", which the GUI leaves blank; a stale number would read as
    % fact.
    for fn = {'unitPeaks','unitPeaksRaw','unitVpp','unitVppRaw', ...
              'unitISIRatio','unitIsoDist','unitLRatio','unitSNR'}
        f = fn{1};
        if isfield(res, f) && numel(res.(f)) >= nOld
            res.(f) = nan(n, 1);
        end
    end
    for fn = {'meanWfGlobal','meanWfGlobalRaw','meanWfLocal','meanWfLocalRaw', ...
              'meanWfRawLow','meanWfRawHigh'}
        f = fn{1};
        if isfield(res, f) && ndims(res.(f)) == 3 && size(res.(f),3) >= nOld
            res.(f) = res.(f)(:,:,keepUnits);
        end
    end
    if isfield(res, 'clusterCentroids') && size(res.clusterCentroids,1) >= nOld
        res.clusterCentroids = res.clusterCentroids(keepUnits, :);
    end
    if isfield(res, 'waveformSim') && all(size(res.waveformSim) >= [nOld nOld])
        res.waveformSim = res.waveformSim(keepUnits, keepUnits);
    end

    % Centres must point at a surviving spike or the GUI refuses to open.
    if isfield(res, 'clusterCenters')
        centers = zeros(n, 1, 'int32');
        useRho = isfield(res, 'spikeRho') && ~isempty(res.spikeRho);
        if useRho, rho = double(res.spikeRho(:)); end
        for iC = 1:n
            idx = res.spikesByCluster{iC};
            if isempty(idx), centers(iC) = 1; continue; end
            if useRho
                [~, best] = max(rho(idx));
                centers(iC) = int32(idx(best));
            else
                centers(iC) = int32(idx(1));
            end
        end
        res.clusterCenters = centers;
    end
end

function r = crossRatio(ta, tb)
    r = NaN;
    ta = sort(ta(:)); tb = sort(tb(:));
    if numel(ta) < 20 || numel(tb) < 20, return; end
    near = countWithin(ta, tb, 0.001);
    far  = countWithin(ta, tb, 0.025) - countWithin(ta, tb, 0.010);
    if far <= 0, return; end
    r = (near / far) * 15;   % 30 ms of far window vs 2 ms of near
end

function c = countWithin(ta, tb, win)
    edges = [-inf; tb; inf];
    [~, ~, hi] = histcounts(ta + win, edges);
    [~, ~, lo] = histcounts(ta - win, edges);
    c = sum(hi - lo);
end

function fs = sampleRateOf(prmFile)
    fs = [];
    txt = fileread(prmFile);
    tok = regexp(txt, '(?m)^\s*sampleRate\s*=\s*([0-9.eE+-]+)', 'tokens', 'once');
    if ~isempty(tok), fs = str2double(tok{1}); end
    if isempty(fs) || ~isfinite(fs) || fs <= 0
        error('jrc_automerge:noSampleRate', 'could not read sampleRate from %s', prmFile);
    end
end
