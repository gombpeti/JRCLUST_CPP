function jrc_split_noise(prmFile, varargin)
%JRC_SPLIT_NOISE Split each unit's non-spike-shaped waveforms into a sibling.
%
%   For every unit, scores each INDIVIDUAL spike against a template derived
%   from the recording itself, and moves the ones that do not look like an
%   action potential into a NEW cluster placed immediately after its parent.
%   The new cluster is annotated 'splitted_to_noise'.
%
%   Nothing is deleted. The suspect spikes sit next to the unit they came
%   from, so in jrc manual you can step from a unit to its sibling, compare
%   them directly, and decide: delete the sibling, merge it back, or keep
%   both. That is the point of doing it this way rather than thresholding -
%   the decision stays with the person who can see both halves, and no cut I
%   choose can quietly lose a real cell.
%
%   Cluster ids after each split shift up by one, so a recording of N units
%   with K splits ends with N+K units and the pairs read as
%
%       ... 41  42  42a  43  44  44a ...   (as 41 42 43 44 45 46 after
%                                           renumbering, parents and
%                                           siblings adjacent)
%
%   jrc_split_noise(prmFile, 'dryRun', true)
%   jrc_split_noise(prmFile, 'minScore', 0.4)
%
%   Options
%     'madFactor'    how far below its unit's OWN median a spike must fall
%                    to be moved out, in median absolute deviations
%                    (default 4). This is the main knob. Judging each spike
%                    against its own unit removes the dependence on spike
%                    SNR: an absolute cut selects the noisiest spikes rather
%                    than the wrongly shaped ones, and on a 755-unit
%                    recording took 17-23% out of the very best units at any
%                    threshold from 0.10 to 0.40.
%     'maxAdaptiveCut'  the adaptive cut is never placed above this
%                    (default 0.5), so a spike that matches the template
%                    well is never called an outlier merely because its unit
%                    is unusually tight.
%     'minScore'     an additional absolute floor, 0 to disable (default 0).
%     'minSplit'     do not bother splitting off fewer than this many spikes
%                    (default 20); a handful of stragglers is not worth a
%                    cluster.
%     'minKeep'      leave the unit alone if splitting would take more than
%                    this fraction of it (default 0.9). Above that the unit
%                    is not a neuron with contamination, it is simply not a
%                    neuron, and relabelling it as a whole is the honest
%                    answer - jrc_relabel already does that.
%     'maxSpikes'    cap spikes scored per unit, for speed (default 0, all)
%     'outputDir'    where the _res.mat lives
%     'dryRun'       report only, write nothing (default false)
%
%   Run jrc_relabel afterwards if you want the usual annotations refreshed;
%   the new clusters change every unit id above them.

    opts = struct('madFactor', 4, 'maxAdaptiveCut', 0.5, 'minScore', 0, ...
                  'minSplit', 20, 'minKeep', 0.9, ...
                  'maxSpikes', 0, 'minTemplateSpikes', 100, 'maxSpread', 0.15, ...
                  'outputDir', '', 'dryRun', false);
    for k = 1:2:numel(varargin)
        name = varargin{k};
        if ~isfield(opts, name)
            error('jrc_split_noise:badOption', 'unknown option ''%s''', name);
        end
        opts.(name) = varargin{k+1};
    end

    resFile = jrc_res_path(prmFile, opts.outputDir);
    fprintf('Loading %s\n', resFile);
    res = load(resFile);

    fs = sampleRateOf(prmFile);
    clu = double(res.spikeClusters(:));
    nClusters = max(clu);
    fprintf('  %d spikes, %d clusters\n\n', numel(clu), nClusters);

    % --- template, derived with no labels and checked against physiology --
    % The raw template is wanted only for its verdict - the scoring template
    % is rebuilt from the filtered waveforms below, since those are what the
    % spikes are scored on and they are a different length.
    [~, ok, why] = spikeTemplate(res, nClusters, fs);
    fprintf('template: %s\n', why);
    if ~ok
        error('jrc_split_noise:badTemplate', ...
              ['the template is not spike-shaped, so there is nothing ' ...
               'trustworthy to score against']);
    end

    % --- per-spike scores -------------------------------------------------
    prefix = regexprep(resFile, '_res\.mat$', '');
    % FILTERED waveforms, not raw.
    %
    % A single raw waveform is mostly noise, so its correlation with any
    % template is dominated by that noise rather than by shape. Measured over
    % the same units: on raw, the 5th percentile of a 608,952-spike unit's
    % own spikes was 0.18 - indistinguishable from a noise cluster at 0.14 -
    % and scoring on raw took 17 to 23% of the spikes out of the best units
    % in the recording at any threshold. On filtered, that same unit's 5th
    % percentile is 0.58 while the noise clusters sit at 0.13 to 0.20, which
    % is a real gap.
    %
    % Filtered is also what the clustering itself used, so judging spikes on
    % it asks the same question the sorter asked.
    filtFile = [prefix '_filt.jrc'];
    if ~isfile(filtFile)
        error('jrc_split_noise:noFilt', '%s not found', filtFile);
    end
    shape = double(res.filtShape);
    m = memmapfile(filtFile, 'Format', {'int16', shape, 'x'});
    % The template verified above is the RAW one, and it has to be: the
    % filtered trace has been through a differentiator, so a half-width
    % measured on it would read far narrower than the membrane really moved
    % and the physiology check would be meaningless. But the spikes being
    % scored are filtered and a different length (21 samples against 41), so
    % the SCORING template is rebuilt from the filtered unit means by the
    % same recipe. Raw verifies; filtered scores.
    tmpl = alignedMedian(res.meanWfGlobal, double(res.clusterSites(:)), nClusters);
    if isempty(tmpl)
        error('jrc_split_noise:noFiltTemplate', ...
              'could not build a filtered template from meanWfGlobal');
    end

    % --- refine the template from the tightest units ----------------------
    % A unit whose spikes all look alike is, by that fact alone, a well
    % formed unit - no labels needed to see it. Measured on filtered
    % waveforms, units confirmed as neurons had an interquartile spread of
    % 0.05 to 0.13 across their own spikes while units confirmed as noise sat
    % at 0.12 to 0.30. So: score once against the rough template, keep the
    % units that are both tight and well populated, and rebuild the template
    % from THEIR individual spikes rather than from unit means.
    %
    % Averaging individual spikes rather than unit means matters: a mean
    % waveform has already had its own outliers averaged in, so a template
    % built from means carries a little of every unit's contamination.
    tc0 = tmpl - mean(tmpl); tc0 = tc0 / norm(tc0);
    spread = nan(nClusters,1);
    for c = 1:nClusters
        idx = find(clu == c);
        if numel(idx) < opts.minTemplateSpikes, continue; end
        idx = idx(round(linspace(1, numel(idx), min(300, numel(idx)))));
        r = scoreSpikes(m, idx, tc0);
        if numel(r) >= 20, spread(c) = prctile(r,75) - prctile(r,25); end
    end
    tight = find(isfinite(spread) & spread <= opts.maxSpread);
    if numel(tight) >= 5
        acc = [];
        for c = tight(:)'
            idx = find(clu == c);
            idx = idx(round(linspace(1, numel(idx), min(200, numel(idx)))));
            W = double(squeeze(m.Data.x(:,1,idx)));
            tr = min(W,[],1); ok = tr < 0;
            if ~any(ok), continue; end
            acc = [acc, W(:,ok) ./ abs(tr(ok))]; %#ok<AGROW>
        end
        if size(acc,2) >= 100
            tmpl = median(acc, 2);
            fprintf('template refined from %d spikes of %d tight units (spread <= %.2f)\n', ...
                    size(acc,2), numel(tight), opts.maxSpread);
        end
    else
        fprintf('too few tight units to refine the template; using the rough one\n');
    end
    tc = tmpl - mean(tmpl); tc = tc / norm(tc);

    isNoise = false(numel(clu), 1);
    nSplit = 0; movedTotal = 0; parents = [];
    fprintf('\nscoring spikes');
    for c = 1:nClusters
        if mod(c, 50) == 0, fprintf('.'); end
        idx = find(clu == c);
        if numel(idx) < opts.minSplit, continue; end

        use = idx;
        if opts.maxSpikes > 0 && numel(use) > opts.maxSpikes
            use = use(round(linspace(1, numel(use), opts.maxSpikes)));
        end

        W = double(squeeze(m.Data.x(:, 1, use)));
        tr = min(W, [], 1);
        r = nan(numel(use), 1);
        live = tr < 0;
        if any(live)
            Wn = W(:, live) ./ abs(tr(live));
            Wc = Wn - mean(Wn, 1);
            nrm = sqrt(sum(Wc.^2, 1));
            good = nrm > 0;
            rr = nan(sum(live), 1);
            rr(good) = (Wc(:, good)' * tc) ./ nrm(good)';
            r(live) = rr;
        end
        % Which spikes are outliers IN THIS UNIT.
        %
        % An absolute cut was tried first and does not work on its own: a
        % single waveform is noisy, so its correlation sits well below its
        % unit's mean however good the unit is. Scoring RAW waveforms, a cut
        % at 0.40 moved 26% of all spikes and one at 0.10 still moved 18% -
        % selecting on spike SNR rather than on spike shape, and taking
        % 17-23% out of the best units in the recording. Filtered waveforms
        % (used here) separate far better, but the spread still varies from
        % unit to unit, so the cut is set per unit rather than globally.
        %
        % Comparing each spike with its OWN unit's distribution removes that
        % dependence. A clean unit has a tight spread and loses almost
        % nothing; a contaminated one has a separable tail, which is exactly
        % the population worth splitting off. MAD rather than the standard
        % deviation, so the contamination cannot inflate the very scale used
        % to detect it.
        valid = isfinite(r);
        bad = ~valid;                     % no negative deflection at all
        if sum(valid) >= 20
            med = median(r(valid));
            mad_ = median(abs(r(valid) - med));
            if mad_ > 0
                cut = med - opts.madFactor * mad_;
                % Never cut above the absolute floor: a spike matching the
                % template well is not an outlier, however tight its unit.
                cut = min(cut, opts.maxAdaptiveCut);
                bad = bad | (valid & r < cut);
            end
        end
        if opts.minScore > 0
            bad = bad | (valid & r < opts.minScore);
        end

        if sum(bad) < opts.minSplit, continue; end
        if sum(bad) > opts.minKeep * numel(use), continue; end

        isNoise(use(bad)) = true;
        nSplit = nSplit + 1;
        movedTotal = movedTotal + sum(bad);
        parents(end+1) = c; %#ok<AGROW>
    end
    fprintf(' done\n\n');

    fprintf('units split                 : %d of %d\n', nSplit, nClusters);
    fprintf('spikes moved to siblings    : %d (%.2f%% of all)\n', ...
            movedTotal, 100*movedTotal/numel(clu));
    fprintf('clusters after the split    : %d\n', nClusters + nSplit);

    if nSplit == 0
        fprintf('\nnothing to split\n');
        return;
    end

    fprintf('\nlargest splits:\n');
    fprintf('%8s %10s %10s %8s\n', 'unit', 'kept', 'moved out', '%% moved');
    moved = arrayfun(@(c) sum(isNoise & clu == c), parents);
    kept  = arrayfun(@(c) sum(~isNoise & clu == c), parents);
    [~, o] = sort(moved, 'descend');
    for k = o(1:min(15, numel(o)))
        fprintf('%8d %10d %10d %7.0f%%\n', parents(k), kept(k), moved(k), ...
                100*moved(k)/(moved(k)+kept(k)));
    end

    if opts.dryRun
        fprintf('\ndry run: nothing written\n');
        return;
    end

    % --- renumber: each sibling goes immediately after its parent ---------
    % Parent c keeps a new id of c + (number of splits at or below c) - 1,
    % and its sibling takes the id straight after. Walking the parents in
    % order and counting how many came before gives both without a search.
    isParent = false(nClusters, 1);
    isParent(parents) = true;
    newParentId = zeros(nClusters, 1);
    shift = 0;
    for c = 1:nClusters
        newParentId(c) = c + shift;
        if isParent(c), shift = shift + 1; end
    end

    newClu = clu;
    pos = clu > 0;
    newClu(pos) = newParentId(clu(pos));
    newClu(isNoise) = newClu(isNoise) + 1;        % sibling sits just after
    res.spikeClusters = int32(newClu);

    nNew = nClusters + nSplit;
    res = reindexClusters(res, nClusters, nNew, newParentId, isParent);

    % Recompute everything the GUI reads but never regenerates.
    %
    % Carrying the parent's mean waveform over to its sibling keeps the file
    % structurally valid, but it is not TRUE of either half any more: the
    % parent has lost spikes and the sibling has never had a waveform of its
    % own. Worse, waveformSim is read straight from the file by the GUI and
    % recomputed only on demand, so leaving it blank kills the similarity
    % matrix and the space key stops finding the nearest unit.
    %
    % So hand the new spike table to JRCLUST's own Clustering object and let
    % it recompute mean waveforms, similarity, centroids and quality scores.
    % Using its methods rather than reimplementing them is what guarantees
    % the result has the same structure the GUI expects.
    fprintf('\nWriting %s\n', resFile);
    jrclust.utils.saveStruct(res, resFile);

    % Written first, then handed to jrc_recompute, so the two paths share one
    % implementation: the same call repairs a file that was split before this
    % step existed.
    jrc_recompute(prmFile, 'outputDir', opts.outputDir);
    fprintf(['\nDone. Open the session and step through the pairs: each\n' ...
             '''splitted_to_noise'' cluster sits immediately after the unit\n' ...
             'it came from.\n\n' ...
             '    jrc manual %s\n\n' ...
             'Then Edit > Recompute all mean waveforms - the new clusters\n' ...
             'have no waveform of their own yet.\n'], prmFile);
end

% ---------------------------------------------------------------------------
function res = reindexClusters(res, nOld, nNew, newParentId, isParent)
    %REINDEXCLUSTERS Grow every per-cluster field to the new numbering.
    %
    %   Each field has to be rebuilt rather than merely padded, because the
    %   siblings are INSERTED among the parents rather than appended. Any
    %   field left at the old length, or shuffled against the others, breaks
    %   deleteUnit later with "Index exceeds the number of array elements" -
    %   which is how a half-finished delete corrupts a clustering.
    map = zeros(nOld, 1);
    for c = 1:nOld, map(c) = newParentId(c); end

    clu = double(res.spikeClusters(:));
    counts = accumarray(clu(clu > 0), 1, [nNew 1]);
    res.unitCount = counts;
    res.spikesByCluster = arrayfun(@(iC) find(clu == iC), (1:nNew)', 'UniformOutput', 0);

    if isfield(res, 'spikeSites')
        ss = double(res.spikeSites(:));
        sites = zeros(nNew, 1);
        for iC = 1:nNew
            idx = res.spikesByCluster{iC};
            if ~isempty(idx), sites(iC) = mode(ss(idx)); end
        end
        res.clusterSites = sites;
        if isfield(res, 'unitPeakSites'), res.unitPeakSites = sites; end
    end

    % Notes: carry the parent's across, mark the sibling.
    old = repmat({''}, nOld, 1);
    if isfield(res, 'clusterNotes') && numel(res.clusterNotes) >= nOld
        old = res.clusterNotes(1:nOld);
    end
    notes = repmat({''}, nNew, 1);
    for c = 1:nOld
        notes{map(c)} = old{c};
        if isParent(c), notes{map(c)+1} = 'splitted_to_noise'; end
    end
    res.clusterNotes = notes(:);

    % Quality scores no longer describe either half. NaN reads as unknown in
    % the GUI; a stale number would read as fact.
    for fn = {'unitPeaks','unitPeaksRaw','unitVpp','unitVppRaw', ...
              'unitISIRatio','unitIsoDist','unitLRatio','unitSNR', ...
              'nSitesOverThresh'}
        f = fn{1};
        if isfield(res, f) && numel(res.(f)) >= nOld
            res.(f) = nan(nNew, 1);
        end
    end

    % Mean waveforms: give the sibling a copy of its parent's for now, so the
    % array is the right size and the GUI can open. Recompute all mean
    % waveforms replaces both.
    for fn = {'meanWfGlobal','meanWfGlobalRaw','meanWfLocal','meanWfLocalRaw', ...
              'meanWfRawLow','meanWfRawHigh'}
        f = fn{1};
        if ~isfield(res, f) || ndims(res.(f)) ~= 3 || size(res.(f),3) < nOld
            continue;
        end
        v = res.(f);
        out = zeros(size(v,1), size(v,2), nNew, 'like', v);
        for c = 1:nOld
            out(:,:,map(c)) = v(:,:,c);
            if isParent(c), out(:,:,map(c)+1) = v(:,:,c); end
        end
        res.(f) = out;
    end

    if isfield(res, 'clusterCentroids') && size(res.clusterCentroids,1) >= nOld
        v = res.clusterCentroids;
        out = nan(nNew, size(v,2));
        for c = 1:nOld
            out(map(c),:) = v(c,:);
            if isParent(c), out(map(c)+1,:) = v(c,:); end
        end
        res.clusterCentroids = out;
    end

    % Waveform similarity is nClusters x nClusters; rebuilding it properly
    % needs the new waveforms, so blank it and let the GUI recompute.
    if isfield(res, 'waveformSim')
        res.waveformSim = nan(nNew, nNew, 'like', res.waveformSim);
    end

    % Centres must point at a spike the cluster actually owns, or the GUI
    % fails with "Array indices must be positive integers" on load.
    if isfield(res, 'clusterCenters')
        centers = zeros(nNew, 1, 'int32');
        useRho = isfield(res, 'spikeRho') && ~isempty(res.spikeRho);
        if useRho, rho = double(res.spikeRho(:)); end
        for iC = 1:nNew
            idx = res.spikesByCluster{iC};
            if isempty(idx), centers(iC) = 1; continue; end
            if useRho
                [~, b] = max(rho(idx));
                centers(iC) = int32(idx(b));
            else
                centers(iC) = int32(idx(1));
            end
        end
        res.clusterCenters = centers;
    end
end

function tmpl = alignedMedian(wf, sites, nClusters)
    %ALIGNEDMEDIAN Median shape over the units already aligned on their trough.
    %   Each unit is normalised to its own trough depth, so amplitude, gain
    %   and probe divide out and only shape remains.
    tmpl = [];
    if isempty(wf), return; end
    nS = size(wf,1);
    raw = nan(nS, nClusters); at = nan(nClusters,1);
    for c = 1:min(nClusters, size(wf,3))
        pk = sites(c);
        if pk < 1 || pk > size(wf,2), continue; end
        w = double(wf(:,pk,c));
        [tr, it] = min(w);
        if tr >= 0, continue; end
        raw(:,c) = w / abs(tr);
        at(c) = it;
    end
    if all(isnan(at)), return; end
    centre = mode(at(~isnan(at)));
    sel = at == centre;
    if sum(sel) < 20, return; end
    tmpl = median(raw(:, sel), 2);
end

function r = scoreSpikes(m, idx, tc)
    %SCORESPIKES Correlation of each spike's shape with the template.
    %   Each spike is normalised to its own trough first, so the comparison
    %   is of shape rather than of size.
    W = double(squeeze(m.Data.x(:, 1, idx)));
    tr = min(W, [], 1);
    ok = tr < 0;
    r = [];
    if ~any(ok), return; end
    Wn = W(:, ok) ./ abs(tr(ok));
    Wc = Wn - mean(Wn, 1);
    nrm = sqrt(sum(Wc.^2, 1));
    good = nrm > 0;
    r = (Wc(:, good)' * tc) ./ nrm(good)';
end

function [tmpl, ok, why] = spikeTemplate(res, nClusters, fs)
    %SPIKETEMPLATE The recording's own spike shape, verified by physiology.
    %   See jrc_relabel/templateScores for the reasoning; this is the same
    %   derivation, kept here so the splitter stands alone.
    tmpl = []; ok = false; why = 'no waveforms available';
    if ~isfield(res,'meanWfGlobalRaw') || ~isfield(res,'clusterSites'), return; end
    wf = res.meanWfGlobalRaw;
    nS = size(wf,1);
    sites = double(res.clusterSites(:));

    raw = nan(nS, nClusters); troughAt = nan(nClusters,1);
    for c = 1:min(nClusters, size(wf,3))
        pk = sites(c);
        if pk < 1 || pk > size(wf,2), continue; end
        w = double(wf(:,pk,c));
        [tr, it] = min(w);
        if tr >= 0, continue; end
        raw(:,c) = w / abs(tr);
        troughAt(c) = it;
    end
    if all(isnan(troughAt)), return; end

    centre = mode(troughAt(~isnan(troughAt)));
    aligned = troughAt == centre;
    if sum(aligned) < 20, why = 'too few aligned units'; return; end
    tmpl = median(raw(:, aligned), 2);

    [tr, it] = min(tmpl);
    reb = max(tmpl(it:end));
    half = tr/2;
    L=it; while L>1 && tmpl(L-1)<=half, L=L-1; end
    R=it; while R<nS && tmpl(R+1)<=half, R=R+1; end
    halfW = (R-L+1)*1000/fs;
    big = abs(tmpl) > 0.25;
    nExc = sum(diff([false; big(:)]) == 1);

    ok = nExc <= 2 && abs(tr) > abs(reb) && halfW >= 0.1 && halfW <= 1.5 && ...
         reb > 0.05 && abs(median(tmpl(min(nS,it+14):end))) < 0.35;
    if ok
        why = sprintf('spike-shaped (half-width %.3f ms, rebound %.2f, from %d units)', ...
                      halfW, reb, sum(aligned));
    else
        why = sprintf('NOT spike-shaped (%d excursions, half-width %.3f ms)', nExc, halfW);
    end
end

function fs = sampleRateOf(prmFile)
    fs = [];
    txt = fileread(prmFile);
    tok = regexp(txt, '(?m)^\s*sampleRate\s*=\s*([0-9.eE+-]+)', 'tokens', 'once');
    if ~isempty(tok), fs = str2double(tok{1}); end
    if isempty(fs) || ~isfinite(fs) || fs <= 0
        error('jrc_split_noise:noSampleRate', 'could not read sampleRate from %s', prmFile);
    end
end
