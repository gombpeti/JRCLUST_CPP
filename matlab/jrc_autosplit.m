function out = jrc_autosplit(prmFile, varargin)
%JRC_AUTOSPLIT Over-split suspect units, drop the rubbish, merge back the rest.
%
%   Automates the curation step that cannot be done by thresholding a whole
%   unit: a cluster that fires below the refractory period is not necessarily
%   worthless. Two independent neurons produce intervals of half a millisecond
%   between them without either cell violating its own refractory period, so
%   the violation is evidence that the cluster holds more than one source, not
%   evidence that the source is noise. Deleting it loses two real neurons.
%
%   The procedure mirrors what curation does by hand:
%
%     1. over-split the unit into maxK groups in its own feature space
%     2. decide, for each PAIR of groups, what their relationship is
%     3. merge the pairs that are one cell, drop duplicates and noise
%     4. keep whatever survives as separate units
%
%   Step 2 is the part that makes this decidable, and it uses the
%   cross-correlogram between the two groups. Counting pairs of spikes at
%   short lag against pairs at 10-25 ms, scaled by the window widths:
%
%     ratio near 0    the two groups AVOID each other at short lag, so they
%                     share one refractory period: one neuron, split in two by
%                     drift or by the clustering being too eager -> MERGE
%
%     ratio near 1    the groups are independent, which is what two different
%                     neurons look like -> KEEP BOTH
%
%     ratio far above the groups fire at the same instant, which no pair of
%     1               neurons does: the same event detected twice, or a shared
%                     artifact -> DROP the weaker group
%
%   The same statistic therefore answers all three questions, which is why the
%   method does not need a separate rule per case.
%
%   out = jrc_autosplit(prmFile)
%   out = jrc_autosplit(prmFile, 'units', [875 536], 'maxK', 4)
%
%   Options
%     'units'         units to examine (default: every unit whose violating
%                     fraction exceeds minViolRatio)
%     'maxK'          how far to over-split before merging back (default 4).
%                     Over-splitting then merging is deliberate: it is far
%                     easier to recognise that two groups belong together than
%                     to guess the right number of groups up front.
%     'refractoryMs'  interval defining a violation (default 1.0; see
%                     jrc_relabel for why 1 ms and not 2)
%     'minViolRatio'  which units to examine (default 0.15)
%     'mergeRatio'    cross-correlogram ratio below which two groups are
%                     judged one neuron (default 0.5)
%     'dupRatio'      ratio above which one group is judged a duplicate of the
%                     other (default 5.0)
%     'minGroup'      smallest group worth keeping (default 30)
%     'maxSpikes'     subsample above this many spikes per unit, for speed
%                     (default 60000; clustering 200k points to decide a split
%                     is not worth the wait, and the modes are visible in a
%                     sample)
%     'outputDir'     where the _res.mat and .jrc files live
%
%   Reports only. Nothing is written: apply the splits in the GUI once the
%   recommendations look right.

    opts = struct('units', [], 'maxK', 4, 'refractoryMs', 1.0, ...
                  'minViolRatio', 0.15, 'mergeRatio', 0.5, 'dupRatio', 5.0, ...
                  'minGroup', 30, 'maxSpikes', 60000, 'outputDir', '');
    for k = 1:2:numel(varargin)
        name = varargin{k};
        if ~isfield(opts, name)
            error('jrc_autosplit:badOption', 'unknown option ''%s''', name);
        end
        opts.(name) = varargin{k+1};
    end

    resFile = jrc_res_path(prmFile, opts.outputDir);
    fprintf('Loading %s\n', resFile);
    res = load(resFile, 'spikeTimes','spikeClusters','spikeAmps','spikeSites', ...
                        'clusterSites','unitCount','featuresShape');
    site = double(res.spikeSites(:));
    if ~isfield(res, 'featuresShape')
        error('jrc_autosplit:noShape', ...
              ['this _res.mat has no featuresShape, so the feature file ' ...
               'cannot be read; re-import with jrc_import_cpp']);
    end

    fs = sampleRateOf(prmFile);
    t    = double(res.spikeTimes(:)) / fs;
    clu  = double(res.spikeClusters(:));
    amp  = abs(double(res.spikeAmps(:)));
    nClusters = max(clu);

    featFile = [regexprep(resFile, '_res\.mat$', '') '_features.jrc'];
    if ~isfile(featFile)
        error('jrc_autosplit:noFeatures', '%s not found', featFile);
    end
    featShape = double(res.featuresShape);
    fmap = memmapfile(featFile, 'Format', {'single', featShape, 'x'});
    fprintf('  features %s\n\n', mat2str(featShape));

    % ------------------------------------------------------- pick candidates
    units = opts.units(:)';
    if isempty(units)
        units = [];
        for c = 1:nClusters
            tk = sort(t(clu == c));
            if numel(tk) < 2*opts.minGroup, continue; end
            isi = diff(tk) * 1000;
            if sum(isi < opts.refractoryMs) / max(sum(isi <= 20),1) > opts.minViolRatio
                units(end+1) = c; %#ok<AGROW>
            end
        end
    end
    fprintf('examining %d units\n\n', numel(units));

    out = struct('unit', {}, 'nSpikes', {}, 'nNeurons', {}, 'nDropped', {}, ...
                 'groupSizes', {}, 'verdict', {});

    for c = units
        idx = find(clu == c);
        n0 = numel(idx);
        if n0 < 2*opts.minGroup, continue; end

        % Subsample for the clustering decision, but keep the full index list
        % so the reported group sizes refer to the real unit.
        useIdx = idx;
        if n0 > opts.maxSpikes
            useIdx = idx(round(linspace(1, n0, opts.maxSpikes)));
        end

        tk = t(useIdx);
        ak = amp(useIdx);
        sk = site(useIdx);
        [tk, ord] = sort(tk);
        ak = ak(ord);
        sk = sk(ord);
        sortedUse = useIdx(ord);

        % --- first split: by DETECTION SITE ---------------------------------
        % Done before any feature clustering, because the commonest reason a
        % unit fires below its refractory period is the same event crossing
        % threshold on two sites that were never merged - and that split is
        % invisible in the PCA features, which are computed per site and so
        % put both copies in the same place. Clustering on features alone
        % spreads the duplicates evenly over the groups, and the between-group
        % test then sees nothing: unit 875, which is 89% duplicates, came out
        % as four clean "neurons" that way.
        lbl = zeros(numel(tk), 1);
        sites = unique(sk);
        counts = arrayfun(@(s) sum(sk == s), sites);
        bigSites = sites(counts >= opts.minGroup);
        for g = 1:numel(bigSites)
            lbl(sk == bigSites(g)) = g;
        end
        if isempty(bigSites), lbl(:) = 1; end
        % Spikes on sites too rare to judge join the largest group rather than
        % being discarded: they are a handful of stragglers, not evidence.
        if any(lbl == 0)
            [~, biggest] = max(arrayfun(@(g) sum(lbl == g), 1:max(lbl)));
            lbl(lbl == 0) = biggest;
        end

        [lbl, notes] = resolveGroups(tk, ak, lbl, opts);

        % --- second split: by feature, within each surviving group ----------
        % Now that duplicates are gone, look for two genuine cells sharing one
        % site. Each group is split and resolved on its own, so a clean group
        % cannot be disturbed by a messy one.
        nextLabel = max(lbl) + 1;
        for g = unique(lbl(lbl > 0))'
            inG = find(lbl == g);
            if numel(inG) < 2*opts.minGroup, continue; end

            F = double(squeeze(fmap.Data.x(:, 1, sortedUse(inG))))';
            if size(F,1) ~= numel(inG), continue; end
            F = (F - mean(F,1)) ./ max(std(F,0,1), eps);

            k = min(opts.maxK, max(2, floor(numel(inG)/opts.minGroup)));
            warnState = warning('off','stats:kmeans:FailedToConvergeRep');
            sub = kmeans(F, k, 'Replicates', 3, 'MaxIter', 200, 'Start', 'plus');
            warning(warnState);

            subLbl = zeros(numel(inG),1);
            for j = 1:k, subLbl(sub == j) = nextLabel + j - 1; end
            [subLbl, subNotes] = resolveGroups(tk(inG), ak(inG), subLbl, opts);

            lbl(inG) = subLbl;
            nextLabel = nextLabel + k;
            notes = [notes, subNotes]; %#ok<AGROW>
        end

        alive = unique(lbl(lbl > 0));
        sizes = arrayfun(@(g) sum(lbl == g), alive);
        keep = sizes >= opts.minGroup;
        alive = alive(keep); sizes = sizes(keep);

        scale = n0 / numel(useIdx);
        verdict = sprintf('%d neuron(s)', numel(alive));
        if numel(alive) == 0, verdict = 'all noise/duplicate'; end

        fprintf('unit %4d  site %3d  %7d spikes -> %s\n', ...
                c, res.clusterSites(c), n0, verdict);
        for g = 1:numel(alive)
            fprintf('      group %d : %7d spikes (%.0f%%)\n', g, ...
                    round(sizes(g)*scale), 100*sizes(g)/numel(useIdx));
        end
        for s = notes, fprintf('      %s\n', s{1}); end

        out(end+1) = struct('unit', c, 'nSpikes', n0, 'nNeurons', numel(alive), ...
                            'nDropped', sum(lbl == 0), ...
                            'groupSizes', round(sizes*scale), ...
                            'verdict', verdict); %#ok<AGROW>
    end

    if isempty(out), fprintf('\nnothing to split\n'); return; end

    fprintf('\n---------------------------------------------\n');
    fprintf('units examined      : %d\n', numel(out));
    fprintf('  resolve to 2+ cells : %d\n', sum([out.nNeurons] >= 2));
    fprintf('  resolve to 1 cell   : %d\n', sum([out.nNeurons] == 1));
    fprintf('  entirely rubbish    : %d\n', sum([out.nNeurons] == 0));
    fprintf(['\nNothing written. Split the multi-cell units in the GUI, and\n' ...
             'let jrc_relabel condemn the ones that came out empty.\n']);
end

% ---------------------------------------------------------------------------
function [lbl, notes] = resolveGroups(tk, ak, lbl, opts)
    %RESOLVEGROUPS Merge, drop or keep each pair of groups, repeatedly.
    %
    %   Runs to a fixed point rather than in one pass, because merging two
    %   groups changes that group's relationship with every other one.
    notes = {};
    for iter = 1:10
        groups = unique(lbl(lbl > 0))';
        if numel(groups) < 2, break; end

        bestPair = []; bestRatio = inf; action = '';
        for a = 1:numel(groups)
            for b = a+1:numel(groups)
                ta = tk(lbl == groups(a));
                tb = tk(lbl == groups(b));
                if numel(ta) < opts.minGroup || numel(tb) < opts.minGroup, continue; end

                r = crossRatio(ta, tb, opts.refractoryMs);
                if ~isfinite(r), continue; end

                if r < opts.mergeRatio && r < bestRatio
                    bestRatio = r; bestPair = [groups(a) groups(b)]; action = 'merge';
                elseif r > opts.dupRatio && isempty(bestPair)
                    bestRatio = r; bestPair = [groups(a) groups(b)]; action = 'drop';
                end
            end
        end

        if isempty(bestPair), break; end

        if strcmp(action, 'merge')
            lbl(lbl == bestPair(2)) = bestPair(1);
            notes{end+1} = sprintf(...
                'merged two groups (cross-correlogram dips, ratio %.2f): one cell', ...
                bestRatio); %#ok<AGROW>
        else
            % Drop the weaker of the two: the duplicate is the copy with the
            % smaller amplitude, which is how JRCLUST resolves duplicates
            % within its merge radius.
            m1 = median(ak(lbl == bestPair(1)));
            m2 = median(ak(lbl == bestPair(2)));
            loser = bestPair(1); if m1 > m2, loser = bestPair(2); end
            notes{end+1} = sprintf(...
                'dropped %d spikes firing at zero lag with another group (ratio %.0f): duplicate', ...
                sum(lbl == loser), bestRatio); %#ok<AGROW>
            lbl(lbl == loser) = 0;
        end
    end
end

function r = crossRatio(ta, tb, refracMs)
    %CROSSRATIO Cross-correlogram density at short lag against 10-25 ms.
    %
    %   1.0 means the two trains are independent. Near 0 means they avoid each
    %   other at short lag, which only happens when they share a refractory
    %   period and are therefore one cell. Far above 1 means they fire
    %   together, which means one event detected twice.
    r = NaN;
    ta = sort(ta(:)); tb = sort(tb(:));
    if numel(ta) < 20 || numel(tb) < 20, return; end

    lo = refracMs / 1000;
    near = countWithin(ta, tb, lo) ;
    far  = countWithin(ta, tb, 0.025) - countWithin(ta, tb, 0.010);
    if far <= 0, return; end

    % near window is +/- refracMs wide, far window is +/- 15 ms wide
    r = (near / far) * (15 / refracMs);
end

function c = countWithin(ta, tb, win)
    %COUNTWITHIN Pairs (a,b) with |a - b| <= win. Both inputs sorted.
    edges = [-inf; tb; inf];
    [~, ~, hiBin] = histcounts(ta + win, edges);
    [~, ~, loBin] = histcounts(ta - win, edges);
    c = sum(hiBin - loBin);
end

function fs = sampleRateOf(prmFile)
    fs = [];
    txt = fileread(prmFile);
    tok = regexp(txt, '(?m)^\s*sampleRate\s*=\s*([0-9.eE+-]+)', 'tokens', 'once');
    if ~isempty(tok), fs = str2double(tok{1}); end
    if isempty(fs) || ~isfinite(fs) || fs <= 0
        error('jrc_autosplit:noSampleRate', 'could not read sampleRate from %s', prmFile);
    end
end
