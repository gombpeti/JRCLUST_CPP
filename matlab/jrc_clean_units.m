function jrc_clean_units(prmFile, varargin)
%JRC_CLEAN_UNITS Walk every cluster and strip the spikes that are not its own.
%
%   jrc_clean_units('mysession.prm', 'dryRun', true)
%   jrc_clean_units('mysession.prm')
%
%   Population-level artifact removal only catches events big enough to move
%   the rate of the WHOLE probe. An artifact that lands in a handful of units
%   stays invisible there, yet it can be a hundred times those units' own
%   firing rate. This works per unit, against each unit's own baseline.
%
%   Two passes:
%
%   1. BURST REMOVAL (per spike). A neuron's firing rate wanders with
%      behaviour but does not leap by orders of magnitude for a few seconds.
%      Bins where a unit fires far above its own median are artifact, and the
%      spikes in them are set to cluster 0.
%
%   2. UNIT VERDICT (per cluster). After cleaning, units whose spike train
%      still cannot come from one neuron are annotated 'to_delete':
%        - refractory violations too frequent to be physical
%        - a PEAK at zero lag in the autocorrelogram where a neuron must dip
%        - nothing much left after burst removal
%
%   Nothing is deleted here. Run jrc_prune afterwards to apply it.
%
%   Options:
%     'binSeconds'    bin width for the per-unit rate (default 10)
%     'madThresh'     bins above median + this many MADs are bursts (default 10)
%     'minBurstRatio' a bin must also exceed this multiple of the unit's median
%                     (default 8), so a near-silent unit cannot flag itself
%     'maxISIRatio'   #(ISI<=2ms)/#(ISI<=20ms) above this condemns a unit
%                     (default 0.45; a clean unit is near 0.02, a strong
%                     burster stays under about 0.2)
%     'maxZeroLag'    autocorrelogram near/far ratio above this condemns a
%                     unit (default 3.0). The ratio is the height of the
%                     ACG's centre bar against its 10-25 ms shoulder: 0 is a
%                     clean refractory dip, 1.0 is a flat ACG, 200 is an
%                     artifact. The default is CALIBRATED, not theoretical -
%                     measured over a 2.6 h 512-site recording, units clean by
%                     every other measure (ISI<=0.2, SNR>=4, >=5000 spikes,
%                     n=239) reach 2.05 at p99, while clearly contaminated
%                     units (ISI>=0.6, n=57) start at 3.19 at p1. The gap
%                     between those two populations is empty, so the cut goes
%                     there. The previous default of 1.0 sat inside the good
%                     population and condemned 342 of 1200 units and 23% of
%                     all spikes, including a 197,474-spike unit firing
%                     steadily at 21 Hz that scored 1.074.
%     'minSpikes'     units left with fewer spikes than this are condemned
%                     (default 30)
%     'outputDir'     where the _res.mat lives
%     'dryRun'        report only (default false)
%     'clearPrevious' clear existing 'to_delete' labels before re-judging
%                     (default true), so that re-running after loosening a
%                     threshold can give units back rather than only ever
%                     condemning more. Set false to keep labels you set by
%                     hand in the GUI and only add to them.

    opts = struct('binSeconds', 10, 'madThresh', 10, 'minBurstRatio', 8, ...
                  'maxISIRatio', 0.45, 'maxZeroLag', 3.0, 'minSpikes', 30, ...
                  'outputDir', '', 'dryRun', false, 'clearPrevious', true);
    for k = 1:2:numel(varargin)
        name = varargin{k};
        if ~isfield(opts, name)
            error('jrc_clean_units:badOption', 'unknown option ''%s''', name);
        end
        opts.(name) = varargin{k+1};
    end

    resFile = jrc_res_path(prmFile, opts.outputDir);
    fprintf('Loading %s\n', resFile);
    res = load(resFile);

    hCfg = jrclust.Config(prmFile);
    t = double(res.spikeTimes(:)) / hCfg.sampleRate;
    clu = double(res.spikeClusters(:));
    nClusters = max(clu);

    fprintf('  %d spikes, %d clusters\n\n', numel(t), nClusters);

    if ~isfield(res, 'clusterNotes') || numel(res.clusterNotes) < nClusters
        res.clusterNotes = repmat({''}, nClusters, 1);
    end
    notes = res.clusterNotes(:);

    % Clear any previous 'to_delete' before judging again.
    %
    % Without this the function can only ever condemn more units: a unit
    % marked under an earlier, stricter threshold would keep its label however
    % well it scores now, so re-running after correcting a threshold could not
    % give anything back. Re-judging from a clean slate makes the result depend
    % only on the current settings.
    %
    % This also clears 'to_delete' set by hand or by jrc_mark_noise. Pass
    % 'clearPrevious', false to keep those and only add to them.
    if opts.clearPrevious
        wasDeleted = strcmp(notes, 'to_delete');
        if any(wasDeleted)
            fprintf(['cleared %d existing ''to_delete'' labels before ' ...
                     're-judging\n  (pass ''clearPrevious'', false to keep ' ...
                     'them)\n'], sum(wasDeleted));
            notes(wasDeleted) = {''};
        end
    end

    drop = false(numel(t), 1);
    verdict = repmat({''}, nClusters, 1);

    nBurstUnits = 0; nBurstSpikes = 0;
    reasonISI = 0; reasonACG = 0; reasonEmpty = 0;

    for c = 1:nClusters
        idx = find(clu == c);
        if isempty(idx), continue; end

        tc = sort(t(idx));

        % ---- pass 1: bins where the unit far exceeds its own baseline -----
        span = tc(end) - tc(1);
        if span > opts.binSeconds && numel(tc) >= 20
            edges = tc(1):opts.binSeconds:(tc(end) + opts.binSeconds);
            counts = histcounts(tc, edges);

            med = median(counts(counts > 0));
            if isempty(med) || med <= 0, med = 1; end
            mad = median(abs(counts - med)) * 1.4826;
            if mad <= 0, mad = 1; end

            cutoff = max(med + opts.madThresh*mad, opts.minBurstRatio*med);
            hot = counts > cutoff;

            if any(hot)
                lo = edges(hot);
                hi = lo + opts.binSeconds;
                inBurst = false(size(idx));
                for b = 1:numel(lo)
                    inBurst = inBurst | (t(idx) >= lo(b) & t(idx) < hi(b));
                end

                if any(inBurst)
                    drop(idx(inBurst)) = true;
                    nBurstUnits = nBurstUnits + 1;
                    nBurstSpikes = nBurstSpikes + sum(inBurst);
                end
            end
        end

        % ---- pass 2: can what remains come from one neuron? ---------------
        keptIdx = idx(~drop(idx));
        tk = sort(t(keptIdx));

        if numel(tk) < opts.minSpikes
            verdict{c} = 'to_delete'; reasonEmpty = reasonEmpty + 1;
            continue;
        end

        isi = diff(tk) * 1000;
        isiRatio = sum(isi <= 2) / max(sum(isi <= 20), 1);
        if isiRatio > opts.maxISIRatio
            verdict{c} = 'to_delete'; reasonISI = reasonISI + 1;
            continue;
        end

        % Zero-lag autocorrelogram. Counts pairs within 2 ms against pairs in
        % 10-25 ms, scaled for the unequal window widths. A neuron's own
        % refractory period empties the near lags.
        zl = zeroLagRatio(tk);
        if isfinite(zl) && zl > opts.maxZeroLag
            verdict{c} = 'to_delete'; reasonACG = reasonACG + 1;
        end
    end

    % ------------------------------------------------------------- report
    condemned = find(strcmp(verdict, 'to_delete'));
    fprintf('pass 1 - burst removal (per unit, against its own baseline):\n');
    fprintf('  units with bursts   : %d / %d\n', nBurstUnits, nClusters);
    fprintf('  spikes removed      : %d (%.2f%% of all)\n', nBurstSpikes, ...
            100*nBurstSpikes/numel(t));

    fprintf('\npass 2 - unit verdicts on what remains:\n');
    fprintf('  condemned           : %d / %d\n', numel(condemned), nClusters);
    fprintf('    refractory violations too frequent : %d\n', reasonISI);
    fprintf('    autocorrelogram peaks at zero lag  : %d\n', reasonACG);
    fprintf('    too few spikes left                : %d\n', reasonEmpty);

    spikesInCondemned = sum(ismember(clu, condemned) & ~drop);
    fprintf('  their surviving spikes: %d\n', spikesInCondemned);
    fprintf('\n  total spikes affected : %d (%.2f%%)\n', ...
            nBurstSpikes + spikesInCondemned, ...
            100*(nBurstSpikes + spikesInCondemned)/numel(t));

    if opts.dryRun
        fprintf('\ndry run: nothing written\n');
        return;
    end

    % --------------------------------------------------------------- apply
    clu(drop) = 0;
    res.spikeClusters = int32(clu);
    notes(condemned) = {'to_delete'};
    res.clusterNotes = notes;

    nC = max(clu);
    if isfield(res, 'spikesByCluster')
        res.spikesByCluster = arrayfun(@(iC) find(clu == iC), 1:nC, ...
                                       'UniformOutput', 0)';
        res.unitCount = cellfun(@numel, res.spikesByCluster);
    end

    fprintf('\nWriting %s\n', resFile);
    jrclust.utils.saveStruct(res, resFile);
    fprintf('\nDone. To apply:\n    jrc_prune(''%s'')\n', prmFile);
end

%% ------------------------------------------------------------------ helpers

function r = zeroLagRatio(tk)
    %ZEROLAGRATIO Near-lag versus mid-lag pair counts in the autocorrelogram.
    %
    %   A real neuron cannot fire twice inside its refractory period, so pairs
    %   separated by under 2 ms should be far rarer, per unit of lag, than
    %   pairs at 10-25 ms. Scaling by the window widths makes 1.0 the
    %   "no refractory period at all" point.
    r = NaN;
    n = numel(tk);
    if n < 20, return; end

    near = 0; far = 0;
    stop = searchsortedRight(tk, tk + 0.025);
    for i = 1:n
        for j = i+1:stop(i)
            d = (tk(j) - tk(i)) * 1000;
            if d <= 2
                near = near + 1;
            elseif d >= 10 && d <= 25
                far = far + 1;
            end
        end
    end

    if far == 0, return; end
    r = (near / far) * (15 / 2);   % 15 ms of far window vs 2 ms of near
end

function idx = searchsortedRight(sorted, values)
    %SEARCHSORTEDRIGHT Index of the first element strictly greater than value.
    idx = zeros(numel(values), 1);
    j = 1;
    n = numel(sorted);
    for i = 1:numel(values)
        if j < i, j = i; end
        while j <= n && sorted(j) <= values(i)
            j = j + 1;
        end
        idx(i) = j - 1;
    end
end
