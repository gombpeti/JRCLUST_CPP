function jrc_relabel(prmFile, varargin)
%JRC_RELABEL Re-judge which units are noise, without touching any spike.
%
%   Recomputes the two spike-train tests on every unit and rewrites
%   clusterNotes from scratch. Nothing else in the result is modified: spike
%   assignments, waveforms and features are left exactly as they are.
%
%   This exists separately from jrc_clean_units because that function also
%   removes burst spikes, moving them to cluster 0 permanently. Re-running it
%   to correct a threshold would zero a further slice of the recording as a
%   side effect of fixing a label. Relabelling has to be free to run as often
%   as a threshold is adjusted, so it is kept free of side effects.
%
%   Previous 'to_delete' labels are cleared before judging. Without that the
%   result could only ever grow: a unit condemned under a stricter setting
%   would keep its label however well it scores now, and loosening a threshold
%   would give nothing back.
%
%   jrc_relabel(prmFile)
%   jrc_relabel(prmFile, 'dryRun', true)
%   jrc_relabel(prmFile, 'maxZeroLag', 3.0, 'maxISIRatio', 0.45)
%
%   Options
%     'maxZeroLag'    autocorrelogram near/far ratio above which a unit is
%                     condemned (default 3.0).
%
%                     The ratio is the height of the ACG's centre bar against
%                     its 10-25 ms shoulder: pairs of spikes under 2 ms apart
%                     against pairs 10-25 ms apart, scaled by the two window
%                     widths so that 1.0 means a perfectly flat ACG.
%
%                     0        a clean refractory dip - one neuron
%                     1 - 2    a neuron, or several neurons mixed: independent
%                              cells have a flat cross-correlation, so merging
%                              them flattens the ACG towards 1 but no further
%                     > 10     the same event counted more than once at the
%                              same instant - duplicate detection or artifact
%
%                     The default is CALIBRATED, not theoretical. Measured on
%                     a 2.6 h 512-site recording, units clean by every other
%                     measure (ISI <= 0.2, SNR >= 4, >= 5000 spikes, n = 239)
%                     reach 2.05 at p99, while clearly contaminated units
%                     (ISI >= 0.6, n = 57) start at 3.19 at p1. The gap
%                     between the two populations is empty, so the cut goes
%                     there. A cut at the theoretical 1.0 sits inside the good
%                     population and condemns real neurons.
%
%     'refractoryMs'  interval below which a second spike is impossible
%                     (default 1.0 ms).
%
%                     Set at the ABSOLUTE refractory period, the floor imposed
%                     by sodium-channel inactivation, rather than at the 1.5-2
%                     ms window used for published quality metrics (Hill,
%                     Mehta & Kleinfeld 2011 use 2 ms; Allen Institute and
%                     Kilosort/Phy use 1.5 ms). The wider windows overlap what
%                     real cells do: hippocampal complex-spike bursts reach
%                     intra-burst intervals of 3-6 ms, thalamic bursts 2-4 ms,
%                     and fast-spiking interneurons briefly reach 1.25-2 ms.
%                     Condemning a unit on intervals a burster can legitimately
%                     produce is not defensible; condemning it on intervals
%                     below 1 ms is, because no cell can fire there at all.
%
%                     Measured on this recording the distinction turned out not
%                     to change the verdicts much - of 178 units condemned by a
%                     2 ms window, the median had 61% of its intervals under
%                     1 ms - but the 1 ms cut is the one that cannot be argued
%                     with.
%
%     'maxISIRatio'   #(ISI<refractoryMs)/#(ISI<=20ms) above this condemns a
%                     unit (default 0.30). A clean unit is near 0; a unit above
%                     this is firing below the refractory period so often that
%                     it must hold more than one source.
%     'minSpikes'     units with fewer spikes than this are condemned
%                     (default 30)
%     'clearPrevious' clear existing 'to_delete' before judging (default
%                     true). Set false to keep labels you set by hand in the
%                     GUI and only add to them.
%     'rejudgeNoise'  re-judge the units the C++ classifier called 'noise'
%                     instead of inheriting its verdict (default true). Those
%                     that pass become 'review', so they are kept but
%                     flagged. Set false to trust its labels - but note that
%                     converting them straight to 'to_delete' on this
%                     recording removed 106 units and 1,138,750 spikes that
%                     pass every calibrated test, including one with a
%                     sub-millisecond rate of 0.001.
%     'outputDir'     where the _res.mat lives
%     'dryRun'        report only, write nothing (default false)
%
%   After a real run, apply it with:
%       jrc_prune(prmFile)

    opts = struct('maxZeroLag', 3.0, 'maxISIRatio', 0.30, 'refractoryMs', 1.0, ...
                  'minSpikes', 30, 'marginFactor', 1.5, 'clearPrevious', true, ...
                  'rejudgeNoise', false, 'minHalfWidthMs', 0.06, ...
                  'maxShankCoverage', 0, 'loudSNR', 10, 'loudMaxSpikes', 1000, 'loudMinCrossings', 2, ...
                  'maxCrossings', 12, 'maxSNR', 50, 'minTemplateScore', 0, ...
                  'outputDir', '', 'dryRun', false);
    for k = 1:2:numel(varargin)
        name = varargin{k};
        if ~isfield(opts, name)
            error('jrc_relabel:badOption', 'unknown option ''%s''', name);
        end
        opts.(name) = varargin{k+1};
    end

    resFile = jrc_res_path(prmFile, opts.outputDir);
    fprintf('Loading %s\n', resFile);
    res = load(resFile);

    % Probe geometry, for the shank-coverage test below. Read from the .prm
    % text rather than through jrclust.Config, for the reason given in
    % sampleRateOf: one malformed entry anywhere in the file would otherwise
    % stop the whole config constructing.
    prmText = fileread(prmFile);
    shankMap = readVector(prmText, 'shankMap');
    ignored  = readVector(prmText, 'ignoreSites');

    fs = sampleRateOf(prmFile);
    t   = double(res.spikeTimes(:)) / fs;
    clu = double(res.spikeClusters(:));
    nClusters = max(clu);
    fprintf('  %d spikes, %d clusters, %.1f s at %g Hz\n\n', ...
            numel(t), nClusters, max(t), fs);

    if ~isfield(res, 'clusterNotes') || numel(res.clusterNotes) < nClusters
        res.clusterNotes = repmat({''}, nClusters, 1);
    end
    notes = res.clusterNotes(:);
    notes = notes(1:nClusters);

    if opts.clearPrevious
        wasDeleted = strcmp(notes, 'to_delete');
        if any(wasDeleted)
            fprintf(['Clearing %d existing ''to_delete'' labels before ' ...
                     're-judging.\n  (pass ''clearPrevious'', false to keep ' ...
                     'them)\n\n'], sum(wasDeleted));
            notes(wasDeleted) = {''};
        end
    end

    % The C++ classifier's own 'noise' verdicts are re-judged here rather than
    % inherited, and the units it condemned are demoted to 'review' instead of
    % to 'to_delete' when they pass these tests.
    %
    % They cannot simply be trusted. That classifier's thresholds were set
    % from first principles, and the ones checked against real data have been
    % wrong in the same direction every time. Converting its 96 'noise'
    % verdicts straight to 'to_delete' on this recording deleted 106 units
    % holding 1,138,750 spikes - 85% of everything removed - and among them
    % were units with a sub-millisecond rate of 0.001 and a zero-lag of 0.019,
    % cleaner than units being kept. Its burst test in particular condemns
    % episodic firing, which a place cell or a task-modulated neuron does
    % legitimately.
    %
    % Its shape and spread tests do see things a spike train cannot, so a unit
    % it flagged is still worth a look - hence 'review' rather than silence.
    cppNoise = false(size(notes));
    if opts.rejudgeNoise
        cppNoise = strcmp(notes, 'noise');
        if any(cppNoise)
            fprintf(['Re-judging %d units the C++ classifier called ''noise''.\n' ...
                     '  Those that pass will be marked ''review'', not deleted.\n' ...
                     '  (pass ''rejudgeNoise'', false to keep its verdicts)\n\n'], ...
                    sum(cppNoise));
            notes(cppNoise) = {''};
        end
    end

    % ------------------------------------------------------------- measure
    % Group the spikes once rather than masking the full 17-million-element
    % vector per cluster, which would otherwise dominate the runtime.
    byCluster = cell(nClusters, 1);
    order = clu > 0;
    [sortedClu, ix] = sort(clu(order));
    tt = t(order);
    tt = tt(ix);
    edges = [0; find(diff(sortedClu)); numel(sortedClu)];
    for k = 1:numel(edges)-1
        c = sortedClu(edges(k)+1);
        byCluster{c} = tt(edges(k)+1 : edges(k+1));
    end

    % --- shape, from the mean raw waveform -------------------------------
    % These answer what a spike train cannot. A one-sample trough has a
    % perfectly respectable autocorrelogram, so the refractory and ACG tests
    % pass it; only its WIDTH gives it away. Likewise a deflection appearing
    % across a whole shank at one amplitude can have clean statistics while
    % being common-mode pickup rather than a cell.
    %
    % Both tests were set from units judged by eye in the GUI: clusters with
    % a 0.050 ms trough (one sample at 20 kHz) or covering half their shank
    % were called noise, while a 0.300 ms trough on five sites was a real
    % neuron the earlier settings had wrongly deleted.
    halfWidthMs = nan(nClusters, 1);
    shankCov    = nan(nClusters, 1);
    crossings   = nan(nClusters, 1);
    if isfield(res, 'meanWfGlobalRaw') && isfield(res, 'clusterSites')
        wf = res.meanWfGlobalRaw;
        sites = double(res.clusterSites(:));
        for c = 1:min(nClusters, size(wf,3))
            pk = sites(c);
            if pk < 1 || pk > size(wf,2), continue; end
            w = double(wf(:,:,c));
            wpk = w(:,pk);
            [trough, it] = min(wpk);
            if trough >= 0, continue; end

            half = trough / 2;
            L = it;  while L > 1 && wpk(L-1) <= half, L = L - 1; end
            Rr = it; while Rr < numel(wpk) && wpk(Rr+1) <= half, Rr = Rr + 1; end
            halfWidthMs(c) = (Rr - L + 1) * 1000 / fs;

            % Smoothness: how often the mean waveform crosses its own mean.
            %
            % A real unit leaves baseline once, troughs, and returns - two
            % crossings - and stays smooth even when built from few spikes,
            % because the spikes align. A cluster of unrelated transients
            % averages to a jagged trace that wanders across its own mean
            % repeatedly: measured here, units judged noise by eye scored 4
            % to 20 while a genuine 88-spike neuron scored 2.
            %
            % Counted against the plain mean with no dead band. A band of
            % 10% of Vpp was tried and destroyed the measurement, because Vpp
            % is set by one large excursion while the ringing being counted
            % is far smaller - every unit then read 1-2, noise and neurons
            % alike.
            crossings(c) = sum(abs(diff(sign(wpk - mean(wpk)))) > 0);

            if ~isempty(shankMap) && pk <= numel(shankMap)
                amp = max(w,[],1) - min(w,[],1);
                onShank = find(shankMap == shankMap(pk));
                onShank = setdiff(onShank, ignored);
                if ~isempty(onShank) && max(amp) > 0
                    shankCov(c) = sum(amp(onShank) >= 0.5*max(amp)) / numel(onShank);
                end
            end
        end
    end

    % --- template match, derived from this recording with no labels -------
    % See buildTemplate below. Reported always; condemns only if
    % minTemplateScore is set above 0, which it is NOT by default.
    [tmplScore, tmplOK, tmplWhy] = templateScores(res, nClusters, fs);

    nSpk = zeros(nClusters, 1);
    isiR = nan(nClusters, 1);
    zlR  = nan(nClusters, 1);

    fprintf('Measuring %d units', nClusters);
    for c = 1:nClusters
        if mod(c, 100) == 0, fprintf('.'); end
        tk = sort(byCluster{c});
        nSpk(c) = numel(tk);
        if nSpk(c) < 2, continue; end

        isi = diff(tk) * 1000;                       % ms
        isiR(c) = sum(isi < opts.refractoryMs) / max(sum(isi <= 20), 1);

        % The autocorrelogram test keeps its own 2 ms near window rather than
        % following refractoryMs. Its threshold was calibrated against that
        % window, and the scaling factor that makes 1.0 mean "flat" depends on
        % the window width, so narrowing it here would silently invalidate the
        % calibration. The two tests answer different questions anyway: this
        % one asks whether the unit fires where it physically cannot, the ACG
        % asks whether the same event has been counted twice.
        zlR(c)  = zeroLagRatio(tk);
    end
    fprintf(' done\n\n');

    % ------------------------------------------------------------- judge
    tooFew   = nSpk < opts.minSpikes;
    badISI   = ~tooFew & isiR > opts.maxISIRatio;
    badACG   = ~tooFew & ~badISI & isfinite(zlR) & zlR > opts.maxZeroLag;
    % Trough narrower than a membrane potential can move and return.
    %
    % The default of 0.06 ms catches a trough spanning ONE sample at 20 kHz
    % and nothing wider. At 0.1 ms the same test fired on 190 of 755 units, a
    % quarter of the recording, which means it was measuring something other
    % than what it claims to: a real trough spans 0.2-0.4 ms, but the 5th
    % percentile here is already 0.100, so a cut at 0.1 sits inside the
    % ordinary population. 0.06 sits below all of it.
    badWidth = opts.minHalfWidthMs > 0 & ~tooFew & ...
               isfinite(halfWidthMs) & halfWidthMs <= opts.minHalfWidthMs;

    % Shank coverage cannot work on this data, and stays off (threshold 0).
    % meanWfGlobalRaw is nSamples x nSites x nClusters, but JRCLUST fills
    % only the unit's own site neighbourhood and leaves every other site at
    % zero - 9 sites of 512 here. No coverage or spread figure computed from
    % it can exceed the neighbourhood size, whatever threshold is set. The
    % C++ classifier's maxSpreadSites = 20 was unreachable for exactly this
    % reason, which is why it never fired once on a real recording.
    badShank = opts.maxShankCoverage > 0 & ~tooFew & ...
               isfinite(shankCov) & shankCov > opts.maxShankCoverage;

    % Loud but nearly silent.
    %
    % SNR and firing rate are loosely coupled for a real neuron: a cell close
    % enough to the electrode to produce a large SNR is close enough to be
    % detected often. A unit with SNR 22 that fires 164 times in 2.6 hours is
    % not a quiet neuron - it is a handful of large transients that clustered
    % together, and the waveform says so: individual spikes from 95 to
    % 3046 uV, a mean waveform of only 48 uV because they do not align, and
    % the largest amplitude sitting on a NEIGHBOURING site rather than its
    % own.
    %
    % This is the test that catches what the spike-train tests cannot. A
    % cluster of large transients separated by more than a millisecond has a
    % clean refractory period and a clean autocorrelogram; only its shape and
    % its implausible rate give it away.
    %
    % Deliberately narrow: across 755 units it condemns 13 holding 2374
    % spikes, 0.015% of the recording, so it removes the degenerate cases
    % without touching genuinely sparse units of ordinary amplitude.
    if isfield(res, 'unitSNR') && numel(res.unitSNR) >= nClusters
        snr = double(res.unitSNR(:));
        snr = snr(1:nClusters);
    else
        snr = nan(nClusters, 1);
    end
    % The third condition is what keeps real quiet neurons. A large, rarely
    % firing unit whose waveform is CLEAN is simply a sparse cell, and
    % deleting it would be the same error as before: an 88-spike unit at SNR
    % 11.7 with a 0.200 ms trough and two crossings was condemned by the
    % first two conditions alone, and it is a neuron. Requiring the waveform
    % to be non-smooth as well separates it from units of similar amplitude
    % and rate that score 4 to 8 crossings.
    badLoud = ~tooFew & isfinite(snr) & isfinite(crossings) & ...
              snr > opts.loudSNR & nSpk < opts.loudMaxSpikes & ...
              crossings > opts.loudMinCrossings;

    % Ringing, on its own. Well above the 95th percentile of this recording
    % (10), so it reaches only the clearly oscillatory waveforms rather than
    % the ordinary spread.
    badRing = opts.maxCrossings > 0 & ~tooFew & ...
              isfinite(crossings) & crossings > opts.maxCrossings;

    % Amplitude beyond anything biological, decisive on its own.
    %
    % Needed because the smoothness condition on badLoud, which protects
    % genuinely quiet neurons, also shields artifacts whose mean waveform
    % happens to be smooth. A unit with SNR 160 and 81 spikes scored only 2
    % crossings and walked through every other test.
    %
    % Safe because the two populations are far apart: measured over this
    % recording, units clean by every independent measure reach SNR 29.4 at
    % most, while the units above 50 pin near 1000 uV - the amplifier rail -
    % or sit at 100 to 270. The cut sits in empty space between them, so it
    % cannot reach a real cell however quiet or however close to the
    % electrode.
    badSNR = opts.maxSNR > 0 & ~tooFew & isfinite(snr) & snr > opts.maxSNR;

    % Template match. OFF by default (threshold 0), because one template
    % assumes one waveform family and the populations still overlap - see
    % templateScores. Reported either way so the score can be judged.
    badTmpl = opts.minTemplateScore > 0 & tmplOK & ~tooFew & ...
              isfinite(tmplScore) & tmplScore < opts.minTemplateScore;

    failed   = tooFew | badISI | badACG | badWidth | badShank | badLoud | badRing | badSNR | badTmpl;

    % Separate units that fail decisively from units that fail by a hair.
    %
    % Any hard threshold has a boundary, and units sitting on it are the ones a
    % small change in the cut flips either way - so they are precisely the ones
    % that should not be deleted unseen. On this recording 49 of 295 condemned
    % units failed by under half a threshold, and those 49 held 46% of all the
    % spikes being removed, because the marginal cases are large units rather
    % than fragments. Deleting them silently is how a 197,474-spike neuron was
    % lost at a zero-lag of 1.074 against a cut of 1.0.
    %
    % Decisive failures are annotated 'to_delete' and jrc_prune removes them.
    % Marginal ones are annotated 'review' instead: jrc_prune leaves those
    % alone, and the GUI's "display annotated units only" brings them up as a
    % short worklist.
    % A shape failure is decisive on its own. There is no "marginally one
    % sample wide": the trough either spans more than a couple of samples or
    % it does not, and a deflection on half a shank is not a near miss.
    decisive = tooFew | badWidth | badShank | badLoud | badRing | badSNR | badTmpl | ...
               (isiR > opts.marginFactor * opts.maxISIRatio) | ...
               (isfinite(zlR) & zlR > opts.marginFactor * opts.maxZeroLag);

    condemn = failed & decisive;
    review  = (failed & ~decisive) | (cppNoise(:) & ~failed);

    spikesIn = @(m) sum(nSpk(m));
    total = sum(nSpk);

    fprintf(['Verdicts (refractoryMs %.2f, maxISIRatio %.2f, ' ...
             'maxZeroLag %.2f, minSpikes %d)\n'], ...
            opts.refractoryMs, opts.maxISIRatio, opts.maxZeroLag, opts.minSpikes);
    fprintf('  too few spikes                     : %4d units\n', sum(tooFew));
    fprintf('  fires below the refractory period  : %4d units\n', sum(badISI));
    fprintf('  autocorrelogram peaks at zero lag  : %4d units\n', sum(badACG));
    fprintf('  trough too narrow to be a spike    : %4d units\n', sum(badWidth));
    fprintf('  covers too much of its own shank   : %4d units\n', sum(badShank));
    fprintf('  large amplitude but almost silent  : %4d units\n', sum(badLoud));
    fprintf('  waveform rings instead of spiking  : %4d units\n', sum(badRing));
    fprintf('  amplitude beyond anything biologic : %4d units\n', sum(badSNR));
    if opts.minTemplateScore > 0
        fprintf('  waveform unlike a spike (template) : %4d units\n', sum(badTmpl));
    end

    % Always reported, whether or not it condemns anything: the score is the
    % most informative single number about a unit's shape, and seeing its
    % distribution is how the threshold gets chosen.
    if tmplOK
        good = tmplScore(isfinite(tmplScore));
        fprintf(['\n  template match (no labels used; %s)\n' ...
                 '    pct      5     25     50     75     95\n' ...
                 '    val %6.2f %6.2f %6.2f %6.2f %6.2f\n'], ...
                tmplWhy, prctile(good,[5 25 50 75 95]));
        if opts.minTemplateScore <= 0
            fprintf(['    not condemning on this - pass ''minTemplateScore'', 0.85\n' ...
                     '    to enable it once the cut has been judged\n']);
        end
    else
        fprintf('\n  template match unavailable: %s\n', tmplWhy);
    end
    fprintf('  ------------------------------------------------\n');
    fprintf('  failed a test                      : %4d of %d units\n\n', ...
            sum(failed), nClusters);

    fprintf('  -> ''to_delete'' (decisive)           : %4d units, %d spikes (%.2f%%)\n', ...
            sum(condemn), spikesIn(condemn), 100*spikesIn(condemn)/max(total,1));
    fprintf('  -> ''review''    (within %.2gx of a cut): %4d units, %d spikes (%.2f%%)\n', ...
            opts.marginFactor, sum(review), spikesIn(review), ...
            100*spikesIn(review)/max(total,1));
    fprintf('  -> kept unannotated                 : %4d units, %d spikes\n', ...
            sum(~failed), spikesIn(~failed));
    fprintf(['\n  jrc_prune removes only ''to_delete''. The ''review'' units stay,\n' ...
             '  and the GUI''s "display annotated units only" lists them.\n']);

    reviewList = find(review);
    [~, ord] = sort(nSpk(reviewList), 'descend');
    show = reviewList(ord(1:min(15, numel(ord))));
    if ~isempty(show)
        fprintf('\n  marked ''review'', largest first:\n');
        fprintf('  %6s %10s %8s %10s %10s  %s\n', ...
                'unit', 'spikes', 'site', 'zeroLag', 'isiRatio', 'fails');
        for c = show'
            site = NaN;
            if isfield(res, 'clusterSites') && numel(res.clusterSites) >= c
                site = res.clusterSites(c);
            end
            f = '';
            if badISI(c),   f = [f 'ISI ']; end
            if badACG(c),   f = [f 'ACG ']; end
            if badWidth(c), f = [f 'narrow ']; end
            if badRing(c),  f = [f 'ringing ']; end
            if badSNR(c),   f = [f 'SNR ']; end
            if badLoud(c),  f = [f 'loud-but-silent']; end
            fprintf('  %6d %10d %8d %10.3f %10.3f  %s\n', ...
                    c, nSpk(c), site, zlR(c), isiR(c), f);
        end
    end

    % The units nearest the cut, so the threshold can be judged by eye rather
    % than taken on trust.
    nearCut = find(~failed & isfinite(zlR));
    [~, ord] = sort(zlR(nearCut), 'descend');
    show = nearCut(ord(1:min(8, numel(ord))));
    if ~isempty(show)
        fprintf('\n  closest kept units to the zero-lag cut:\n');
        fprintf('  %6s %10s %8s %10s %10s\n', 'unit', 'spikes', 'site', 'zeroLag', 'isiRatio');
        for c = show'
            site = NaN;
            if isfield(res, 'clusterSites') && numel(res.clusterSites) >= c
                site = res.clusterSites(c);
            end
            fprintf('  %6d %10d %8d %10.3f %10.3f\n', c, nSpk(c), site, zlR(c), isiR(c));
        end
    end

    condemnedList = find(condemn);
    [~, ord] = sort(nSpk(condemnedList), 'descend');
    show = condemnedList(ord(1:min(8, numel(ord))));
    if ~isempty(show)
        fprintf('\n  largest condemned units:\n');
        fprintf('  %6s %10s %8s %10s %10s\n', 'unit', 'spikes', 'site', 'zeroLag', 'isiRatio');
        for c = show'
            site = NaN;
            if isfield(res, 'clusterSites') && numel(res.clusterSites) >= c
                site = res.clusterSites(c);
            end
            fprintf('  %6d %10d %8d %10.3f %10.3f\n', c, nSpk(c), site, zlR(c), isiR(c));
        end
    end

    if opts.dryRun
        fprintf('\ndry run: nothing written\n');
        return;
    end

    % ------------------------------------------------------------- apply
    notes(condemn) = {'to_delete'};
    notes(review)  = {'review'};
    res.clusterNotes = notes;

    fprintf('\nWriting %s\n', resFile);
    jrclust.utils.saveStruct(res, resFile);
    fprintf('\nDone. Inspect in the GUI with\n    jrc manual %s\n', prmFile);
    fprintf('or apply the deletions with\n    jrc_prune(''%s'')\n', prmFile);
end

% ---------------------------------------------------------------------------
function [score, ok, why] = templateScores(res, nClusters, fs)
    %TEMPLATESCORES How spike-shaped each unit is, judged without labels.
    %
    %   Answers the objection "how can you call something not-a-neuron if you
    %   were never shown a neuron" in two steps, neither of which needs
    %   labelled examples:
    %
    %     1. A CANDIDATE template comes from the recording itself. JRCLUST
    %        centres each spike window on the detected peak, so most units
    %        already carry their trough at the same sample. Normalise each
    %        unit's mean waveform to its own trough depth - leaving shape
    %        only, with amplitude, gain and probe divided out - and take the
    %        per-sample median of the aligned ones. Real units outnumber
    %        noise by better than ten to one, so the median is theirs.
    %
    %     2. PHYSIOLOGY decides whether that candidate is a spike, using
    %        constraints from membrane biophysics rather than precedent: one
    %        dominant excursion, a duration no membrane can beat, a
    %        repolarisation, a return to baseline. If the candidate fails,
    %        the recording is majority noise and this refuses to score rather
    %        than silently inverting itself.
    %
    %   Because every recording builds its own template, nothing transfers
    %   between sessions or animals - which is exactly why it generalises. A
    %   different probe, filter or preparation produces a different template,
    %   correctly, since its own neurons define what normal looks like there.
    %
    %   KNOWN LIMITATION: one template assumes one waveform family. A
    %   legitimate but atypical cell - a fast-spiking interneuron with a
    %   narrow spike, say - scores low for being unusual rather than for
    %   being fake. On a measured recording, units of 85,000 to 220,000
    %   spikes at SNR 6.7 to 12.9 scored 0.72 to 0.76 while hand-judged noise
    %   scored 0.25 to 0.77, so the two overlap. That is why nothing is
    %   condemned on this score by default.
    score = nan(nClusters,1); ok = false; why = 'no waveforms available';
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

    centre  = mode(troughAt(~isnan(troughAt)));
    aligned = troughAt == centre;
    if sum(aligned) < 20, why = 'too few aligned units'; return; end
    tmpl = median(raw(:, aligned), 2);

    % -- physiology check on the candidate ------------------------------
    [tr, it] = min(tmpl);
    reb  = max(tmpl(it:end));
    half = tr/2;
    L=it; while L>1 && tmpl(L-1)<=half, L=L-1; end
    R=it; while R<nS && tmpl(R+1)<=half, R=R+1; end
    halfW = (R-L+1)*1000/fs;

    % Significant excursions, NOT baseline crossings. The raw window spans
    % about 2 ms, so most of it is flat tail; counting crossings there
    % measures sample noise and fails good templates.
    big = abs(tmpl) > 0.25;
    nExc = sum(diff([false; big(:)]) == 1);

    ok = nExc <= 2 && abs(tr) > abs(reb) && halfW >= 0.1 && halfW <= 1.5 && ...
         reb > 0.05 && abs(median(tmpl(min(nS,it+14):end))) < 0.35;
    if ~ok
        why = sprintf(['template is not spike-shaped (%d excursions, ' ...
                       'half-width %.3f ms) - recording may be majority noise'], ...
                      nExc, halfW);
        return;
    end
    why = sprintf('template OK: half-width %.3f ms, rebound %.2f, from %d units', ...
                  halfW, reb, sum(aligned));

    % -- score every unit, best match over a shift ----------------------
    tc = tmpl - mean(tmpl);
    maxShift = 8;
    for c = 1:nClusters
        if all(isnan(raw(:,c))), continue; end
        best = -1;
        for s = -maxShift:maxShift
            v = circshift(raw(:,c), s);
            keep = (1+max(0,s)) : (nS+min(0,s));
            vv = v(keep) - mean(v(keep));
            tt = tc(keep) - mean(tc(keep));
            d = norm(vv)*norm(tt);
            if d > 0, best = max(best, (vv'*tt)/d); end
        end
        score(c) = best;
    end
end

function v = readVector(prmText, name)
    %READVECTOR A bracketed numeric vector from the .prm text, or empty.
    v = [];
    tok = regexp(prmText, ['(?s)' name '\s*=\s*\[(.*?)\]'], 'tokens', 'once');
    if isempty(tok), return; end
    try
        v = str2num(tok{1}); %#ok<ST2NM>
    catch
        v = [];
    end
    v = v(:)';
end

function fs = sampleRateOf(prmFile)
    %SAMPLERATEOF Sample rate, read straight from the .prm text.
    %
    %   Deliberately avoids jrclust.Config: a malformed entry anywhere in the
    %   file (an outputDir written as a cell, for instance) makes the whole
    %   config fail to construct, and relabelling needs nothing from it but
    %   this one number.
    fs = [];
    txt = fileread(prmFile);
    tok = regexp(txt, '(?m)^\s*sampleRate\s*=\s*([0-9.eE+-]+)', 'tokens', 'once');
    if ~isempty(tok), fs = str2double(tok{1}); end
    if isempty(fs) || ~isfinite(fs) || fs <= 0
        error('jrc_relabel:noSampleRate', ...
              'could not read sampleRate from %s', prmFile);
    end
end

function r = zeroLagRatio(tk)
    %ZEROLAGRATIO Near-lag versus mid-lag pair counts in the autocorrelogram.
    %
    %   Pairs under 2 ms apart against pairs 10-25 ms apart, scaled by the two
    %   window widths so that 1.0 is a flat ACG. A neuron's refractory period
    %   empties the near lags and drives this towards 0.
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
    %SEARCHSORTEDRIGHT Index of the last element not greater than value.
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
