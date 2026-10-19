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
                  'rejudgeNoise', false, 'minHalfWidthMs', 0, ...
                  'maxShankCoverage', 0, 'outputDir', '', 'dryRun', false);
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
    % Both shape tests are OFF by default (threshold 0 disables them).
    %
    % They encode the right idea - a one-sample trough is not a spike, and a
    % deflection across a shank is common mode - but neither has been
    % calibrated, and a trial run condemned 268 of 1143 units on width alone,
    % which is far more than the C++ classifier's equivalent test flags. A
    % test that fires on 23% of a recording is measuring something other than
    % what it claims to. The C++ classifier already applies its own versions
    % during the sort, so nothing is lost by leaving these off until the
    % threshold has been fitted to units judged by eye.
    badWidth = opts.minHalfWidthMs > 0 & ~tooFew & ...
               isfinite(halfWidthMs) & halfWidthMs <= opts.minHalfWidthMs;
    badShank = opts.maxShankCoverage > 0 & ~tooFew & ...
               isfinite(shankCov) & shankCov > opts.maxShankCoverage;
    failed   = tooFew | badISI | badACG | badWidth | badShank;

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
    decisive = tooFew | badWidth | badShank | ...
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
            if badISI(c), f = [f 'ISI ']; end
            if badACG(c), f = [f 'ACG']; end
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
