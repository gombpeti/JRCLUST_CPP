function out = jrc_split_candidates(prmFile, varargin)
%JRC_SPLIT_CANDIDATES Tell units that should be SPLIT from units to delete.
%
%   A cluster firing below the absolute refractory period is not necessarily
%   rubbish. Two independent neurons can produce an interval of half a
%   millisecond between them - cell A fires, cell B fires shortly after -
%   without either cell violating its own refractory period. What is violated
%   is the assumption that the cluster holds one cell. Deleting such a unit
%   throws away two real neurons; splitting it recovers both.
%
%   The same refractory violation can equally come from one event counted
%   twice, when a spike crosses threshold on two sites that were never merged,
%   or from an artifact appearing across the probe at one instant. Those are
%   genuinely worthless and should go.
%
%   The two cases are separable by looking at WHAT is paired:
%
%     two neurons     the two spikes of a violating pair have DIFFERENT
%                     amplitudes, because they come from cells at different
%                     distances from the site, and the unit's amplitude
%                     distribution is bimodal
%
%     one event twice the two spikes have the SAME amplitude, because they are
%                     the same deflection, and the distribution is unimodal
%
%   This function measures both and reports a recommendation per unit. It
%   changes nothing: splitting is done in the GUI, where the result can be
%   seen before it is accepted.
%
%   out = jrc_split_candidates(prmFile)
%   out = jrc_split_candidates(prmFile, 'refractoryMs', 1.0)
%
%   Options
%     'refractoryMs'  interval below which a pair counts as a violation
%                     (default 1.0; see jrc_relabel for why 1 ms rather than
%                     the 1.5-2 ms used for published quality metrics)
%     'minViolRatio'  only examine units whose violating fraction exceeds this
%                     (default 0.30, matching jrc_relabel's maxISIRatio)
%     'minAsymmetry'  median amplitude difference within violating pairs, as a
%                     fraction of the larger amplitude, above which the pair
%                     is judged to hold two different cells (default 0.30)
%     'minSeparation' separation of the two amplitude modes in pooled standard
%                     deviations, above which a split is likely to succeed
%                     (default 2.0 - two sd apart is visibly distinct rather
%                     than a smooth tail being cut arbitrarily)
%     'outputDir'     where the _res.mat lives
%
%   Returns a struct array with one entry per examined unit.

    opts = struct('refractoryMs', 1.0, 'minViolRatio', 0.30, ...
                  'minAsymmetry', 0.30, 'minSeparation', 2.0, 'outputDir', '');
    for k = 1:2:numel(varargin)
        name = varargin{k};
        if ~isfield(opts, name)
            error('jrc_split_candidates:badOption', 'unknown option ''%s''', name);
        end
        opts.(name) = varargin{k+1};
    end

    resFile = jrc_res_path(prmFile, opts.outputDir);
    fprintf('Loading %s\n', resFile);
    res = load(resFile, 'spikeTimes','spikeClusters','spikeAmps','spikeSites', ...
                        'clusterSites','unitCount');

    fs = sampleRateOf(prmFile);
    t    = double(res.spikeTimes(:)) / fs;
    clu  = double(res.spikeClusters(:));
    amp  = abs(double(res.spikeAmps(:)));   % trough depth, sign removed
    site = double(res.spikeSites(:));
    nClusters = max(clu);

    fprintf('  %d spikes, %d clusters\n\n', numel(t), nClusters);

    out = struct('unit', {}, 'nSpikes', {}, 'violRatio', {}, 'asymmetry', {}, ...
                 'separation', {}, 'crossSite', {}, 'recommendation', {});

    fprintf('%6s %9s %8s %9s %9s %9s %8s  %s\n', ...
            'unit','spikes','site','viol','asym','sep','xsite','recommendation');

    for c = 1:nClusters
        idx = find(clu == c);
        if numel(idx) < 30, continue; end

        [tk, ord] = sort(t(idx));
        ak = amp(idx(ord));
        sk = site(idx(ord));

        isi = diff(tk) * 1000;
        viol = isi < opts.refractoryMs;
        nViol = sum(viol);
        violRatio = nViol / max(sum(isi <= 20), 1);
        if violRatio <= opts.minViolRatio, continue; end
        if nViol < 10, continue; end

        % --- are the two spikes of a violating pair the same size? ---------
        a1 = ak(viol);
        a2 = ak([false; viol(1:end)] == 1);     % the second spike of each pair
        n = min(numel(a1), numel(a2));
        a1 = a1(1:n); a2 = a2(1:n);
        denom = max(max(a1, a2), 1);
        asym = median(abs(a1 - a2) ./ denom);

        % Were the two detected on different sites? The same event crossing
        % threshold on two unmerged sites is a duplicate, not two cells.
        crossSite = mean(sk(viol) ~= sk([false; viol(1:end)] == 1));

        % --- is the unit's amplitude distribution bimodal? ------------------
        sep = twoModeSeparation(ak);

        if asym >= opts.minAsymmetry && sep >= opts.minSeparation
            rec = 'SPLIT - two cells';
        elseif crossSite > 0.5 && asym < opts.minAsymmetry
            rec = 'delete - same event on two sites';
        elseif asym < opts.minAsymmetry
            rec = 'delete - one event counted twice';
        else
            rec = 'inspect - mixed evidence';
        end

        fprintf('%6d %9d %8d %9.3f %9.3f %9.3f %8.2f  %s\n', ...
                c, numel(idx), res.clusterSites(c), violRatio, asym, sep, ...
                crossSite, rec);

        out(end+1) = struct('unit', c, 'nSpikes', numel(idx), ...
                            'violRatio', violRatio, 'asymmetry', asym, ...
                            'separation', sep, 'crossSite', crossSite, ...
                            'recommendation', rec); %#ok<AGROW>
    end

    if isempty(out)
        fprintf('\nno unit exceeds the violation threshold\n');
        return;
    end

    isSplit = strncmp({out.recommendation}, 'SPLIT', 5);
    fprintf('\n%d units examined: %d to split, %d to delete/inspect\n', ...
            numel(out), sum(isSplit), sum(~isSplit));
    fprintf('spikes in the split candidates: %d\n', sum([out(isSplit).nSpikes]));
    fprintf(['\nSplit them in the GUI: open jrc manual, select the unit, and\n' ...
             'use Split. Nothing has been changed here.\n']);
end

% ---------------------------------------------------------------------------
function sep = twoModeSeparation(a)
    %TWOMODESEPARATION Distance between two amplitude modes, in pooled sd.
    %
    %   Fits two clusters to the amplitudes by 1-D k-means and reports how far
    %   apart their centres are relative to their spread. Below about 2 there
    %   is one mode with a tail; above it there are two populations.
    sep = 0;
    a = a(isfinite(a));
    if numel(a) < 30, return; end

    lo = prctile(a, 25);
    hi = prctile(a, 75);
    if hi <= lo, return; end

    for iter = 1:50
        mid = 0.5 * (lo + hi);
        g1 = a <= mid;
        if ~any(g1) || all(g1), return; end
        newLo = mean(a(g1));
        newHi = mean(a(~g1));
        if abs(newLo - lo) < 1e-6 && abs(newHi - hi) < 1e-6
            lo = newLo; hi = newHi; break;
        end
        lo = newLo; hi = newHi;
    end

    g1 = a <= 0.5*(lo+hi);
    v1 = var(a(g1));
    v2 = var(a(~g1));
    pooled = sqrt(0.5 * (v1 + v2));
    if pooled > 0, sep = abs(hi - lo) / pooled; end
end

function fs = sampleRateOf(prmFile)
    %SAMPLERATEOF Sample rate, read straight from the .prm text.
    fs = [];
    txt = fileread(prmFile);
    tok = regexp(txt, '(?m)^\s*sampleRate\s*=\s*([0-9.eE+-]+)', 'tokens', 'once');
    if ~isempty(tok), fs = str2double(tok{1}); end
    if isempty(fs) || ~isfinite(fs) || fs <= 0
        error('jrc_split_candidates:noSampleRate', ...
              'could not read sampleRate from %s', prmFile);
    end
end
