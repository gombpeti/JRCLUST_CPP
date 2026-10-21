function updateFigCorr(obj)
    %UPDATEFIGCORR Plot auto- and cross-correlograms for the selection
    if isempty(obj.selected) || ~obj.hasFig('FigCorr')
        return;
    end

    plotFigCorr(obj.hFigs('FigCorr'), obj.hClust, obj.hCfg, obj.selected);
end

%% LOCAL FUNCTIONS
function hFigCorr = plotFigCorr(hFigCorr, hClust, hCfg, selected)
    %PLOTFIGCORR Correlogram for one unit, a pair, or an N x N matrix
    %
    %   One unit  : its autocorrelogram.
    %   Two units : their cross-correlogram, as before.
    %   N units   : an N x N grid, autocorrelograms down the diagonal and
    %               cross-correlograms off it.
    %
    %   The matrix is the useful case when deciding whether several units are
    %   one cell. A refractory dip in a CROSS-correlogram means the two
    %   trains never fire together within a millisecond, which only one
    %   refractory period explains - so they are one neuron split in two.
    %   A flat cross-correlogram means they are independent cells, however
    %   alike their waveforms. Waveform similarity cannot separate those two
    %   cases; the correlogram can, and seeing every pair at once is what
    %   makes a group of four or five decidable in one look.
    jitterMs = 0.5;   % bin width
    nLagsMs  = 25;    % show +/- 25 ms

    jitterSamp = round(jitterMs*hCfg.sampleRate/1000);
    nLags = round(nLagsMs/jitterMs);
    lagTime = (-nLags:nLags)*jitterMs;

    nSel = numel(selected);

    % Each unit's own colour, as the waveform view draws it, so a diagonal
    % panel can be matched to its unit at a glance. Only the diagonal is
    % coloured: those are the autocorrelograms and each belongs to exactly
    % one unit. An off-diagonal panel belongs to a PAIR, so colouring it
    % after either one would be misleading - it stays neutral.
    colors = jrclust.views.selectionColors(selected, hClust.nClusters);

    % Spike times once per unit, in jitter bins.
    binned = cell(nSel, 1);
    for k = 1:nSel
        binned{k} = int32(double(hClust.spikeTimes(hClust.spikesByCluster{selected(k)}))/jitterSamp);
    end

    % ---- one or two units: the original single plot --------------------
    if nSel <= 2
        hFigCorr.clf();

        if nSel == 1
            counts = correlogram(binned{1}, binned{1}, nLags, true);
            titleStr = sprintf('Unit %d', selected(1));
        else
            counts = correlogram(binned{1}, binned{2}, nLags, false);
            titleStr = sprintf('Unit %d vs. Unit %d', selected(1), selected(2));
        end

        hFigCorr.addAxes('default');
        hFigCorr.addPlot('hBar', @bar, lagTime, counts, 1);
        hFigCorr.axApply('default', @xlabel, 'Time (ms)');
        hFigCorr.axApply('default', @ylabel, 'Counts');
        hFigCorr.axApply('default', @grid, 'on');
        hFigCorr.axApply('default', @set, 'YScale', 'log');
        hFigCorr.axApply('default', @title, titleStr);
        hFigCorr.axApply('default', @set, 'XLim', jitterMs*[-nLags, nLags]);
        return;
    end

    % ---- N x N matrix ---------------------------------------------------
    % Rebuilt from scratch each time rather than updated in place: the grid
    % changes shape whenever the number of selected units changes, so there
    % is nothing stable to update.
    hFigCorr.clf();

    pad = 0.055;
    w = (1 - pad*1.6) / nSel;
    h = (1 - pad*1.6) / nSel;

    for r = 1:nSel
        for c = 1:nSel
            key = sprintf('cc%d_%d', r, c);
            pos = [pad + (c-1)*w, 1 - pad*0.75 - r*h, w*0.84, h*0.74];

            % addAxes returns the handle, and addPlot needs it as its FIRST
            % argument: without it addPlot falls back to obj.gca(), so every
            % bar lands in one default axes spanning the whole figure and the
            % panels draw on top of each other.
            hAx = hFigCorr.addAxes(key, 'Position', pos);

            isAuto = (r == c);
            counts = correlogram(binned{r}, binned{c}, nLags, isAuto);

            hFigCorr.addPlot(['b' key], @bar, hAx, lagTime, counts, 1);

            % The diagonal belongs to exactly one unit, so its bars take that
            % unit's colour - the same one the waveform view uses. An
            % off-diagonal panel belongs to a PAIR, so colouring it after
            % either unit would be misleading; those stay neutral.
            if isAuto
                hFigCorr.plotApply(['b' key], @set, ...
                                   'FaceColor', colors(r, :), ...
                                   'EdgeColor', colors(r, :));
            end

            hFigCorr.axApply(key, @set, 'XLim', jitterMs*[-nLags, nLags], ...
                             'YScale', 'log', 'FontSize', 7);
            hFigCorr.axApply(key, @grid, 'on');

            % Every panel says which pair it is, so a cell can be read
            % without counting rows and columns back to an edge label.
            if isAuto
                hFigCorr.axApply(key, @title, sprintf('Unit %d', selected(r)), ...
                                 'FontSize', 8, 'FontWeight', 'bold');
                % tinted, so self-comparisons are obvious at a glance
                hFigCorr.axApply(key, @set, 'Color', [0.95 0.95 0.88]);
            else
                hFigCorr.axApply(key, @title, ...
                                 sprintf('%d vs %d', selected(r), selected(c)), ...
                                 'FontSize', 8);
            end

            if r == nSel
                hFigCorr.axApply(key, @xlabel, 'ms', 'FontSize', 7);
            else
                hFigCorr.axApply(key, @set, 'XTickLabel', []);
            end
        end
    end
end

function counts = correlogram(iTimes, jTimes, nLags, isAuto)
    %CORRELOGRAM Counts of jTimes + lag coinciding with iTimes.
    %
    %   For a cross-correlogram the reference train is widened by one bin
    %   either side, as the original code did, so that a coincidence falling
    %   just across a bin edge is still counted.
    %
    %   The zero-lag bin of an autocorrelogram is skipped: every spike
    %   trivially coincides with itself and would swamp the plot.
    if ~isAuto
        iTimes = [iTimes, iTimes - 1, iTimes + 1];
    end

    lags = -nLags:nLags;
    counts = zeros(size(lags));
    for k = 1:numel(lags)
        if isAuto && lags(k) == 0
            continue;
        end
        counts(k) = numel(intersect(iTimes, jTimes + lags(k)));
    end
end
