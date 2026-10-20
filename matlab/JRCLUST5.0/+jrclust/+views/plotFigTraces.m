function tracesFilt = plotFigTraces(hFigTraces, hCfg, tracesRaw, resetAxis, hClust)
    %PLOTFIGTRACES Plot raw traces view, one subplot panel per visible shank
    hBox = [];
    if hCfg.getOr('showTraceProgressBox', 0)
        hBox = jrclust.utils.qMsgBox('Plotting...', 0, 1);
    end

    hFigTraces.wait(1);

    sampleRate = hCfg.sampleRate / hCfg.nSkip;
    viSamples1 = 1:hCfg.nSkip:size(tracesRaw, 2);
    evtWindowSamp = round(hCfg.evtWindowSamp / hCfg.nSkip); %show 2x of range

    if strcmpi(hFigTraces.figData.filter, 'on')
        % back up old settings
        sampleRateOld = hCfg.sampleRate;
        filterTypeOld = hCfg.filterType;

        % temporarily alter settings to send traces through a filter
        hCfg.sampleRate = sampleRate;
        hCfg.useGPU = 0;
        hCfg.filterType = hCfg.dispFilter;

        if hCfg.fftThresh > 0
            % fftClean expects samples x channels; tracesRaw here is sites x samples, so transpose in and back out
            tracesRaw = jrclust.filters.fftClean(tracesRaw', hCfg.fftThresh, hCfg)';
        end

        % filtCAR/applyCAR expect samples x channels; tracesRaw here is sites x samples, so transpose in and back out
        tracesFilt = jrclust.filters.filtCAR(tracesRaw(:, viSamples1)', [], [], 0, hCfg);
        tracesFilt = jrclust.utils.bit2uV(tracesFilt', hCfg);
        filterToggle = hCfg.filterType;

        % restore old settings
        hCfg.sampleRate = sampleRateOld;
        hCfg.useGPU = 1;
        hCfg.filterType = filterTypeOld;
    else
        tracesFilt = jrclust.utils.meanSubtract(single(tracesRaw(:, viSamples1))) * hCfg.bitScaling;
        filterToggle = 'off';
    end

    multiEdges = [];
    if hCfg.nSegmentsTraces == 1
        XData = ((hFigTraces.figData.windowBounds(1):hCfg.nSkip:hFigTraces.figData.windowBounds(end))-1) / hCfg.sampleRate;
        XLabel = 'Time (s)';
    else
        XData = (0:(size(tracesFilt, 2) - 1)) / (hCfg.sampleRate / hCfg.nSkip) + (hFigTraces.figData.windowBounds(1)-1) / hCfg.sampleRate;
        [multiBounds, multiRange, multiEdges] = jrclust.views.sampleSkip(hFigTraces.figData.windowBounds, hFigTraces.figData.nSamplesTotal, hCfg.nSegmentsTraces);

        tlim_show = (cellfun(@(x) x(1), multiBounds([1, end]))) / hCfg.sampleRate;
        XLabel = sprintf('Time (s), %d segments merged (%0.1f ~ %0.1f s, %0.2f s each)', hCfg.nSegmentsTraces, tlim_show, diff(hCfg.dispTimeLimits));

        mrX_edges = XData(repmat(multiEdges(:)', [3, 1]));
        csTime_bin = cellfun(@(x) sprintf('%0.1f', x(1)/hCfg.sampleRate), multiBounds, 'UniformOutput', 0);
    end

    traceColor = [1 1 1]*.88;
    textColor = [1 1 1]*.9;
    if isfield(hFigTraces.figData, 'traceColor')
        traceColor = hFigTraces.figData.traceColor;
    end
    if isfield(hFigTraces.figData, 'textColor')
        textColor = hFigTraces.figData.textColor;
    end

    hFigTraces.figApply(@set, 'Color', 'k');

    shankOrder = hFigTraces.figData.shankOrder;
    shankVisible = hFigTraces.figData.shankVisible;
    shankAxes = hFigTraces.figData.shankAxes;
    visibleFullIdx = find(shankVisible);

    % remove spike overlays from other clusters/threads (hTrace_*/hEdges_* line objects are reused
    % across redraws via updatePlot; they self-clean if a panel rebuild invalidates their axes)
    plotKeys = keys(hFigTraces.hPlots);
    chSpk = plotKeys(startsWith(plotKeys, 'chSpk'));
    if ~isempty(chSpk)
        cellfun(@(pk) hFigTraces.rmPlot(pk), chSpk);
    end

    if strcmpi(hFigTraces.figData.spikes, 'on') && ~isempty(hClust)
        recPos = find(strcmp(hFigTraces.figData.hRec.rawPath, hCfg.rawRecordings));
        if recPos == 1
            offset = 0;
        else
            % find all recordings coming before hRec and sum up nSamples for each
            hRecs = arrayfun(@(iRec) jrclust.detect.newRecording(hCfg.rawRecordings{iRec}, hCfg), 1:(recPos-1), 'UniformOutput', 0);
            offset = sum(cellfun(@(hR) hR.nSamples, hRecs));
        end

        recTimes = double(hClust.spikeTimes) - double(offset);
        tStart = single(hFigTraces.figData.windowBounds(1) - 1)/hCfg.sampleRate;

        if hCfg.nSegmentsTraces > 1
            spikesInRange = inRange(recTimes, multiBounds);
            allSpikeSites = single(hClust.spikeSites(spikesInRange));
            allSpikeTimes = double(recTimes(spikesInRange));
            allSpikeTimes = round(whereMember(allSpikeTimes, multiRange) / hCfg.nSkip);
        else
            spikesInRange = recTimes >= hFigTraces.figData.windowBounds(1) & recTimes < hFigTraces.figData.windowBounds(end);
            allSpikeSites = single(hClust.spikeSites(spikesInRange));
            allSpikeTimes = double(recTimes(spikesInRange));
            allSpikeTimes = round((allSpikeTimes - hCfg.sampleRate*tStart) / hCfg.nSkip); % time offset
        end

        allSpikeClusters = hClust.spikeClusters(spikesInRange);
    end

    for v = 1:numel(visibleFullIdx)
        shankID = shankOrder(visibleFullIdx(v));
        hAx = shankAxes(v);
        siteIdx = find(hCfg.shankMap == shankID);
        nSitesShank = numel(siteIdx);

        panelTraces = tracesFilt(siteIdx, :); % nSitesShank x nSamples

        % site 1 (most dorsal, first entry of siteIdx) should render at the TOP of the panel.
        % Rather than reversing the axes' YDir (which would also flip each waveform's polarity
        % upside down), we flip which Y-offset each site row gets: site row 1 -> offset
        % nSitesShank (top), last row -> offset 1 (bottom). Polarity is untouched either way.
        siteRowOffset = nSitesShank:-1:1;

        plotKeyTrace = sprintf('hTrace_%d', shankID);
        hFigTraces.multiplot(plotKeyTrace, hFigTraces.figData.maxAmp, XData, panelTraces', siteRowOffset, 0, hAx);
        hFigTraces.plotApply(plotKeyTrace, @set, 'Color', traceColor, 'LineWidth', .5, 'Visible', hFigTraces.figData.traces);

        set(hAx, 'Color', 'k', 'XColor', textColor, 'YColor', textColor, ...
                 'GridColor', [1 1 1]*.35, 'MinorGridColor', [1 1 1]*.2, ...
                 'GridAlpha', .35, 'MinorGridAlpha', .25, 'Box', 'on', ...
                 'YTick', 1:nSitesShank, 'YTickLabel', arrayfun(@num2str, flip(siteIdx), 'UniformOutput', 0));
        grid(hAx, hFigTraces.figData.grid);

        title(hAx, sprintf('Shank %d (%0.1f uV)', shankID, hFigTraces.figData.maxAmp), 'Color', textColor, 'FontWeight', 'normal');
        xlabel(hAx, XLabel, 'Color', textColor);
        if v == 1
            ylabel(hAx, 'Site #', 'Color', textColor);
        end

        if hCfg.nSegmentsTraces > 1
            mrY_edges = repmat([0; nSitesShank + 1; nan], 1, numel(multiEdges));
            plotKeyEdges = sprintf('hEdges_%d', shankID);
            if hFigTraces.hasPlot(plotKeyEdges)
                hFigTraces.plotApply(plotKeyEdges, @set, 'XData', mrX_edges(:), 'YData', mrY_edges(:));
            else
                hFigTraces.addPlot(plotKeyEdges, @line, hAx, mrX_edges(:), mrY_edges(:), 'Color', [1 1 1]*.65, 'LineWidth', 1);
            end
            set(hAx, 'XTick', XData(multiEdges), 'XTickLabel', csTime_bin);
        end

        % plot spikes for sites on this shank only
        if strcmpi(hFigTraces.figData.spikes, 'on') && ~isempty(hClust)
            [~, localSite] = ismember(double(allSpikeSites), siteIdx); % 0 if not on this shank
            onShank = localSite > 0;

            if isempty(hClust)
                for iLocal = 1:nSitesShank % deal with subsample factor
                    onSite = find(onShank & localSite == iLocal);
                    if isempty(onSite)
                        continue;
                    end

                    timesOnSite = allSpikeTimes(onSite);
                    [mrY11, mrX11] = vr2mr3_(panelTraces(iLocal, :), timesOnSite, evtWindowSamp); %display purpose x2
                    mrT11 = single(mrX11-1) / sampleRate + tStart;

                    plotKey = sprintf('chSpk%d_%d', shankID, iLocal);
                    hFigTraces.addPlot(plotKey, @line, hAx, nan, nan, 'Color', [1 0 0], 'LineWidth', 1.5);
                    hFigTraces.multiplot(plotKey, hFigTraces.figData.maxAmp, mrT11, mrY11, siteRowOffset(iLocal), 0, hAx);
                end
            else % batch spikes by cluster to keep the trace GUI responsive
                inRangeClusters = allSpikeClusters(onShank);
                onShankLocalSite = localSite(onShank);
                onShankTimes = allSpikeTimes(onShank);
                spikeColors = [lines(hClust.nClusters); 0 0 0];
                uniqueClusters = unique(inRangeClusters(inRangeClusters > 0));
                winLen = evtWindowSamp(end) - evtWindowSamp(1) + 1;

                for iClu = 1:numel(uniqueClusters)
                    iCluster = uniqueClusters(iClu);
                    clusterMask = inRangeClusters == iCluster;
                    clusterTimes = onShankTimes(clusterMask);
                    clusterLocalSites = onShankLocalSite(clusterMask);
                    nClusterSpikes = numel(clusterTimes);

                    % One line object per cluster, with NaNs separating waveforms.
                    batchX = nan(winLen + 1, nClusterSpikes);
                    batchY = nan(winLen + 1, nClusterSpikes);
                    for iSpike = 1:nClusterSpikes
                        iLocal = clusterLocalSites(iSpike);
                        [spikeTrace, spikeRange] = vr2mr3_( ...
                            panelTraces(iLocal, :), clusterTimes(iSpike), evtWindowSamp);
                        nSamples = numel(spikeRange);
                        batchX(1:nSamples, iSpike) = ...
                            double(spikeRange(:) - 1) / sampleRate + tStart;
                        batchY(1:nSamples, iSpike) = ...
                            double(spikeTrace(:)) / hFigTraces.figData.maxAmp + siteRowOffset(iLocal);
                    end

                    plotKey = sprintf('chSpk%d_%d', shankID, iCluster);
                    hFigTraces.addPlot(plotKey, @line, hAx, batchX(:), batchY(:), ...
                                       'Color', spikeColors(iCluster, :), 'LineWidth', 2);
                end
            end
        end
    end

    if resetAxis
        jrclust.views.resetFigTraces(hFigTraces, tracesRaw, hCfg);
    end

    hFigTraces.figApply(@set, 'Name', sprintf('%s: filter: %s', hCfg.configFile, filterToggle));
    hFigTraces.wait(0);

    jrclust.utils.tryClose(hBox);
end

%% LOCAL FUNCTIONS
function isInRange = inRange(vals, ranges)
    isInRange = false(size(vals));

    if ~iscell(ranges)
        ranges = {ranges};
    end

    for iRange = 1:numel(ranges)
        bounds = ranges{iRange};
        isInRange = isInRange | (vals >= bounds(1) & vals <= bounds(2));
    end
end

function loc = whereMember(needle, haystack)
    % needle = int32(needle); % ??
    % haystack = int32(haystack); % ??
    [~, loc] = ismember(int32(needle), int32(haystack));
end

function [mr, ranges] = vr2mr3_(traces, spikeTimes, evtWindow)
    % JJJ 2015 Dec 24
    % vr2mr2: quick version and doesn't kill index out of range
    % assumes vi is within range and tolerates spkLim part of being outside
    % works for any datatype

    % prepare indices
    spikeTimes = spikeTimes(:)';
    ranges = int32(bsxfun(@plus, (evtWindow(1):evtWindow(end))', spikeTimes));

    ranges(ranges < 1) = 1;
    ranges(ranges > numel(traces)) = numel(traces); %keep # sites consistent

    % build spike table
    nSpikes = numel(spikeTimes);
    mr = reshape(traces(ranges(:)), [], nSpikes);
end
