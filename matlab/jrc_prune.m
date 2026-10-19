function jrc_prune(prmFile, varargin)
%JRC_PRUNE Write a copy with noise spikes removed, without the GUI.
%
%   jrc_prune('mydata.prm')
%
%   Drops every spike belonging to a unit annotated 'to_delete', plus anything
%   already at spikeClusters == 0, renumbers the surviving units, and writes a
%   complete pruned session to <prm dir>\pruned\ :
%
%       <name>_pruned.prm
%       <name>_pruned_res.mat
%       <name>_pruned_raw.jrc / _filt.jrc / _features.jrc / _hist.jrc
%
%   which opens with:  jrc manual <name>_pruned.prm
%
%   Options:
%     'outputDir'   where the input _res.mat and .jrc files live
%     'prunedDir'   where to write (default: <prm dir>\pruned)
%     'chunkSpikes' spikes per I/O block (default 50000; lower it if memory
%                   is tight, raise it for speed)
%     'dryRun'      report what would happen without writing
%
%   WHY THIS EXISTS
%   The GUI's Edit > Delete annotated runs out of memory on large sessions.
%   It subsets the waveform arrays with vals(:,:,indices), and spikesRaw for
%   15 million spikes is about 11 GB, so the subset needs another 11 GB
%   contiguous before anything is freed. This script never holds those arrays:
%   it streams the .jrc files a block of spikes at a time, so peak memory is
%   set by 'chunkSpikes' rather than by the size of the recording.

    opts = struct('outputDir', '', 'prunedDir', '', 'chunkSpikes', 50000, ...
                  'dryRun', false);
    for k = 1:2:numel(varargin)
        name = varargin{k};
        if ~isfield(opts, name)
            error('jrc_prune:badOption', 'unknown option ''%s''', name);
        end
        opts.(name) = varargin{k+1};
    end

    if ~exist(prmFile, 'file')
        error('jrc_prune:noPrm', '%s does not exist', prmFile);
    end

    hCfg = jrclust.Config(prmFile);

    [srcDir, sessionName] = fileparts(hCfg.configFile);
    if ~isempty(hCfg.outputDir), srcDir = hCfg.outputDir; end
    if ~isempty(opts.outputDir), srcDir = opts.outputDir; end

    resFile = fullfile(srcDir, [sessionName '_res.mat']);
    if ~exist(resFile, 'file')
        error('jrc_prune:noRes', '%s not found', resFile);
    end

    fprintf('Loading %s\n', resFile);
    res = load(resFile);

    nSpikes = numel(res.spikeClusters);
    spikeClusters = double(res.spikeClusters(:));

    % ------------------------------------------------------- what to drop
    nClustersOld = max(spikeClusters);
    toDelete = [];
    if isfield(res, 'clusterNotes') && ~isempty(res.clusterNotes)
        toDelete = find(strcmp(res.clusterNotes(:), 'to_delete'));
    end

    if isempty(toDelete) && ~any(spikeClusters == 0)
        fprintf('Nothing annotated ''to_delete'' and no spikes at cluster 0.\n');
        fprintf('Run jrc_mark_noise first, or there is nothing to prune.\n');
        return;
    end

    spikeClusters(ismember(spikeClusters, toDelete)) = 0;

    keepSpikes = spikeClusters > 0;
    nKeep = nnz(keepSpikes);
    nDrop = nSpikes - nKeep;

    fprintf('  units annotated to_delete : %d / %d\n', numel(toDelete), nClustersOld);
    fprintf('  spikes dropped            : %d / %d (%.2f%%)\n', ...
            nDrop, nSpikes, 100 * nDrop / nSpikes);
    fprintf('  spikes kept               : %d\n', nKeep);

    if nKeep == 0
        error('jrc_prune:empty', 'no spikes would survive; refusing to write');
    end

    % --------------------------------------------------------- destination
    prunedDir = opts.prunedDir;
    if isempty(prunedDir)
        [prmDir, ~] = fileparts(hCfg.configFile);
        prunedDir = fullfile(prmDir, 'pruned');
    end
    prunedName = [sessionName '_pruned'];

    fprintf('  writing to                : %s\n', prunedDir);

    if opts.dryRun
        fprintf('\ndry run: nothing written\n');
        return;
    end

    if exist(prunedDir, 'dir') ~= 7, mkdir(prunedDir); end

    prunedPrm      = fullfile(prunedDir, [prunedName '.prm']);
    prunedRes      = fullfile(prunedDir, [prunedName '_res.mat']);
    prunedRaw      = fullfile(prunedDir, [prunedName '_raw.jrc']);
    prunedFilt     = fullfile(prunedDir, [prunedName '_filt.jrc']);
    prunedFeatures = fullfile(prunedDir, [prunedName '_features.jrc']);
    prunedHist     = fullfile(prunedDir, [prunedName '_hist.jrc']);

    % ------------------------------------------------ stream the .jrc files
    % One block of spikes at a time, so peak memory does not scale with the
    % recording. This is the part the GUI cannot do.
    res.rawShape      = streamPrune(fullfile(srcDir, [sessionName '_raw.jrc']), ...
                                    prunedRaw, res, 'rawShape', 'int16', ...
                                    keepSpikes, opts.chunkSpikes, 'raw');
    res.filtShape     = streamPrune(fullfile(srcDir, [sessionName '_filt.jrc']), ...
                                    prunedFilt, res, 'filtShape', 'int16', ...
                                    keepSpikes, opts.chunkSpikes, 'filt');
    res.featuresShape = streamPrune(fullfile(srcDir, [sessionName '_features.jrc']), ...
                                    prunedFeatures, res, 'featuresShape', 'single', ...
                                    keepSpikes, opts.chunkSpikes, 'features');

    % ------------------------------------------------------ prune the table
    res.spikeClusters = spikeClusters;
    res = pruneSpikeFields(res, keepSpikes);
    res = compactClusters(res);
    res = rebuildClusterCenters(res);   % must follow compaction

    res.spikeClusters     = int32(res.spikeClusters);
    res.initialClustering = int32(res.spikeClusters);

    nClustersNew = double(max(res.spikeClusters));
    res.spikesByCluster = arrayfun(@(iC) find(res.spikeClusters == iC), ...
                                   1:nClustersNew, 'UniformOutput', 0)';
    res.unitCount = cellfun(@numel, res.spikesByCluster);

    % Rebuild the by-site tables. ORIENTATION MATTERS: these must be 1 x nSites
    % row cell arrays, or getClusterFeaturesSite errors when a cluster is
    % clicked in the GUI.
    nSites = hCfg.nSites;
    if isfield(res, 'spikeSites') && ~isempty(res.spikeSites)
        res.spikesBySite = arrayfun(@(s) int32(find(res.spikeSites == s)), ...
                                    1:nSites, 'UniformOutput', 0);
    end
    if isfield(res, 'spikeSites2') && ~isempty(res.spikeSites2)
        res.spikesBySite2 = arrayfun(@(s) int32(find(res.spikeSites2 == s)), ...
                                     1:nSites, 'UniformOutput', 0);
    end
    if isfield(res, 'spikeSites3') && ~isempty(res.spikeSites3)
        res.spikesBySite3 = arrayfun(@(s) int32(find(res.spikeSites3 == s)), ...
                                     1:nSites, 'UniformOutput', 0);
    else
        res.spikesBySite3 = cell(1, nSites);
    end

    % Fresh history: the pruned copy is a new starting point.
    res.editPos = 1;
    res.history = containers.Map('KeyType', 'int32', 'ValueType', 'char');
    res.history(1) = sprintf('pruned save; removed %d spikeClusters==0 spikes', nDrop);

    fidHist = fopen(prunedHist, 'w');
    if fidHist == -1
        error('jrc_prune:hist', 'cannot write %s', prunedHist);
    end
    cleaner = onCleanup(@() fclose(fidHist));
    fwrite(fidHist, int32(1), 'int32');
    fwrite(fidHist, int32(res.spikeClusters), 'int32');
    clear cleaner;

    % ------------------------------------------------------- pruned config
    hCfg.save(prunedPrm, 1);
    commentInvalidZeroNSitesExcl(prunedPrm);
    appendOutputDir(prunedPrm, prunedDir);

    fprintf('Writing %s\n', prunedRes);
    jrclust.utils.saveStruct(res, prunedRes);

    fprintf('\nDone: %d units, %d spikes (was %d units, %d spikes)\n', ...
            nClustersNew, nKeep, nClustersOld, nSpikes);
    fprintf('Open with:\n    jrc manual %s\n', prunedPrm);
end

%% ------------------------------------------------------------------ helpers

function shape = streamPrune(inFile, outFile, res, shapeField, dtype, ...
                             keepSpikes, chunkSpikes, label)
    %STREAMPRUNE Copy a .jrc file keeping only selected spikes, block by block.
    shape = [];

    if ~isfield(res, shapeField) || isempty(res.(shapeField))
        fprintf('  %-9s : no shape recorded, skipping\n', label);
        return;
    end
    if ~exist(inFile, 'file')
        fprintf('  %-9s : %s not found, skipping\n', label, inFile);
        return;
    end

    oldShape = double(res.(shapeField));
    perSpike = oldShape(1) * oldShape(2);
    nSpikes  = oldShape(3);

    if numel(keepSpikes) ~= nSpikes
        error('jrc_prune:shapeMismatch', ...
              '%s holds %d spikes but the table has %d', label, nSpikes, ...
              numel(keepSpikes));
    end

    fin = fopen(inFile, 'r');
    if fin == -1, error('jrc_prune:open', 'cannot read %s', inFile); end
    fout = fopen(outFile, 'w');
    if fout == -1
        fclose(fin);
        error('jrc_prune:open', 'cannot write %s', outFile);
    end
    cleaner = onCleanup(@() [fclose(fin), fclose(fout)]);

    nKept = 0;
    i = 1;
    while i <= nSpikes
        j = min(i + chunkSpikes - 1, nSpikes);
        count = j - i + 1;

        block = fread(fin, perSpike * count, ['*' dtype]);
        if numel(block) < perSpike * count
            error('jrc_prune:truncated', '%s ended early at spike %d', label, i);
        end

        block = reshape(block, perSpike, count);
        block = block(:, keepSpikes(i:j));

        if ~isempty(block)
            fwrite(fout, block, dtype);
            nKept = nKept + size(block, 2);
        end

        i = j + 1;
    end

    clear cleaner;

    shape = [oldShape(1), oldShape(2), nKept];
    fprintf('  %-9s : %d -> %d spikes  [%s]\n', label, nSpikes, nKept, ...
            mat2str(shape));
end

function res = pruneSpikeFields(res, keepSpikes)
    %PRUNESPIKEFIELDS Subset every field indexed by spike.
    nSpikes = numel(keepSpikes);
    names = fieldnames(res);

    for i = 1:numel(names)
        fn = names{i};
        val = res.(fn);

        % Skipped here and handled below: these hold spike INDICES, so they
        % need renumbering as well as subsetting, and doing the subset twice
        % would index a already-shortened array with a full-length mask.
        if any(strcmp(fn, {'ordRho', 'spikeNeigh', ...
                           'spikesRaw', 'spikesFilt', 'spikeFeatures'}))
            continue;
        elseif isvector(val) && numel(val) == nSpikes
            res.(fn) = val(keepSpikes);
        elseif isnumeric(val) && ismatrix(val) && size(val, 1) == nSpikes
            res.(fn) = val(keepSpikes, :);
        elseif isnumeric(val) && ismatrix(val) && size(val, 2) == nSpikes
            res.(fn) = val(:, keepSpikes);
        end
    end

    % Spike indices stored in other fields must be renumbered, or they point
    % at the wrong spikes once the gaps close up.
    oldToNew = zeros(nSpikes, 1, 'int32');
    oldToNew(keepSpikes) = int32(1:nnz(keepSpikes));

    if isfield(res, 'spikeNeigh') && ~isempty(res.spikeNeigh)
        neigh = double(res.spikeNeigh(:));
        neigh = neigh(keepSpikes);
        valid = neigh > 0 & neigh <= nSpikes;
        valid(valid) = keepSpikes(neigh(valid));

        % A spike whose nearest-denser neighbour was pruned points at itself,
        % the same convention assignClusters uses for a density maximum.
        %
        % It must NOT be left at 0: doAssign evaluates
        % spikeClusters(spikeNeigh), and 0 is not a valid MATLAB index, so
        % jrc recluster would error on the pruned file.
        out = int32(1:numel(neigh))';
        out(valid) = oldToNew(neigh(valid));
        res.spikeNeigh = uint32(out);
    end

    if isfield(res, 'ordRho') && ~isempty(res.ordRho)
        ord = double(res.ordRho(:));
        ord = ord(ord > 0 & ord <= nSpikes);
        ord = ord(keepSpikes(ord));
        res.ordRho = double(oldToNew(ord));
    end

    % clusterCenters is deliberately NOT remapped here. A centre spike may
    % itself have been pruned, and there is no sensible replacement until the
    % final labelling is known, so rebuildClusterCenters does it afterwards.
end

function res = rebuildClusterCenters(res)
    %REBUILDCLUSTERCENTERS Give every surviving cluster a valid centre spike.
    %
    %   JRCLUST checks spikeClusters(clusterCenters) == 1:nClusters, so a
    %   centre must be a spike belonging to its own cluster. Pruning can remove
    %   the original centre, and leaving a 0 there makes the GUI fail with
    %   "Array indices must be positive integers" before it can even open.
    %
    %   The densest surviving spike is the natural replacement, that being what
    %   a density peak means.
    cl = double(res.spikeClusters(:));
    nClusters = max(cl);
    if nClusters < 1
        res.clusterCenters = [];
        return;
    end

    pos = find(cl > 0);
    c = cl(pos);

    if isfield(res, 'spikeRho') && numel(res.spikeRho) == numel(cl)
        rho = double(res.spikeRho(:));
        peak = accumarray(c, rho(pos), [nClusters 1], @max, -inf);
        isPeak = rho(pos) >= peak(c);
        centers = accumarray(c(isPeak), pos(isPeak), [nClusters 1], @(p) p(1), 0);
    else
        centers = accumarray(c, pos, [nClusters 1], @(p) p(1), 0);
    end

    % An empty cluster cannot happen after compaction, but never emit a 0.
    if any(centers == 0)
        error('jrc_prune:centers', ...
              '%d cluster(s) ended up with no spikes', sum(centers == 0));
    end

    res.clusterCenters = centers;
end

function res = compactClusters(res)
    %COMPACTCLUSTERS Renumber surviving units to 1..N and subset per-unit fields.
    good = unique(res.spikeClusters(res.spikeClusters > 0));
    if isempty(good), return; end

    nOld = max(good);
    keepUnits = false(nOld, 1);
    keepUnits(good) = true;

    oldToNew = zeros(nOld, 1, 'int32');
    oldToNew(good) = int32(1:numel(good));

    positive = res.spikeClusters > 0;
    res.spikeClusters(positive) = double(oldToNew(res.spikeClusters(positive)));

    vectorFields = {'clusterNotes', 'clusterSites', 'unitCount', 'unitPeaks', ...
                    'unitPeaksRaw', 'unitPeakSites', 'unitVpp', 'unitVppRaw', ...
                    'unitISIRatio', 'unitIsoDist', 'unitLRatio', 'clusterCenters', ...
                    'nSitesOverThresh', 'unitSNR'};
    for i = 1:numel(vectorFields)
        fn = vectorFields{i};
        if isfield(res, fn) && numel(res.(fn)) >= nOld
            v = res.(fn);
            res.(fn) = v(keepUnits);
        end
    end

    if isfield(res, 'clusterCentroids') && size(res.clusterCentroids, 1) >= nOld
        res.clusterCentroids = res.clusterCentroids(keepUnits, :);
    end

    waveformFields = {'meanWfGlobal', 'meanWfGlobalRaw', 'meanWfLocal', ...
                      'meanWfLocalRaw', 'meanWfRawLow', 'meanWfRawHigh'};
    for i = 1:numel(waveformFields)
        fn = waveformFields{i};
        if isfield(res, fn) && ndims(res.(fn)) == 3 && size(res.(fn), 3) >= nOld
            res.(fn) = res.(fn)(:, :, keepUnits);
        end
    end

    if isfield(res, 'waveformSim') && all(size(res.waveformSim) >= [nOld nOld])
        res.waveformSim = res.waveformSim(keepUnits, keepUnits);
    end
end

function commentInvalidZeroNSitesExcl(configFile)
    %COMMENTINVALIDZERONSITESEXCL nSitesExcl = 0 is rejected by validation.
    lines = regexp(fileread(configFile), '\r\n|\n|\r', 'split');
    for i = 1:numel(lines)
        if ~isempty(regexp(lines{i}, '^\s*nSitesExcl\s*=\s*0\s*;', 'once'))
            lines{i} = ['% ', lines{i}, ' % omitted in pruned config'];
        end
    end

    fid = fopen(configFile, 'w');
    if fid == -1, error('jrc_prune:cfg', 'cannot update %s', configFile); end
    cleaner = onCleanup(@() fclose(fid));
    fprintf(fid, '%s\n', lines{:});
    clear cleaner;
end

function appendOutputDir(configFile, outputDir)
    %APPENDOUTPUTDIR Set outputDir, REPLACING any line already there.
    %
    %   Appending is not safe. hCfg.save has already written an outputDir of
    %   its own, and it writes it as a CELL when the config holds one, which
    %   jrclust.Config then rejects outright:
    %
    %       Could not set "outputDir" ... assertion failed
    %
    %   A second line appended below it does not help, because the first one
    %   is read first. The pruned .prm then cannot be opened at all - which is
    %   how a previous run left a file carrying outputDir = {'f:\...'} on line
    %   5 and outputDir = 'F:\...\pruned' on line 162.
    %
    %   So: replace the first outputDir line, comment out any others, and only
    %   append when the file has none.
    lines = readFileLines(configFile);

    pattern = '^\s*outputDir\s*=';
    hits = find(~cellfun(@isempty, regexp(lines, pattern, 'once')));
    newLine = sprintf('outputDir = ''%s''; %% Pruned copy output folder', ...
                      strrep(outputDir, '''', ''''''));

    if isempty(hits)
        lines{end+1} = '';
        lines{end+1} = newLine;
    else
        lines{hits(1)} = newLine;
        for k = hits(2:end)
            lines{k} = ['% ', lines{k}, ' % superseded by the line above'];
        end
    end

    fid = fopen(configFile, 'w');
    if fid == -1, error('jrc_prune:cfg', 'cannot update %s', configFile); end
    cleaner = onCleanup(@() fclose(fid));
    fprintf(fid, '%s\n', lines{:});
    clear cleaner;
end

function lines = readFileLines(configFile)
    fid = fopen(configFile, 'r');
    if fid == -1, error('jrc_prune:cfg', 'cannot read %s', configFile); end
    cleaner = onCleanup(@() fclose(fid));
    raw = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
    clear cleaner;
    lines = raw{1};
end
