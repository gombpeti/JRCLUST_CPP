function jrc_import_cpp(prmFile, varargin)
%JRC_IMPORT_CPP Import results from the C++/CUDA sorter into JRCLUST.
%
%   jrc_import_cpp('mydata.prm')
%   jrc_import_cpp('mydata.prm', 'autoMerge', false)
%   jrc_import_cpp('mydata.prm', 'rmOutliers', true, 'rmRefrac', true)
%
%   Options (name/value):
%     'autoMerge'  (default true)  merge similar units by mean-waveform
%                                  correlation, the same step
%                                  jrc detect-sort runs. Uses maxUnitSim,
%                                  autoMergeBy and nPassesMerge from the .prm.
%     'rmOutliers' (default false) Mahalanobis outlier removal per unit
%                                  (outlierThresh in the .prm).
%     'rmRefrac'   (default false) drop spikes violating the refractory
%                                  period within each unit (refracInt).
%     'outputDir'  (default '')    where the .jrc files are, and where
%                                  _res.mat is written. Overrides outputDir
%                                  from the .prm, which is useful when that
%                                  path is stale.
%
%   Reads the files written by `jrclust.exe sort <prm> --save`:
%       <session>_cpp.bin       spike table and clustering
%       <session>_filt.jrc      filtered waveforms  (int16)
%       <session>_raw.jrc       raw waveforms       (int16)
%       <session>_features.jrc  PCA features        (single)
%
%   and writes <session>_res.mat in JRCLUST's own format, so that
%
%       jrc manual mydata.prm
%
%   opens it for manual curation exactly as if it had been produced by
%   `jrc detect-sort`.
%
%   The derived quantities (mean waveforms, SNR, waveform similarity, quality
%   scores) are NOT taken from the C++ side. They are computed here by
%   JRCLUST's own DensityPeakClustering constructor, so they are guaranteed
%   consistent with the rest of the toolbox.

    if nargin < 1
        error('jrc_import_cpp:usage', 'Usage: jrc_import_cpp(''mydata.prm'')');
    end
    if ~exist(prmFile, 'file')
        error('jrc_import_cpp:noPrm', '%s does not exist', prmFile);
    end

    opts = struct('autoMerge', true, 'rmOutliers', false, 'rmRefrac', false, ...
                  'outputDir', '');
    for k = 1:2:numel(varargin)
        name = varargin{k};
        if ~isfield(opts, name)
            error('jrc_import_cpp:badOption', 'unknown option ''%s''', name);
        end
        if strcmp(name, 'outputDir')
            opts.outputDir = varargin{k+1};
        else
            opts.(name) = logical(varargin{k+1});
        end
    end

    hCfg = jrclust.Config(prmFile);

    [outDir, session] = fileparts(hCfg.configFile);
    if ~isempty(hCfg.outputDir)
        outDir = hCfg.outputDir;
    end
    if ~isempty(opts.outputDir)
        outDir = opts.outputDir;
    end
    if ~exist(outDir, 'dir')
        error('jrc_import_cpp:noOutputDir', ...
              ['output directory %s does not exist.\n' ...
               'Pass ''outputDir'' if the .prm''s outputDir is stale.'], outDir);
    end
    prefix = fullfile(outDir, session);

    % _res.mat must sit beside the .jrc files, or loadFiles will not find them.
    resFile = [prefix '_res.mat'];

    sidecarFile = [prefix '_cpp.bin'];
    if ~exist(sidecarFile, 'file')
        error('jrc_import_cpp:noSidecar', ...
              ['%s not found.\n' ...
               'Run:  jrclust.exe sort "%s" --save'], sidecarFile, prmFile);
    end

    fprintf('Importing C++ results for %s\n', session);

    % ---------------------------------------------------------------- sidecar
    S = readSidecar(sidecarFile);

    fprintf('  spikes              %d\n', S.nSpikes);
    fprintf('  clusters            %d\n', S.nClusters);
    fprintf('  features            %d x %d per spike\n', S.nFeatures, S.nPeaks);

    % Sanity-check the probe against what the sorter used, since a mismatched
    % .prm would silently produce nonsense downstream.
    if S.nSites ~= hCfg.nSites
        error('jrc_import_cpp:siteMismatch', ...
              'sorter used %d sites but this .prm has %d', S.nSites, hCfg.nSites);
    end

    % ------------------------------------------------------------- waveforms
    filtShape = [S.filtSamples, S.nSiteNeighbors, S.nSpikes];
    rawShape  = [S.rawSamples,  S.nSiteNeighbors, S.nSpikes];
    featShape = [S.nFeatures,   S.nPeaks,         S.nSpikes];

    spikesFilt    = readBin([prefix '_filt.jrc'],     filtShape, '*int16');
    spikesRaw     = readBin([prefix '_raw.jrc'],      rawShape,  '*int16');
    spikeFeatures = readBin([prefix '_features.jrc'], featShape, '*single');

    % ------------------------------------------------------------------ dRes
    nSites = hCfg.nSites;

    dRes = struct();
    dRes.spikeTimes  = S.spikeTimes;
    dRes.spikeAmps   = S.spikeAmps;
    dRes.spikeSites  = S.spikeSites;
    dRes.spikeSites2 = S.spikeSites2;
    dRes.spikeSites3 = [];

    % ORIENTATION MATTERS. These must be 1 x nSites ROW cell arrays holding
    % COLUMN vectors, which is what detect.m produces. getClusterFeaturesSite
    % does neCell2mat(spikesBySite(neighbors)') with `neighbors` a column
    % vector: indexing a row cell array with it yields a row, the transpose
    % makes it a column, and cell2mat then stacks the column vectors along
    % dim 1. Store these as nSites x 1 instead and the transpose produces a
    % row, so cell2mat tries to concatenate unequal-length columns along dim 2
    % and errors with "Dimensions of arrays being concatenated are not
    % consistent" the first time a cluster is clicked in the GUI.
    dRes.spikesBySite  = arrayfun(@(s) find(S.spikeSites  == s), 1:nSites, ...
                                  'UniformOutput', 0);
    if any(S.spikeSites2 > 0)
        dRes.spikesBySite2 = arrayfun(@(s) find(S.spikeSites2 == s), 1:nSites, ...
                                      'UniformOutput', 0);
    else
        dRes.spikesBySite2 = cell(1, nSites);
    end
    dRes.spikesBySite3 = cell(1, nSites);

    dRes.centerSites = [S.spikeSites, S.spikeSites2];

    dRes.spikesRaw     = spikesRaw;
    dRes.spikesFilt    = spikesFilt;
    dRes.spikeFeatures = spikeFeatures;
    dRes.spikesFilt2   = zeros(S.filtSamples, S.nSiteNeighbors, 0, 'int16');
    dRes.spikesFilt3   = zeros(S.filtSamples, S.nSiteNeighbors, 0, 'int16');

    dRes.rawShape      = rawShape;
    dRes.filtShape     = filtShape;
    dRes.featuresShape = featShape;

    dRes.siteThresh     = S.siteThresh;
    dRes.meanSiteThresh = S.siteThresh;
    dRes.siteRMS        = S.siteThresh / max(hCfg.qqFactor, eps);

    dRes.spikePositions = spikePositions(S, hCfg);

    dRes.nSpikes    = S.nSpikes;
    dRes.detectTime = S.detectSeconds;
    dRes.detectedOn = now();

    % ------------------------------------------------------------------ sRes
    sRes = struct();
    sRes.spikeRho   = S.spikeRho;
    sRes.spikeDelta = S.spikeDelta;
    sRes.spikeNeigh = S.spikeNeigh;

    % ordRho: spikes ordered by descending density. JRCLUST propagates cluster
    % labels along this ordering.
    [~, ordRho] = sort(S.spikeRho, 'descend');
    sRes.ordRho = ordRho;

    sRes.spikeClusters     = S.spikeClusters;
    sRes.initialClustering = S.spikeClusters;
    sRes.clusterCenters = S.clusterCenters;
    sRes.clusterSites   = double(S.clusterSites);
    sRes.unitCount      = double(S.unitCount);
    sRes.spikesByCluster = arrayfun(@(c) find(S.spikeClusters == c), ...
                                    1:S.nClusters, 'UniformOutput', 0);

    sRes.rhoCutSite   = S.rhoCutSite;
    sRes.rhoCuts      = S.rhoCutSite;
    sRes.rhoCutGlobal = [];

    sRes.sortTime = S.sortSeconds;
    sRes.sortedOn = now();

    % ---------------------------------------------------- build and commit
    % The constructor runs commit(), which fills in mean waveforms, waveform
    % similarity and the quality scores the GUI displays.
    fprintf('  building clustering object (computing mean waveforms)...\n');
    hClust = jrclust.sort.DensityPeakClustering(hCfg, sRes, dRes);

    nBefore = hClust.nClusters;
    fprintf('  clusters from the C++ sorter: %d\n', nBefore);

    % ------------------------------------------------- post-processing
    % These are JRCLUST's own routines, not reimplementations, so the result
    % matches what jrc detect-sort would have produced from the same spikes.
    if opts.rmRefrac
        fprintf('  removing refractory-period violations...\n');
        total = 0;
        for iC = 1:hClust.nClusters
            total = total + hClust.rmRefracSpikes(iC);
        end
        fprintf('    removed %d spikes\n', total);
    end

    if opts.rmOutliers
        fprintf('  removing outlier spikes (outlierThresh = %g MAD)...\n', ...
                hCfg.getOr('outlierThresh', 7.5));
        hClust.rmOutlierSpikes();
    end

    if opts.autoMerge
        % SortController.m runs exactly this after building the clustering.
        fprintf('  auto-merging (%s, maxUnitSim = %g, up to %d passes)...\n', ...
                hCfg.getOr('autoMergeBy', 'pearson'), ...
                hCfg.getOr('maxUnitSim', 0.98), ...
                hCfg.getOr('nPassesMerge', 10));

        if ~hClust.autoMerge()
            warning('jrc_import_cpp:autoMerge', 'auto-merge did not complete');
        end

        fprintf('    %d -> %d clusters (%d merged away)\n', ...
                nBefore, hClust.nClusters, nBefore - hClust.nClusters);
    end

    % ------------------------------------------------- unit labels from C++
    % noise / single / multi, one per cluster. These are indexed by the
    % clustering the sorter produced, so re-merging here would misalign them.
    if isfield(S, 'clusterNotes') && ~isempty(S.clusterNotes)
        if opts.autoMerge
            warning('jrc_import_cpp:staleNotes', ...
                    ['auto-merge renumbered the clusters, so the sorter''s ' ...
                     'unit labels no longer line up and were dropped. Re-run ' ...
                     'with ''autoMerge'', false to keep them.']);
        elseif numel(S.clusterNotes) ~= hClust.nClusters
            warning('jrc_import_cpp:noteCount', ...
                    'got %d labels for %d clusters; skipping', ...
                    numel(S.clusterNotes), hClust.nClusters);
        else
            applied = 0;
            for iC = 1:hClust.nClusters
                note = S.clusterNotes{iC};
                if ~isempty(note)
                    hClust.addNote(iC, note);
                    applied = applied + 1;
                end
            end

            counts = @(t) sum(strcmp(S.clusterNotes, t));
            fprintf('  unit labels applied: %d  (noise %d, single %d, multi %d, to_delete %d)\n', ...
                    applied, counts('noise'), counts('single'), counts('multi'), ...
                    counts('to_delete'));
            if counts('to_delete') > 0
                fprintf('    -> in the GUI: Edit > Delete annotated removes them\n');
            elseif counts('noise') > 0
                fprintf('    -> in the GUI: View > Display annotated units only to review\n');
            end
        end
    end

    % ----------------------------------------------------------------- save
    % Mirrors @JRC/saveRes.m: flatten hClust's properties into the top level,
    % then drop the bulk arrays that live in the .jrc files. JRC.res is
    % read-only, so we write the file directly instead of going through it.
    res = dRes;
    fn = fieldnames(sRes);
    for i = 1:numel(fn)
        res.(fn{i}) = sRes.(fn{i});
    end

    res.initialClustering = hClust.initialClustering;
    res.spikeClusters     = hClust.spikeClusters;

    hClustFields = fieldnames(hClust);
    for i = 1:numel(hClustFields)
        fname = hClustFields{i};
        try
            res.(fname) = hClust.(fname);
        catch
            % skip properties that error on read (transient/dependent)
        end
    end

    % These are reconstructed from the .jrc files on load, so keep them out of
    % the .mat - they are gigabytes.
    bulk = {'spikesRaw', 'spikesRawVolt', 'spikesFilt', 'spikesFiltVolt', ...
            'spikesFilt2', 'spikesFilt3', 'spikeFeatures', 'nClusters', ...
            'nEdits', 'unitFields', 'hRecs', 'dRes', 'sRes', 'hCfg'};
    for i = 1:numel(bulk)
        if isfield(res, bulk{i})
            res = rmfield(res, bulk{i});
        end
    end

    fprintf('  writing %s\n', resFile);
    jrclust.utils.saveStruct(res, resFile);

    fprintf('\nDone. Now run:\n    jrc manual %s\n', prmFile);
end

%% ------------------------------------------------------------------ helpers

function S = readSidecar(filename)
    %READSIDECAR Read the C++ sorter's spike table.
    fid = fopen(filename, 'r');
    if fid < 0
        error('jrc_import_cpp:open', 'cannot open %s', filename);
    end

    try
        magic = fread(fid, 8, '*char')';
        if ~strcmp(magic, 'JRCPP004')
            error('jrc_import_cpp:magic', ...
                  'unexpected file signature ''%s'' in %s (rebuild the sorter?)', ...
                  magic, filename);
        end

        S.nSpikes = fread(fid, 1, 'int64');

        hdr = fread(fid, 7, 'int32');
        S.nSites         = hdr(1);
        S.nClusters      = hdr(2);
        S.nFeatures      = hdr(3);
        S.nPeaks         = hdr(4);
        S.filtSamples    = hdr(5);
        S.rawSamples     = hdr(6);
        S.nSiteNeighbors = hdr(7);

        reals = fread(fid, 3, 'double');
        S.sampleRate    = reals(1);
        S.detectSeconds = reals(2);
        S.sortSeconds   = reals(3);

        n = S.nSpikes;
        c = S.nClusters;
        s = S.nSites;

        S.spikeTimes    = fread(fid, n, '*int32');
        S.spikeAmps     = fread(fid, n, '*int16');
        S.spikeSites    = fread(fid, n, '*int32');
        S.spikeSites2   = fread(fid, n, '*int32');
        S.spikeRho      = fread(fid, n, '*single');
        S.spikeDelta    = fread(fid, n, '*single');
        S.spikeNeigh    = fread(fid, n, '*int32');
        S.spikeClusters = fread(fid, n, '*int32');

        S.clusterCenters = fread(fid, c, '*int32');
        S.clusterSites   = fread(fid, c, '*int32');
        S.unitCount      = fread(fid, c, '*int32');

        S.siteThresh  = fread(fid, s, '*single');
        S.rhoCutSite  = fread(fid, s, '*single');

        % Feature-weighted centroids, [x,y] interleaved -> nSpikes x 2.
        pos = fread(fid, n * 2, '*single');
        S.spikePositions = reshape(pos, 2, [])';

        % Per-cluster annotations, length-prefixed.
        S.clusterNotes = cell(c, 1);
        for iC = 1:c
            len = fread(fid, 1, 'int32');
            if isempty(len)
                S.clusterNotes{iC} = '';
            elseif len > 0
                S.clusterNotes{iC} = fread(fid, len, '*char')';
            else
                S.clusterNotes{iC} = '';
            end
        end
    catch ME
        fclose(fid);
        rethrow(ME);
    end

    fclose(fid);

    % JRCLUST expects uint32 neighbour indices and column vectors throughout.
    S.spikeNeigh     = uint32(S.spikeNeigh);
    S.clusterCenters = double(S.clusterCenters);
end

function data = readBin(filename, shape, dataType)
    %READBIN Load a .jrc binary and reshape it, or return empty if absent.
    if ~exist(filename, 'file')
        warning('jrc_import_cpp:missing', '%s not found; leaving empty', filename);
        data = [];
        return;
    end

    fid = fopen(filename, 'r');
    data = fread(fid, Inf, dataType);
    fclose(fid);

    expected = prod(shape);
    if numel(data) ~= expected
        error('jrc_import_cpp:shape', ...
              '%s holds %d values but shape %s needs %d', ...
              filename, numel(data), mat2str(shape), expected);
    end

    data = reshape(data, shape);
end

function positions = spikePositions(S, hCfg)
    %SPIKEPOSITIONS Feature-weighted centroid of each spike.
    %
    %   Taken from the sorter, which computes it exactly as spikePos.m does.
    %   These drive drift-aware mean waveforms, so using real centroids rather
    %   than a site-location stand-in matters: with one y per site the drift
    %   bands collapse and merging behaves differently.
    %
    %   Orientation matters: getClusterMean indexes this as (spikes, 2), so it
    %   must be nSpikes x 2, NOT 2 x nSpikes.
    if isfield(S, 'spikePositions') && ~isempty(S.spikePositions) && ...
            any(S.spikePositions(:) ~= 0)
        positions = single(S.spikePositions);
    else
        warning('jrc_import_cpp:noPositions', ...
                'sorter supplied no spike positions; falling back to site locations');
        positions = single(hCfg.siteLoc(S.spikeSites, :));
    end
end
