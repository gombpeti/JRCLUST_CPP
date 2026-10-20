function jrc_recompute(prmFile, varargin)
%JRC_RECOMPUTE Rebuild the derived per-cluster data the GUI needs.
%
%   Recomputes mean waveforms, waveform similarity, centroids and quality
%   scores from the current spike table, using JRCLUST's own methods.
%
%   Needed after anything that moves spikes between clusters outside the GUI
%   - jrc_split_noise, jrc_automerge - because those change what every
%   derived figure describes, and the GUI reads them from the file rather
%   than regenerating them on load. waveformSim matters most: it is what the
%   space key uses to jump to the most similar unit, so a blank or stale one
%   silently breaks navigation and the similarity matrix shows nothing.
%
%   jrc_recompute(prmFile)
%   jrc_recompute(prmFile, 'outputDir', 'F:\...\resort')
%
%   Needs the _filt.jrc and _raw.jrc files beside the result; they are loaded
%   for the computation and dropped again before saving.

    opts = struct('outputDir', '');
    for k = 1:2:numel(varargin)
        name = varargin{k};
        if ~isfield(opts, name)
            error('jrc_recompute:badOption', 'unknown option ''%s''', name);
        end
        opts.(name) = varargin{k+1};
    end

    resFile = jrc_res_path(prmFile, opts.outputDir);
    fprintf('Loading %s\n', resFile);
    res = load(resFile);
    nClusters = max(double(res.spikeClusters(:)));
    fprintf('  %d spikes, %d clusters\n', numel(res.spikeClusters), nClusters);

    if isfield(res, 'waveformSim') && ~isempty(res.waveformSim)
        before = sprintf('%s, %.0f%% finite', mat2str(size(res.waveformSim)), ...
                         100*mean(isfinite(res.waveformSim(:))));
    else
        before = 'empty';
    end
    fprintf('  waveformSim before: %s\n', before);

    res = recomputeClustering(res, prmFile, resFile);

    fprintf('\nWriting %s\n', resFile);
    jrclust.utils.saveStruct(res, resFile);
    fprintf('\nDone. Reopen the session; the similarity matrix and the space\n');
    fprintf('key should work again.\n');
end

% ---------------------------------------------------------------------------
function res = recomputeClustering(res, prmFile, resFile)
    %RECOMPUTECLUSTERING Let JRCLUST regenerate the derived per-cluster data.
    %
    %   Mean waveforms, waveform similarity, centroids and quality scores all
    %   describe a clustering that no longer exists once spikes have moved.
    %   The GUI reads them from the file and does not rebuild them on load -
    %   waveformSim in particular, which is what the space key uses to jump
    %   to the most similar unit.
    %
    %   The bulk arrays live in the .jrc files and have to be attached for
    %   the computation, then dropped again before saving: they are gigabytes
    %   and the GUI reloads them itself.
    fprintf('\nrecomputing waveforms, similarity and quality scores\n');
    hCfg = jrclust.Config(prmFile);

    hClust = jrclust.sort.DensityPeakClustering(hCfg);
    md = ?jrclust.sort.DensityPeakClustering;
    pl = md.PropertyList;
    fn = fieldnames(res);
    for i = 1:numel(fn)
        f = fn{i};
        p = pl(strcmp(f, {pl.Name}));
        if isprop(hClust, f) && ~(~isempty(p) && p.Dependent && isempty(p.SetMethod))
            try hClust.(f) = res.(f); catch, end %#ok<CTCH>
        end
    end

    % The spike waveforms themselves, from the .jrc files beside the result.
    prefix = regexprep(resFile, '_res\.mat$', '');
    attach = {'spikesFilt', '_filt.jrc', 'filtShape', '*int16'; ...
              'spikesRaw',  '_raw.jrc',  'rawShape',  '*int16'; ...
              'spikeFeatures', '_features.jrc', 'featuresShape', '*single'};
    for i = 1:size(attach,1)
        f = [prefix attach{i,2}];
        if ~isfile(f) || ~isfield(res, attach{i,3}), continue; end
        shape = double(res.(attach{i,3}));
        fid = fopen(f, 'r');
        c = onCleanup(@() fclose(fid));
        v = fread(fid, prod(shape), attach{i,4});
        clear c;
        try hClust.(attach{i,1}) = reshape(v, shape); catch, end %#ok<CTCH>
    end

    hClust.updateWaveforms();      % mean waveforms AND waveformSim
    try hClust.computeCentroids(); catch, end %#ok<CTCH>
    try hClust.computeQualityScores(); catch, end %#ok<CTCH>

    % Flatten back, as @JRC/saveRes.m does.
    props = fieldnames(hClust);
    for i = 1:numel(props)
        f = props{i};
        try res.(f) = hClust.(f); catch, end %#ok<CTCH>
    end

    bulk = {'spikesRaw','spikesRawVolt','spikesFilt','spikesFiltVolt', ...
            'spikesFilt2','spikesFilt3','spikeFeatures','nClusters', ...
            'nEdits','unitFields','hRecs','dRes','sRes','hCfg'};
    for i = 1:numel(bulk)
        if isfield(res, bulk{i}), res = rmfield(res, bulk{i}); end
    end

    if isfield(res,'waveformSim') && ~isempty(res.waveformSim)
        fprintf('  waveformSim: %s, %.0f%% finite\n', mat2str(size(res.waveformSim)), ...
                100*mean(isfinite(res.waveformSim(:))));
    else
        warning('jrc_split_noise:noSim', ...
                'waveformSim is still empty - the space key will not work');
    end
end

