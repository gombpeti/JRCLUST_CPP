function jrc_mark_noise(prmFile, varargin)
%JRC_MARK_NOISE Mark noise units for deletion, without opening the GUI.
%
%   jrc_mark_noise('mydata.prm')
%       Converts every unit annotated 'noise' into 'to_delete', which is the
%       annotation JRCLUST's Edit > Delete annotated (and jrc_prune) acts on.
%
%   jrc_mark_noise('mydata.prm', 'units', [3 17 42])
%       Additionally marks those unit IDs, whatever they are annotated now.
%
%   jrc_mark_noise('mydata.prm', 'from', {'noise','multi'})
%       Treats every listed annotation as noise.
%
%   jrc_mark_noise('mydata.prm', 'undo', true)
%       Puts 'to_delete' units back to 'noise', in case you change your mind.
%
%   Options:
%     'outputDir'  where the _res.mat lives, if the .prm's outputDir is stale
%     'dryRun'     report what would change without writing (default false)
%
%   This edits the annotations in the existing _res.mat in place. It does not
%   remove anything - use jrc_prune to actually drop the spikes.

    opts = struct('units', [], 'from', {{'noise'}}, 'undo', false, ...
                  'outputDir', '', 'dryRun', false);
    for k = 1:2:numel(varargin)
        name = varargin{k};
        if ~isfield(opts, name)
            error('jrc_mark_noise:badOption', 'unknown option ''%s''', name);
        end
        opts.(name) = varargin{k+1};
    end
    if ischar(opts.from), opts.from = {opts.from}; end

    resFile = jrc_res_path(prmFile, opts.outputDir);
    fprintf('Loading %s\n', resFile);
    res = load(resFile);

    if ~isfield(res, 'clusterNotes') || isempty(res.clusterNotes)
        nClusters = double(max(res.spikeClusters));
        res.clusterNotes = repmat({''}, nClusters, 1);
        fprintf('  no annotations present; created %d empty ones\n', nClusters);
    end

    notes = res.clusterNotes(:);
    nClusters = numel(notes);

    before = sum(strcmp(notes, 'to_delete'));

    if opts.undo
        hit = strcmp(notes, 'to_delete');
        notes(hit) = {'noise'};
        fprintf('  reverted %d units from to_delete back to noise\n', sum(hit));
    else
        % Annotations that count as noise.
        hit = false(nClusters, 1);
        for i = 1:numel(opts.from)
            match = strcmp(notes, opts.from{i});
            fprintf('  %-10s : %d units\n', opts.from{i}, sum(match));
            hit = hit | match;
        end

        % Explicit unit IDs, regardless of their current annotation.
        if ~isempty(opts.units)
            bad = opts.units(opts.units < 1 | opts.units > nClusters);
            if ~isempty(bad)
                error('jrc_mark_noise:badUnit', ...
                      'unit id(s) out of range 1..%d: %s', nClusters, mat2str(bad));
            end
            hit(opts.units) = true;
            fprintf('  explicit  : %d units\n', numel(opts.units));
        end

        notes(hit) = {'to_delete'};
    end

    after = sum(strcmp(notes, 'to_delete'));

    % How many spikes this actually affects, so the cost is visible up front.
    marked = find(strcmp(notes, 'to_delete'));
    nSpikesAffected = sum(ismember(double(res.spikeClusters), marked));

    fprintf('\n  to_delete: %d -> %d units (%d spikes, %.1f%% of %d)\n', ...
            before, after, nSpikesAffected, ...
            100 * nSpikesAffected / numel(res.spikeClusters), ...
            numel(res.spikeClusters));

    if opts.dryRun
        fprintf('\ndry run: nothing written\n');
        return;
    end

    res.clusterNotes = notes;
    fprintf('\nWriting %s\n', resFile);
    jrclust.utils.saveStruct(res, resFile);

    fprintf('Done. Next: jrc_prune(''%s'')\n', prmFile);
end
