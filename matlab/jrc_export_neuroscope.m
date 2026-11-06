function groups = jrc_export_neuroscope(prmFile, varargin)
%JRC_EXPORT_NEUROSCOPE Write a JRCLUST result as Neuroscope .res/.clu files.
%
%   Splits the sort back into its recordings and, next to each recording's
%   Neuroscope .xml, writes a .res.N / .clu.N pair per spike group N (and,
%   with 'klusters', the .spk.N / .fet.N Klusters needs). It rewrites that
%   .xml's spike groups so Neuroscope draws every unit's spikes, in the
%   unit's own colour, on the channels around it.
%
%   jrc_export_neuroscope(prm)                  every recording in the .prm
%   jrc_export_neuroscope                       asks for the .prm first
%   jrc_export_neuroscope(prm, 'dryRun', true)  prints the groups, writes nothing
%
%   Everything comes from the .prm: each entry of rawRecordings is one
%   session, and that session's files go to its own folder, beside its .dat,
%   named after it - <folder>\<dat name>.xml/.nrs/.res.N/.clu.N/.spk.N/.fet.N.
%
%   Options
%     'grouping'     'radius' (default) | 'jrclust' | 'shank' | a number of
%                    sites per block - see SPIKE GROUPS
%     'radiusUm'     for 'radius': a unit is drawn on every channel of its
%                    shank at most this far above or below its peak site,
%                    by siteLoc in the .prm (default 100)
%     'minSites'     smallest group for 'jrclust' grouping (default 8)
%     'maxSites'     largest group for 'jrclust' grouping (default 24)
%     'xmlFiles'     cell array, one .xml per entry of rawRecordings, to use
%                    instead of the one beside each .dat
%     'ask'          true: always prompt for each .xml; false: never prompt.
%                    By default it prompts only when the .xml is missing.
%     'nrs'          register the .clu files in each session's .nrs, so
%                    Neuroscope opens them with the .dat (default true)
%     'useXmlSkip'   also keep units off the channels each session's .xml
%                    marks skip="1" (default true). The .prm's ignoreSites
%                    are always kept clear: they were never sorted.
%     'selectUnits'  'all' (default): every unit starts selected, so its
%                    spikes are coloured as soon as the .dat opens;
%                    'none': tick them in the cluster palette yourself
%     'raster'       also keep Neuroscope's raster below the traces
%                    (default false; with every unit selected it has a row
%                    per unit and channel file)
%     'skipNotes'    units annotated with any of these are left out,
%                    e.g. {'to_delete'} (default {}: export every unit)
%     'includeNoise' also export spikes JRCLUST assigned to no unit, as
%                    cluster 0 (default false)
%     'klusters'     also write .spk.N and .fet.N, the waveform and feature
%                    files only Klusters needs (default false). Reads every
%                    .dat once, front to back - minutes per session, where
%                    Neuroscope's .res/.clu/.xml/.nrs take seconds.
%     'highPass'     cutoff in Hz of the filter applied to .spk waveforms
%                    (default 300)
%     'spikeColors'  colour each spike group's channels in Neuroscope's
%                    spike-group palette (default true)
%     'outputDir'    where the _res.mat is, if the .prm's outputDir is stale
%     'dryRun'       compute and report, write nothing (default false)
%     'checkNeuroscope'  refuse to run while Neuroscope is open (default
%                    true): it writes its own .xml/.nrs back when it closes
%
%   CHANNELS. JRCLUST numbers sites 1..nSites - spikeSites, clusterSites,
%   siteNeighbors, shankMap and siteLoc all index sites. The recording stores
%   site s as channel siteMap(s), counted from 1. Neuroscope counts channels
%   from 0, so every site is written to the .xml as siteMap(s) - 1, the same
%   numbering as the existing anatomical groups.
%
%   SPIKE GROUPS. Neuroscope gives each channel to exactly one spike group,
%   and draws a unit's spikes on a channel only if the .clu file of that
%   channel's group holds the unit.
%
%   'radius' (default) makes every channel its own spike group and writes
%   each unit into the .clu.N of every channel of its shank within radiusUm
%   above or below its peak site - the site where its amplitude is largest
%   (clusterSites). How many channels that is follows from the probe
%   geometry in the .prm: 2-3 at ~100 um pitch, ~13 at 15.5 um. Every unit
%   is centred on its peak, and nothing is drawn on channels too far away
%   to see it. A unit therefore appears in several .clu files - once per
%   channel - and is ticked per channel in the cluster palette.
%
%   'jrclust' cuts each shank into contiguous blocks of neighbouring sites
%   instead, placed where the fewest units' mean waveforms spread across a
%   boundary, and puts each unit in the one block holding its peak site.
%
%   CLUSTER IDS are JRCLUST's unit numbers, unchanged, so unit 57 in
%   Neuroscope is unit 57 in jrc manual, in every .clu file it appears in.
%
%   TIMES. A concatenated sort numbers samples continuously across all
%   rawRecordings, from 1. Each .res holds 0-based sample numbers within its
%   own recording, which is what Neuroscope reads against that .dat.
%
%   The .xml and .nrs are backed up next to themselves before they change.
%   In the .xml only <spikeDetection>, <units> and, with 'spikeColors', the
%   spikeColor of each channel are rewritten; the anatomical groups are left
%   as they are. In the .nrs the cluster-file entries of this export are
%   replaced, each unit in its jrc manual colour; every display gets
%   waveform display (Units > Waveforms) and, with 'selectUnits' 'all',
%   every unit selected. Event files and the rest of the display settings
%   are kept. Opening the .dat in Neuroscope then shows each unit's spikes
%   coloured on the channels around its peak.

    if nargin < 1 || isempty(prmFile)
        [f, p] = uigetfile('*.prm', 'JRCLUST parameter file');
        if isequal(f, 0), return; end
        prmFile = fullfile(p, f);
    end

    opts = struct('grouping', 'radius', 'radiusUm', 100, 'minSites', 8, 'maxSites', 24, ...
                  'xmlFiles', {{}}, 'ask', [], 'nrs', true, ...
                  'selectUnits', 'all', 'raster', false, 'useXmlSkip', true, ...
                  'skipNotes', {{}}, 'includeNoise', false, ...
                  'klusters', false, 'highPass', 300, ...
                  'spikeColors', true, 'outputDir', '', 'dryRun', false, ...
                  'checkNeuroscope', true);
    for k = 1:2:numel(varargin)
        name = varargin{k};
        if ~isfield(opts, name)
            error('jrc_export_neuroscope:badOption', 'unknown option ''%s''', name);
        end
        opts.(name) = varargin{k+1};
    end
    if ischar(opts.skipNotes), opts.skipNotes = {opts.skipNotes}; end

    % Neuroscope keeps a session's .xml and .nrs in memory and writes both
    % back when it saves or closes, so a session open during the export
    % silently loses the new spike groups and display settings afterwards.
    if ~opts.dryRun && opts.checkNeuroscope && neuroscopeRunning()
        error('jrc_export_neuroscope:neuroscopeOpen', ...
              ['Neuroscope is running. Close every Neuroscope window first: on closing it ' ...
               'writes its own copy of the .xml and .nrs over the export. ' ...
               '(''checkNeuroscope'', false skips this check.)']);
    end

    % ------------------------------------------------------------ load --
    hCfg = jrclust.Config(prmFile);
    resFile = jrc_res_path(prmFile, opts.outputDir);
    fprintf('Loading %s\n', resFile);
    want = {'spikeTimes', 'spikeClusters', 'spikeSites', 'clusterSites', ...
            'clusterNotes', 'meanWfGlobal', 'unitIsoDist', 'unitSNR'};
    have = who('-file', resFile);
    res = load(resFile, want{ismember(want, have)});

    siteMap = double(hCfg.siteMap(:));
    ignoredSites = double(hCfg.ignoreSites(:));          % site indices, never sorted
    shank = double(hCfg.shankMap(:));
    siteLoc = double(hCfg.siteLoc);
    nClusters = double(max(res.spikeClusters));
    clusterSites = double(res.clusterSites(:));
    notes = repmat({''}, nClusters, 1);
    if isfield(res, 'clusterNotes'), notes(1:numel(res.clusterNotes)) = res.clusterNotes(:); end
    notes(cellfun(@(n) ~ischar(n) && ~isstring(n), notes)) = {''};

    spikeClusters = double(res.spikeClusters(:));
    spikeTimes = double(res.spikeTimes(:));
    nPerUnit = accumarray(spikeClusters(spikeClusters > 0), 1, [nClusters 1]);

    units = find(nPerUnit > 0 & ~ismember(notes, opts.skipNotes));
    fprintf('  %d units, %d spikes (%d not in any unit)\n', numel(units), ...
            sum(nPerUnit(units)), sum(spikeClusters <= 0));
    if numel(units) < nClusters
        fprintf('  leaving out %d units (empty or annotated %s)\n', ...
                nClusters - numel(units), strjoin(cellstr(opts.skipNotes), '/'));
    end

    % ---------------------------------------------------- spike groups --
    [siteGroup, mode] = buildGroups(opts, res, hCfg, shank, siteLoc, clusterSites, units);
    nGroups = max(siteGroup);
    groupSites = arrayfun(@(g) orderSites(find(siteGroup == g), siteLoc), 1:nGroups, ...
                          'UniformOutput', false);

    % unitGroup: the group of a unit's peak site. unitGroups: every group the
    % unit is drawn on - with 'radius', the group of each site of its shank
    % within radiusUm above or below its peak site; otherwise just its own.
    unitGroup = zeros(nClusters, 1);
    unitGroup(units) = siteGroup(clusterSites(units));
    unitGroups = cell(nClusters, 1);
    for u = units(:)'
        p = clusterSites(u);
        if strcmp(mode, 'radius')
            near = shank == shank(p) & abs(siteLoc(:, 2) - siteLoc(p, 2)) <= opts.radiusUm;
            unitGroups{u} = unique(siteGroup(near));
        else
            unitGroups{u} = unitGroup(u);
        end
    end

    reportGroups(mode, opts, groupSites, shank, siteLoc, siteMap, unitGroups, units, res, siteGroup);
    if strcmp(mode, 'radius')
        colorKey = num2cell(cellfun(@(s) shank(s(1)), groupSites));   % colour by shank
    else
        colorKey = num2cell(1:nGroups);
    end
    groups = struct('sites', groupSites, ...
                    'channels', cellfun(@(s) siteMap(s)' - 1, groupSites, 'UniformOutput', false), ...
                    'units', arrayfun(@(g) units(cellfun(@(c) any(c == g), unitGroups(units)))', ...
                                      1:nGroups, 'UniformOutput', false), ...
                    'colorKey', colorKey);

    % ------------------------------------------------ recordings & times --
    recs = cellstr(hCfg.rawRecordings);
    nRecs = numel(recs);
    bytesPerSample = numel(typecast(cast(0, hCfg.dataType), 'uint8'));
    nSamples = zeros(nRecs, 1);
    for r = 1:nRecs
        d = dir(recs{r});
        if isempty(d)
            error('jrc_export_neuroscope:noDat', ...
                  '%s not found - the .prm''s rawRecordings must point at the recordings', recs{r});
        end
        nSamples(r) = (d.bytes - hCfg.headerOffset) / (hCfg.nChans * bytesPerSample);
    end
    offsets = [0; cumsum(nSamples)];
    if max(spikeTimes) > offsets(end)
        error('jrc_export_neuroscope:timesBeyondEnd', ...
              ['spike time %d is beyond the %d samples of rawRecordings - was the ' ...
               'sort run on a different set or order of files?'], max(spikeTimes), offsets(end));
    end
    spikeRec = discretize(spikeTimes - 1, offsets);         % times are 1-based

    % which spikes go out, and with what cluster id
    out = ismember(spikeClusters, units);
    clu = spikeClusters;
    if opts.includeNoise
        noise = spikeClusters <= 0;
        clu(noise) = 0;
        out = out | noise;
    end

    xmlFiles = resolveXmlFiles(opts, recs);

    for r = 1:nRecs
        if isempty(xmlFiles{r})
            fprintf('\n[%d/%d] %s: no .xml chosen, skipped\n', r, nRecs, recs{r});
            continue;
        end
        [xmlDir, base] = fileparts(xmlFiles{r});
        inRec = out & spikeRec == r;
        fprintf('\n[%d/%d] %s\n  %d spikes, samples %d..%d of the concatenation\n', ...
                r, nRecs, xmlFiles{r}, sum(inRec), offsets(r) + 1, offsets(r + 1));
        if ~isfile(fullfile(xmlDir, [base '.dat']))
            warning('jrc_export_neuroscope:noDatBesideXml', ...
                    ['no %s.dat beside %s - Neuroscope only uses an .xml with the ' ...
                     'same name and folder as the .dat it opens'], base, xmlFiles{r});
        end
        checkXml(xmlFiles{r}, hCfg, siteMap);

        % Channels nothing is drawn on in this session: the .prm's ignored
        % sites, which were never sorted, and the channels this session's
        % .xml marks skip="1" (dead or noisy in this recording). They leave
        % every spike group, so no unit's window reaches them.
        excluded = ignoredSites;
        nFromXml = 0;
        if opts.useXmlSkip
            fromXml = find(ismember(siteMap - 1, xmlSkippedChannels(xmlFiles{r})));
            nFromXml = numel(setdiff(fromXml, excluded));
            excluded = union(excluded, fromXml);
        end
        [sGroups, sUnitGroups, sUnitGroup, sUnits, emptyGroups] = ...
            excludeSites(groups, unitGroups, unitGroup, units, excluded);
        lost = setdiff(units, sUnits);
        fprintf(['  %d channels left out: %d ignored in the .prm, %d more skipped in this .xml; ' ...
                 '%d unit-channel pairs dropped'], numel(excluded), numel(ignoredSites), nFromXml, ...
                sum(cellfun(@numel, unitGroups(units))) - sum(cellfun(@numel, sUnitGroups(sUnits))));
        if ~isempty(lost)
            fprintf('; %d units have no channel left and are not shown: %s', numel(lost), compactList(lost));
        end
        fprintf('\n');
        if opts.dryRun, continue; end

        % one entry per (spike, group it is drawn in): a unit's spike is
        % written into the .clu/.res of every group in sUnitGroups
        [idx, spikeGroup] = expandToGroups(find(inRec & ~ismember(clu, lost)), clu, ...
                                           sUnitGroups, siteGroup, res);
        keep = ~ismember(spikeGroup, emptyGroups);       % unsorted spikes on excluded sites
        idx = idx(keep);
        spikeGroup = spikeGroup(keep);
        localTimes = spikeTimes(idx) - 1 - offsets(r);
        written = writeResClu(xmlDir, base, nGroups, spikeGroup, localTimes, clu(idx));
        if opts.klusters
            [nS, peak] = spikeWindow(xmlFiles{r});
            written = [written, writeSpkFet(recs{r}, hCfg, nSamples(r), xmlDir, base, ...
                                            sGroups, spikeGroup, localTimes, nS, peak, ...
                                            opts.highPass)]; %#ok<AGROW>
        end
        removeStale(xmlDir, base, written);
        updateXml(xmlFiles{r}, sGroups, sUnits, sUnitGroup, notes, nPerUnit, res, ...
                  clusterSites, siteMap, opts.spikeColors);
        if opts.nrs
            cluFiles = written(~cellfun(@isempty, regexp(written, '\.clu\.\d+$', 'once')));
            updateNrs(fullfile(xmlDir, [base '.nrs']), base, cluFiles, nClusters, ...
                      strcmpi(opts.selectUnits, 'all'), opts.raster);
        end
    end

    if opts.dryRun
        fprintf('\nDry run: nothing written.\n');
    else
        fprintf(['\nDone. Open a session''s .dat in Neuroscope: its .clu files load with it,\n' ...
                 'each unit''s spikes coloured on the channels around its peak\n' ...
                 '(Units > Waveforms; untick units in the cluster palette to hide them).\n']);
    end
end

% =====================================================================
function tf = neuroscopeRunning()
    if ispc
        [~, out] = system('tasklist /FI "IMAGENAME eq neuroscope.exe" /NH');
        tf = contains(lower(out), 'neuroscope.exe');
    else
        [status, ~] = system('pgrep -xi neuroscope');
        tf = status == 0;
    end
end

% =====================================================================
function ch = xmlSkippedChannels(xmlFile)
    %XMLSKIPPEDCHANNELS 0-based channels an .xml's anatomical groups mark skip="1"
    doc = xmlread(xmlFile);
    nodes = doc.getElementsByTagName('channel');
    ch = [];
    for i = 0:nodes.getLength() - 1
        n = nodes.item(i);
        if strcmp(char(n.getAttribute('skip')), '1')
            ch(end + 1) = str2double(char(n.getTextContent())); %#ok<AGROW>
        end
    end
end

function [groups, unitGroups, unitGroup, units, emptyGroups] = ...
        excludeSites(groups, unitGroups, unitGroup, units, excluded)
    %EXCLUDESITES Take sites out of every spike group, and units out of
    %   groups left empty
    %
    %   A group keeps its number, so .clu.N means the same channel in every
    %   session; with 'radius' an excluded channel's group is simply empty.
    %   A unit whose peak group is gone keeps its first remaining group as
    %   its home in <units>; one with none left is dropped for the session.
    for g = 1:numel(groups)
        keep = ~ismember(groups(g).sites, excluded);
        groups(g).sites = groups(g).sites(keep);
        groups(g).channels = groups(g).channels(keep);
    end
    emptyGroups = find(arrayfun(@(g) isempty(g.sites), groups));
    for u = units(:)'
        unitGroups{u} = unitGroups{u}(~ismember(unitGroups{u}, emptyGroups));
        if isempty(unitGroups{u})
            unitGroup(u) = 0;
        elseif ~ismember(unitGroup(u), unitGroups{u})
            unitGroup(u) = unitGroups{u}(1);
        end
    end
    units = units(unitGroup(units) > 0);
    for g = 1:numel(groups)
        groups(g).units = units(cellfun(@(c) any(c == g), unitGroups(units)))';
    end
end

% =====================================================================
function [idx, group] = expandToGroups(spikes, clu, unitGroups, siteGroup, res)
    %EXPANDTOGROUPS One entry per spike and spike group it is drawn in
    %
    %   A unit's spikes go into every group in unitGroups{unit}; spikes of
    %   no unit (cluster 0, with 'includeNoise') into the group of the site
    %   they were detected on. Entries come unit by unit; writeResClu and
    %   writeSpkFet both sort by time, stably, so they agree on the order.
    c = clu(spikes);
    [cs, o] = sort(c);
    spikes = spikes(o);
    edges = [0; find(diff(cs)); numel(cs)];
    idxParts = cell(numel(edges) - 1, 1);
    grpParts = idxParts;
    for k = 1:numel(edges) - 1
        s = spikes(edges(k) + 1:edges(k + 1));
        u = cs(edges(k) + 1);
        if u <= 0
            g = siteGroup(double(res.spikeSites(s)));
            idxParts{k} = s;
            grpParts{k} = g(:);
        else
            g = unitGroups{u}(:)';
            idxParts{k} = repmat(s, numel(g), 1);
            grpParts{k} = repelem(g(:), numel(s), 1);       % column even for one group
        end
    end
    idx = vertcat(idxParts{:});
    group = vertcat(grpParts{:});
end

% =====================================================================
function [siteGroup, mode] = buildGroups(opts, res, hCfg, shank, siteLoc, clusterSites, units)
    %BUILDGROUPS Assign every site to one spike group, contiguous along a shank
    %
    %   'radius' gives every site its own group; which of them a unit is
    %   drawn on is decided per unit afterwards, from the probe geometry.
    nSites = numel(shank);
    siteGroup = zeros(nSites, 1);

    if ischar(opts.grouping) || isstring(opts.grouping)
        mode = lower(char(opts.grouping));
    else
        mode = 'block';
        blockSize = opts.grouping;
    end
    fprintf('Spike groups: %s', mode);
    switch mode
        case 'radius',  fprintf(' (one per channel; each unit on every channel of its shank within %g um above or below its peak site)\n', opts.radiusUm);
        case 'jrclust', fprintf(' (%d..%d sites)\n', opts.minSites, opts.maxSites);
        case 'block',   fprintf(' of %d sites\n', blockSize);
        case 'shank',   fprintf(' (one per shank)\n');
        otherwise
            error('jrc_export_neuroscope:badGrouping', ...
                  'grouping must be ''radius'', ''jrclust'', ''shank'' or a number of sites');
    end

    if strcmp(mode, 'jrclust')
        footprint = unitFootprints(res, hCfg, units);
    end

    g = 0;
    for sh = unique(shank)'
        idx = orderSites(find(shank == sh), siteLoc);
        n = numel(idx);
        switch mode
            case 'radius'
                segs = [(1:n)' (1:n)'];
            case 'shank'
                segs = [1 n];
            case 'block'
                starts = (1:blockSize:n)';
                segs = [starts min(starts + blockSize - 1, n)];
            case 'jrclust'
                onShank = units(shank(clusterSites(units)) == sh);
                [~, pos] = ismember(clusterSites(onShank), idx);
                segs = cutShank(footprint(idx, onShank), pos, ...
                                min(opts.minSites, n), min(opts.maxSites, n));
        end
        for s = 1:size(segs, 1)
            g = g + 1;
            siteGroup(idx(segs(s, 1):segs(s, 2))) = g;
        end
    end
end

function fp = unitFootprints(res, hCfg, units)
    %UNITFOOTPRINTS Share of each unit's signal energy on each site
    %
    %   From the mean waveform JRCLUST stores on the unit's neighbourhood; a
    %   result without it falls back to the neighbourhood itself, equally
    %   weighted.
    nSites = numel(hCfg.siteMap);
    fp = zeros(nSites, max(units));
    if isfield(res, 'meanWfGlobal') && size(res.meanWfGlobal, 3) >= max(units)
        wf = double(res.meanWfGlobal(:, :, units));
        vpp = squeeze(max(wf, [], 1) - min(wf, [], 1));
        fp(:, units) = reshape(vpp, nSites, []).^2;
    else
        for u = units(:)'
            fp(hCfg.siteNeighbors(:, res.clusterSites(u)), u) = 1;
        end
    end
    total = sum(fp, 1);
    total(total == 0) = 1;
    fp = fp ./ total;
end

function segs = cutShank(E, pos, minK, maxK)
    %CUTSHANK Split n ordered sites into contiguous runs of minK..maxK sites
    %
    %   E(i,u) is unit u's share of energy on site i, pos(u) the position of
    %   its peak site. A unit belongs to the run holding its peak; the cost of
    %   a split is the energy every unit loses outside its own run. Solved
    %   exactly by dynamic programming over the cut positions, with a tiny
    %   preference for runs near the middle of the allowed size where
    %   several splits lose the same.
    n = size(E, 1);
    cE = [zeros(1, size(E, 2)); cumsum(E, 1)];
    target = (minK + maxK) / 2;
    best = inf(n + 1, 1);
    best(1) = 0;
    prev = zeros(n + 1, 1);
    for b = 1:n
        for a = max(1, b - maxK + 1):(b - minK + 1)
            if ~isfinite(best(a)), continue; end
            in = pos >= a & pos <= b;
            c = best(a) + sum(1 - (cE(b + 1, in) - cE(a, in))) + 1e-6 * abs(b - a + 1 - target);
            if c < best(b + 1)
                best(b + 1) = c;
                prev(b + 1) = a;
            end
        end
    end
    if ~isfinite(best(n + 1))       % shank too short for minK: one run
        segs = [1 n];
        return;
    end
    segs = zeros(0, 2);
    b = n;
    while b > 0
        a = prev(b + 1);
        segs = [a b; segs]; %#ok<AGROW>
        b = a - 1;
    end
end

function idx = orderSites(idx, siteLoc)
    %ORDERSITES Sites from the top of the shank down, left to right at a depth
    [~, o] = sortrows([-siteLoc(idx, 2), siteLoc(idx, 1)]);
    idx = idx(o);
end

function reportGroups(mode, opts, groupSites, shank, siteLoc, siteMap, unitGroups, units, res, siteGroup)
    nGroups = numel(groupSites);
    fprintf('  %d spike groups\n', nGroups);
    nUnitsIn = arrayfun(@(g) sum(cellfun(@(c) any(c == g), unitGroups(units))), 1:nGroups);
    if strcmp(mode, 'radius')
        % one group per channel: summarise rather than list 512 rows
        nCh = cellfun(@numel, unitGroups(units));
        fprintf('  each unit is drawn on %d..%d channels (median %g) within %g um of its peak site\n', ...
                min(nCh), max(nCh), median(nCh), opts.radiusUm);
        fprintf('  %d channels show at least one unit; the busiest shows %d units\n', ...
                sum(nUnitsIn > 0), max(nUnitsIn));
        pitch = arrayfun(@(sh) median(diff(unique(siteLoc(shank == sh, 2)))), unique(shank)');
        fprintf('  vertical distance between neighbouring depths, per shank (median): %s um\n', ...
                mat2str(round(pitch, 1)));
    else
        fprintf('  group  shank  sites  channels (0-based)            units\n');
        for g = 1:nGroups
            s = groupSites{g};
            ch = siteMap(s) - 1;
            chText = sprintf('%d ', ch(1:min(5, end)));
            if numel(ch) > 5, chText = [chText '...']; end %#ok<AGROW>
            fprintf('  %5d  %5d  %5d  %-30s %5d\n', g, shank(s(1)), numel(s), chText, nUnitsIn(g));
        end
    end

    % how much of each unit's own neighbourhood the channels it is drawn on cover
    if isfield(res, 'meanWfGlobal') && size(res.meanWfGlobal, 3) >= max(units)
        kept = zeros(numel(units), 1);
        for k = 1:numel(units)
            u = units(k);
            wf = double(res.meanWfGlobal(:, :, u));
            e = (max(wf, [], 1) - min(wf, [], 1)).^2;
            kept(k) = sum(e(ismember(siteGroup, unitGroups{u}))) / max(sum(e), eps);
        end
        fprintf(['  share of each unit''s signal on the channels it is drawn on: median %.2f; ' ...
                 '%d units below 0.8, %d below 0.5\n'], median(kept), sum(kept < 0.8), sum(kept < 0.5));
    end
end

% =====================================================================
function xmlFiles = resolveXmlFiles(opts, recs)
    nRecs = numel(recs);
    if ~isempty(opts.xmlFiles)
        xmlFiles = cellstr(opts.xmlFiles);
        if numel(xmlFiles) ~= nRecs
            error('jrc_export_neuroscope:xmlCount', ...
                  'xmlFiles has %d entries but rawRecordings has %d', numel(xmlFiles), nRecs);
        end
        return;
    end

    % Each recording's own .xml: same folder, same name as its .dat, which is
    % the only .xml Neuroscope reads when it opens that .dat. Asked for only
    % when it is missing (or always, with 'ask', true).
    xmlFiles = cell(nRecs, 1);
    for r = 1:nRecs
        [d, base] = fileparts(recs{r});
        guess = fullfile(d, [base '.xml']);
        if isfile(guess) && ~isequal(opts.ask, true)
            xmlFiles{r} = guess;
        elseif usejava('desktop') && ~isequal(opts.ask, false)
            [f, p] = uigetfile('*.xml', sprintf('Neuroscope .xml for recording %d of %d: %s', ...
                                                r, nRecs, recs{r}), guess);
            if ~isequal(f, 0), xmlFiles{r} = fullfile(p, f); end
        else
            warning('jrc_export_neuroscope:noXml', '%s not found, recording skipped', guess);
        end
    end
end

function checkXml(xmlFile, hCfg, siteMap)
    %CHECKXML The .xml must describe the same recording the sort was run on
    doc = xmlread(xmlFile);
    nCh = str2double(textOf(doc, 'nChannels'));
    rate = str2double(textOf(doc, 'samplingRate'));
    if nCh ~= hCfg.nChans
        error('jrc_export_neuroscope:nChannels', '%s has %d channels, the .prm %d', ...
              xmlFile, nCh, hCfg.nChans);
    end
    if abs(rate - hCfg.sampleRate) > 1e-6
        error('jrc_export_neuroscope:rate', '%s samples at %g Hz, the .prm at %g Hz', ...
              xmlFile, rate, hCfg.sampleRate);
    end
    if max(siteMap) > nCh
        error('jrc_export_neuroscope:siteMap', 'siteMap refers to channel %d of %d', max(siteMap), nCh);
    end
end

function t = textOf(doc, tag)
    nodes = doc.getElementsByTagName(tag);
    if nodes.getLength() == 0
        t = '';
    else
        t = char(nodes.item(0).getTextContent());
    end
end

% =====================================================================
function written = writeResClu(dirName, base, nGroups, spikeGroup, times, clu)
    %WRITERESCLU One .res.N / .clu.N pair per group that holds any unit
    %
    %   .res: one 0-based sample number per line. .clu: the number of
    %   clusters in the file, then one cluster id per spike, in .res order.
    %   Neuroscope rejects a pair whose line counts differ, and only reads a
    %   number once it sees the character after it, so every line ends in a
    %   newline, the last one included.
    written = {};
    for g = 1:nGroups
        in = spikeGroup == g;
        if ~any(in), continue; end
        [t, o] = sort(times(in));
        c = clu(in);
        c = c(o);
        resFile = fullfile(dirName, sprintf('%s.res.%d', base, g));
        cluFile = fullfile(dirName, sprintf('%s.clu.%d', base, g));
        writeColumn(resFile, '', int64(t));
        writeColumn(cluFile, sprintf('%d\n', numel(unique(c))), int64(c));
        written(end + 1:end + 2) = {resFile, cluFile};
        fprintf('  group %2d: %8d spikes, units %s\n', g, numel(t), compactList(unique(c)));
    end
end

function writeColumn(fileName, header, values)
    fid = fopen(fileName, 'w');
    if fid < 0, error('jrc_export_neuroscope:write', 'cannot write %s', fileName); end
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, '%s', header);
    fprintf(fid, '%d\n', values);
end

function removeStale(dirName, base, written)
    %REMOVESTALE Delete .res/.clu files an earlier export wrote and this one did not
    %
    %   Only files listed in this script's own manifest are touched, so .clu
    %   files from Klusters or anything else in the folder are left alone.
    manifest = fullfile(dirName, [base '.jrc_neuroscope.txt']);
    if isfile(manifest)
        old = strtrim(splitlines(fileread(manifest)));
        old = old(~cellfun(@isempty, old));
        for f = setdiff(old, written)'
            if isfile(f{1})
                delete(f{1});
                fprintf('  removed %s (left over from a previous export)\n', f{1});
            end
        end
    end
    fid = fopen(manifest, 'w');
    fprintf(fid, '%s\n', written{:});
    fclose(fid);
end

function [nS, peak] = spikeWindow(xmlFile)
    %SPIKEWINDOW Waveform length and peak position Neuroscope uses, default 32/16
    doc = xmlread(xmlFile);
    nS = str2double(textOf(doc, 'nSamples'));
    peak = str2double(textOf(doc, 'peakSampleIndex'));
    if ~isfinite(nS) || nS < 2, nS = 32; end
    if ~isfinite(peak) || peak < 1 || peak > nS, peak = 16; end
end

function written = writeSpkFet(datFile, hCfg, nFrames, dirName, base, groups, ...
                               spikeGroup, times, nS, peak, highPass)
    %WRITESPKFET Klusters' waveform (.spk.N) and feature (.fet.N) files
    %
    %   .spk: int16 waveforms, nS samples per spike with every channel of the
    %   group at each sample (channel index fastest), the .res sample at
    %   position `peak` - the layout ndManager's own .spk files have.
    %   Waveforms are cut from the .dat, high-pass filtered, read once front
    %   to back in chunks.
    %
    %   .fet: the number of features, then per spike 3 principal components
    %   per channel (channel by channel) and the spike time last, as integers.
    nGroups = numel(groups);
    nCh = hCfg.nChans;
    has = arrayfun(@(g) any(spikeGroup == g), 1:nGroups);
    written = {};
    if ~any(has), return; end

    % spikes in .res order: by time, ties in their original order
    [times, o] = sort(times);
    spikeGroup = spikeGroup(o);
    offs = (1:nS) - peak;                                  % window around the spike

    fids = zeros(nGroups, 1);
    spkFiles = cell(nGroups, 1);
    for g = find(has)
        spkFiles{g} = fullfile(dirName, sprintf('%s.spk.%d', base, g));
        fids(g) = fopen(spkFiles{g}, 'w');
    end
    cleanup = onCleanup(@() arrayfun(@(f) fclose(f), fids(fids > 0)));

    [b, a] = butter(3, highPass / (hCfg.sampleRate / 2), 'high');
    pad = 2048;                                            % filter settling, > window
    chunk = 2^18;
    fid = fopen(datFile, 'r');
    datCleanup = onCleanup(@() fclose(fid));
    dataType = ['*' char(hCfg.dataType)];
    bytesPerFrame = nCh * numel(typecast(cast(0, hCfg.dataType), 'uint8'));

    fprintf('  .spk: reading %s (%.1f GB)\n', datFile, nFrames * bytesPerFrame / 2^30);
    tic;
    lastShown = 0;
    first = 1;                                             % next spike to place
    for s = 0:chunk:nFrames - 1
        e = min(s + chunk, nFrames);
        last = first - 1 + find(times(first:end) < e, 1, 'last');
        if isempty(last) || last < first, continue; end
        k = first:last;
        first = last + 1;

        % frames s-pad .. e+pad, as far as the file goes. No zero padding at
        % the ends: the .dat carries a DC offset, and a step from zero to it
        % rings through the filter right where the first spikes sit.
        r0 = max(s - pad, 0);
        r1 = min(e + pad, nFrames);
        fseek(fid, hCfg.headerOffset + r0 * bytesPerFrame, 'bof');
        X = fread(fid, [nCh, r1 - r0], dataType);
        X = filtfilt(b, a, double(X)');                    % frames x channels

        % a window running off either end of the file repeats the edge sample
        rows = min(max((times(k) - r0) + 1 + offs, 1), r1 - r0);   % nSpikes x nS
        for g = unique(spikeGroup(k))'
            in = spikeGroup(k) == g;
            ch = groups(g).channels + 1;
            W = X(reshape(rows(in, :)', [], 1), ch);       % (nS*nIn) x nCh
            W = reshape(W, nS, sum(in), numel(ch));        % nS x nIn x nCh
            fwrite(fids(g), int16(round(permute(W, [3 1 2]))), 'int16');
        end

        if toc - lastShown > 10
            lastShown = toc;
            fprintf('    %3.0f%%  %.0f s\n', 100 * e / nFrames, toc);
        end
    end
    clear cleanup datCleanup;
    fprintf('    done in %.0f s\n', toc);

    for g = find(has)
        nChG = numel(groups(g).channels);
        fid = fopen(spkFiles{g}, 'r');
        W = fread(fid, inf, '*int16');
        fclose(fid);
        nSpk = numel(W) / (nChG * nS);
        W = reshape(single(W), nChG, nS, nSpk);
        t = times(spikeGroup == g);

        F = zeros(nSpk, 3 * nChG, 'single');
        sub = randperm(nSpk, min(nSpk, 20000));
        for c = 1:nChG
            Y = reshape(W(c, :, :), nS, nSpk)';             % nSpk x nS
            mu = mean(Y(sub, :), 1);
            [~, ~, V] = svd(double(Y(sub, :) - mu), 'econ');
            nPc = min(3, size(V, 2));
            F(:, 3*(c-1) + (1:nPc)) = (Y - mu) * single(V(:, 1:nPc));
        end
        scale = 1000 / max(prctile(abs(F(:)), 99.9), eps);
        fetFile = fullfile(dirName, sprintf('%s.fet.%d', base, g));
        fid = fopen(fetFile, 'w');
        fprintf(fid, '%d\n', 3 * nChG + 1);
        fprintf(fid, [repmat('%d ', 1, 3 * nChG) '%d\n'], ...
                [int64(round(F * scale)), int64(t(:))]');
        fclose(fid);
        written(end + 1:end + 2) = {spkFiles{g}, fetFile};
    end
    fprintf('  .fet: %d groups, 3 PCs per channel\n', sum(has));
end

function insertAt(parent, el, anchor)
    if isempty(anchor)
        parent.appendChild(el);
    else
        parent.insertBefore(el, anchor);
    end
end

function updateNrs(nrsFile, base, cluFiles, nClusters, selectAll, raster)
    %UPDATENRS List this export's .clu files in the Neuroscope session file
    %
    %   Neuroscope loads every cluster file a session lists when the .dat is
    %   opened. Entries for <base>.clu.N files are replaced; event files,
    %   displays and everything else stay. Each unit gets its jrc manual
    %   colour, and each entry the file's own modification time - a missing
    %   colour or a different time makes Neuroscope warn that the file
    %   changed since the session was saved.
    if ~isfile(nrsFile)
        fprintf('  no %s - open the .clu files once in Neuroscope and save the session\n', nrsFile);
        return;
    end
    backup = sprintf('%s.%s.bak', nrsFile, datestr(now, 'yyyymmdd_HHMMSS')); %#ok<TNOW1,DATST>
    copyfile(nrsFile, backup);

    doc = xmlread(nrsFile);
    root = doc.getDocumentElement();
    stripWhitespace(root);
    files = directChild(root, 'files');
    if isempty(files)
        files = doc.createElement('files');
        root.insertBefore(files, root.getFirstChild());
    end

    % drop entries for this base's .clu files from an earlier export
    pattern = ['(^|[\\/])' regexptranslate('escape', base) '\.clu\.\d+$'];
    kids = files.getChildNodes();
    for i = kids.getLength() - 1:-1:0
        f = kids.item(i);
        if f.getNodeType() ~= f.ELEMENT_NODE, continue; end
        url = directChild(f, 'url');
        if ~isempty(url) && ~isempty(regexp(strrep(char(url.getTextContent()), '%5C', '\'), pattern, 'once'))
            files.removeChild(f);
        end
    end

    wholeSeconds = truncateModTimes(cluFiles);
    colors = lines(max(nClusters, 1));
    urls = strrep(cluFiles, '\', '/');
    fileIds = cell(size(cluFiles));
    for i = 1:numel(cluFiles)
        f = cluFiles{i};
        d = dir(f);
        stamp = datestr(d.datenum, 'yyyy-mm-ddTHH:MM:SS'); %#ok<DATST>

        fid = fopen(f, 'r');
        fgetl(fid);
        ids = unique(fscanf(fid, '%d'));
        fclose(fid);
        fileIds{i} = ids;

        entry = add(doc, files, 'file');
        add(doc, entry, 'type', '0');
        add(doc, entry, 'url', urls{i});
        add(doc, entry, 'modificationDate', stamp);
        items = add(doc, entry, 'items');
        for id = ids'
            d = add(doc, items, 'itemDescription');
            add(doc, d, 'item', sprintf('%d', id));
            if id >= 1
                rgb = round(255 * colors(mod(id - 1, size(colors, 1)) + 1, :));
            else
                rgb = [255 0 0];                           % unsorted spikes, red
            end
            add(doc, d, 'color', sprintf('#%02x%02x%02x', rgb));
        end
    end

    % Every display: spikes drawn as waveforms on the traces, and the units
    % of this export selected, as in an ndManager session.
    %
    % Each <spikePresentation> is one mode that is on: 0 raster, 1 waveforms
    % (coloured on the traces), 2 vertical lines. Neuroscope starts with only
    % the raster. With every unit selected the raster would hold a row per
    % unit and channel file - thousands - so it is on only with 'raster'.
    %
    % Selections are stored per .clu file; those of an earlier export point
    % at files whose contents have changed, so they are replaced.
    modes = {'1'};
    if raster, modes{end + 1} = '0'; end
    displays = doc.getElementsByTagName('display');
    for i = 0:displays.getLength() - 1
        disp_ = displays.item(i);
        kids = disp_.getChildNodes();
        anchor = [];
        for j = kids.getLength() - 1:-1:0
            k = kids.item(j);
            if k.getNodeType() ~= k.ELEMENT_NODE, continue; end
            tag = char(k.getNodeName());
            if strcmp(tag, 'spikePresentation')
                anchor = k.getNextSibling();
                disp_.removeChild(k);
            elseif any(strcmp(tag, {'clustersSelected', 'clustersSkipped'}))
                u = directChild(k, 'fileUrl');
                if ~isempty(u) && ~isempty(regexp(strrep(char(u.getTextContent()), '%5C', '\'), pattern, 'once'))
                    anchor = k.getNextSibling();
                    disp_.removeChild(k);
                end
            end
        end
        if isempty(anchor), anchor = directChild(disp_, 'channelPositions'); end
        for m = modes
            el = doc.createElement('spikePresentation');
            el.appendChild(doc.createTextNode(m{1}));
            insertAt(disp_, el, anchor);
        end
        if selectAll
            for f = 1:numel(urls)
                sel = doc.createElement('clustersSelected');
                add(doc, sel, 'fileUrl', urls{f});
                for id = fileIds{f}'
                    add(doc, sel, 'cluster', sprintf('%d', id));
                end
                insertAt(disp_, sel, anchor);
            end
        end
    end

    xmlwrite(nrsFile, doc);
    fprintf('  registered %d cluster files in %s; waveform display on', numel(cluFiles), nrsFile);
    if selectAll, fprintf(', all units selected'); end
    fprintf('\n');
    if ~wholeSeconds
        warning('jrc_export_neuroscope:modTime', ...
                ['could not round the .clu modification times to whole seconds; ' ...
                 'Neuroscope may report each file as modified since the session was saved']);
    end
end

function ok = truncateModTimes(files)
    %TRUNCATEMODTIMES Drop the sub-second part of each file's modification time
    %
    %   The .nrs stores a file's time to the second and Neuroscope compares
    %   it with the file's own time to the millisecond, so a file written at
    %   12:27:26.431 always looks modified. MATLAB cannot set file times, so
    %   this asks the operating system.
    ok = true;
    batch = 16;                                            % keeps each command line short
    for i0 = 1:batch:numel(files)
        part = files(i0:min(i0 + batch - 1, end));
        if ispc
            quoted = strjoin(cellfun(@(f) ['''' strrep(f, '''', '''''') ''''], part, ...
                                     'UniformOutput', false), ',');
            cmd = sprintf(['powershell -NoProfile -NonInteractive -Command "foreach ($f in @(%s)) ' ...
                           '{ $i = Get-Item -LiteralPath $f; $t = $i.LastWriteTime; ' ...
                           '$i.LastWriteTime = $t.AddTicks(-($t.Ticks %% 10000000)) }"'], quoted);
        else
            cmd = strjoin(cellfun(@(f) sprintf('touch -d "$(date -r ''%s'' ''+%%Y-%%m-%%d %%H:%%M:%%S'')" ''%s''', ...
                                               f, f), part, 'UniformOutput', false), ' && ');
        end
        [status, ~] = system(cmd);
        ok = ok && status == 0;
    end
    for i = 1:numel(files)
        d = dir(files{i});
        ok = ok && abs(d.datenum * 86400 - round(d.datenum * 86400)) < 1e-3;
    end
end

function s = compactList(v)
    v = v(:)';
    if numel(v) > 12
        s = [sprintf('%d ', v(1:12)) sprintf('... (%d)', numel(v))];
    else
        s = strtrim(sprintf('%d ', v));
    end
end

% =====================================================================
function updateXml(xmlFile, groups, units, unitGroup, notes, nPerUnit, res, ...
                   clusterSites, siteMap, setSpikeColors)
    %UPDATEXML Rewrite <spikeDetection> and <units>, keep everything else
    backup = sprintf('%s.%s.bak', xmlFile, datestr(now, 'yyyymmdd_HHMMSS')); %#ok<TNOW1,DATST>
    copyfile(xmlFile, backup);
    fprintf('  backed up the .xml to %s\n', backup);

    doc = xmlread(xmlFile);
    root = doc.getDocumentElement();
    stripWhitespace(root);

    % nSamples / peakSampleIndex per group, from Neuroscope's own setting
    nSamples = textOf(doc, 'nSamples');
    peak = textOf(doc, 'peakSampleIndex');
    if isempty(nSamples), nSamples = '32'; end
    if isempty(peak), peak = '16'; end

    spike = replaceSection(doc, root, 'spikeDetection', {'units', 'neuroscope', 'programs'});
    cg = add(doc, spike, 'channelGroups');
    for g = 1:numel(groups)
        grp = add(doc, cg, 'group');
        chs = add(doc, grp, 'channels');
        for ch = groups(g).channels
            add(doc, chs, 'channel', sprintf('%d', ch));
        end
        add(doc, grp, 'nSamples', nSamples);
        add(doc, grp, 'peakSampleIndex', peak);
        add(doc, grp, 'nFeatures', '3');
    end

    unitsNode = replaceSection(doc, root, 'units', {'neuroscope', 'programs'});
    for u = units(:)'
        un = add(doc, unitsNode, 'unit');
        add(doc, un, 'group', sprintf('%d', unitGroup(u)));
        add(doc, un, 'cluster', sprintf('%d', u));
        add(doc, un, 'structure', '');
        add(doc, un, 'type', char(notes{u}));
        iso = '';
        if isfield(res, 'unitIsoDist') && u <= numel(res.unitIsoDist) && isfinite(res.unitIsoDist(u))
            iso = sprintf('%.2f', res.unitIsoDist(u));
        end
        add(doc, un, 'isolationDistance', iso);
        q = '';
        if isfield(res, 'unitSNR') && u <= numel(res.unitSNR) && isfinite(res.unitSNR(u))
            q = sprintf('SNR %.1f', res.unitSNR(u));
        end
        add(doc, un, 'quality', q);
        add(doc, un, 'notes', sprintf('JRCLUST unit %d, peak site %d = channel %d, %d spikes', ...
                                      u, clusterSites(u), siteMap(clusterSites(u)) - 1, nPerUnit(u)));
    end

    if setSpikeColors
        colorSpikeGroups(doc, groups);
    end

    xmlwrite(xmlFile, doc);
    fprintf('  wrote %d spike groups and %d units into %s\n', numel(groups), numel(units), xmlFile);
end

function node = replaceSection(doc, root, tag, before)
    %REPLACESECTION Empty <tag> under root, creating it ahead of `before` if absent
    node = directChild(root, tag);
    if isempty(node)
        node = doc.createElement(tag);
        ref = [];
        for b = before
            ref = directChild(root, b{1});
            if ~isempty(ref), break; end
        end
        if isempty(ref)
            root.appendChild(node);
        else
            root.insertBefore(node, ref);
        end
    else
        while node.hasChildNodes()
            node.removeChild(node.getFirstChild());
        end
    end
end

function colorSpikeGroups(doc, groups)
    %COLORSPIKEGROUPS Give each spike group's channels one spikeColor
    %
    %   Neighbouring groups get clearly different colours, so the spike-group
    %   palette shows where one group ends and the next begins. With one
    %   group per channel ('radius') the colour follows the shank instead.
    palette = {'#ff3b30', '#34c759', '#007aff', '#ff9500', '#af52de', '#00c7be', ...
               '#ffcc00', '#ff2d55', '#5856d6', '#a2845e'};
    chColor = containers.Map('KeyType', 'double', 'ValueType', 'char');
    for g = 1:numel(groups)
        for ch = groups(g).channels
            chColor(ch) = palette{mod(groups(g).colorKey - 1, numel(palette)) + 1};
        end
    end
    entries = doc.getElementsByTagName('channelColors');
    for i = 0:entries.getLength() - 1
        e = entries.item(i);
        chNode = directChild(e, 'channel');
        if isempty(chNode), continue; end
        ch = str2double(char(chNode.getTextContent()));
        if ~isKey(chColor, ch), continue; end
        sc = directChild(e, 'spikeColor');
        if isempty(sc), sc = add(doc, e, 'spikeColor'); end
        sc.setTextContent(chColor(ch));
    end
end

function child = directChild(parent, tag)
    child = [];
    kids = parent.getChildNodes();
    for i = 0:kids.getLength() - 1
        k = kids.item(i);
        if k.getNodeType() == k.ELEMENT_NODE && strcmp(char(k.getNodeName()), tag)
            child = k;
            return;
        end
    end
end

function el = add(doc, parent, tag, text)
    el = doc.createElement(tag);
    if nargin > 3 && ~isempty(text)
        el.appendChild(doc.createTextNode(text));
    end
    parent.appendChild(el);
end

function stripWhitespace(node)
    %STRIPWHITESPACE Drop formatting-only text so xmlwrite indents cleanly
    kids = node.getChildNodes();
    for i = kids.getLength() - 1:-1:0
        k = kids.item(i);
        if k.getNodeType() == k.TEXT_NODE && isempty(strtrim(char(k.getNodeValue())))
            node.removeChild(k);
        elseif k.getNodeType() == k.ELEMENT_NODE
            stripWhitespace(k);
        end
    end
end
