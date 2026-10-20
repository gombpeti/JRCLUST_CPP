function addMenu(obj, hFig)
    drawnow;

    outerPosition = hFig.outerPosition;
    hFig.figApply(@set, 'MenuBar','None');

    obj.hMenus('FileMenu') = hFig.uimenu('Label', 'File');
    uimenu(obj.hMenus('FileMenu'), 'Label', 'Save', 'Callback', @(hO, hE) obj.saveFiles());
    uimenu(obj.hMenus('FileMenu'), 'Label', 'Save As Pruned', 'Callback', @(hO, hE) saveAsPruned(obj));
    uimenu(obj.hMenus('FileMenu'), 'Label', 'Save figures as .fig', 'Callback', @(hO, hE) obj.saveFigures('.fig'));
    uimenu(obj.hMenus('FileMenu'), 'Label', 'Save figures as .png', 'Callback', @(hO, hE) obj.saveFigures('.png'));
    uimenu(obj.hMenus('FileMenu'), 'Label', 'Export units to csv', 'Callback', @(hO, hE) obj.hClust.exportToCSV(), 'Separator', 'on');
    uimenu(obj.hMenus('FileMenu'), 'Label', 'Export unit qualities to csv', 'Callback', @(hO, hE) obj.hClust.exportQualityScores());
    uimenu(obj.hMenus('FileMenu'), 'Label', 'Export hClust to workspace', 'Callback', @(hO, hE) assignin('base', 'hClust', obj.hClust));
    uimenu(obj.hMenus('FileMenu'), 'Label', 'Export all mean unit waveforms to workspace', 'Callback', @(hO, hE) obj.exportMeanWf(1));
    uimenu(obj.hMenus('FileMenu'), 'Label', 'Export selected mean unit waveforms to workspace', 'Callback', @(hO, hE) obj.exportMeanWf(0));
    uimenu(obj.hMenus('FileMenu'), 'Label', 'Export all traces from the selected unit', 'Callback', @(hO, hE) obj.exportTraces());
    uimenu(obj.hMenus('FileMenu'), 'Label', 'Export firing rate for all units', 'Callback', @(hO, hE) obj.exportFiringRate());
    uimenu(obj.hMenus('FileMenu'), 'Label', 'Print summary', 'Callback', @(hO, hE) obj.summarize(), 'Separator', 'on');
    uimenu(obj.hMenus('FileMenu'), 'Label', 'Exit', 'Callback', @(hO, hE) obj.endSession(), 'Separator', 'on', 'Accelerator', 'Q');

    obj.hMenus('EditMenu') = hFig.uimenu('Label', 'Edit');
    uimenu(obj.hMenus('EditMenu'), 'Label', '[M]erge', 'Callback', @(hO, hE) obj.mergeSelected());
    uimenu(obj.hMenus('EditMenu'), 'Label', 'Merge auto', 'Callback', @(hO, hE) obj.autoMerge());
    uimenu(obj.hMenus('EditMenu'), 'Label', 'Merge auto (cross-correlogram checked)', 'Callback', @(hO, hE) autoMergeChecked(obj));
    uimenu(obj.hMenus('EditMenu'), 'Label', '[D]elete', 'Callback', @(hO, hE) obj.deleteClusters(), 'Separator', 'on');
    uimenu(obj.hMenus('EditMenu'), 'Label', 'Delete auto', 'Callback', @(hO, hE) obj.autoDelete());
    uimenu(obj.hMenus('EditMenu'), 'Label', 'Delete annotated', 'Callback', @(hO, hE) obj.deleteAnnotated()); % TW
    uimenu(obj.hMenus('EditMenu'), 'Label', 'Relabel units (noise / to_delete / review)', 'Callback', @(hO, hE) relabelUnits(obj), 'Separator', 'on');
    uimenu(obj.hMenus('EditMenu'), 'Label', '[S]plit', 'Callback', @(hO, hE) obj.autoSplit(1), 'Separator', 'on');
    uimenu(obj.hMenus('EditMenu'), 'Label', 'Auto split max-chan', 'Callback', @(hO, hE) obj.autoSplit(0));
    uimenu(obj.hMenus('EditMenu'), 'Label', 'Auto split multi-chan', 'Callback', @(hO, hE) obj.autoSplit(1));
    uimenu(obj.hMenus('EditMenu'), 'Label', 'Recompute selected cluster mean waveform', 'Callback', @(hO, hE) obj.recomputeWaveforms(), 'Separator', 'on');
    uimenu(obj.hMenus('EditMenu'), 'Label', 'Recompute all mean waveforms', 'Callback', @(hO, hE) obj.recomputeWaveforms(1));
    uimenu(obj.hMenus('EditMenu'), 'Label', 'Reorder clusters by center site', 'Callback', @(hO, hE) obj.reorderClusters('clusterSites'), 'Separator', 'on');
    uimenu(obj.hMenus('EditMenu'), 'Label', 'Reorder clusters by centroid', 'Callback', @(hO, hE) obj.reorderClusters('Y + X'));
    uimenu(obj.hMenus('EditMenu'), 'Label', 'Reorder sites by channel number', 'Callback', @(h0, hE) obj.reorderSites(1:obj.hCfg.nSites),'Separator','on');    
    uimenu(obj.hMenus('EditMenu'), 'Label', 'Reorder sites by location dimension 1 (X)', 'Callback', @(h0, hE) obj.reorderSites(obj.hCfg.siteLoc(:,1)));
    uimenu(obj.hMenus('EditMenu'), 'Label', 'Reorder sites by location dimension 2 (Y)', 'Callback', @(h0, hE) obj.reorderSites(obj.hCfg.siteLoc(:,2)));

    obj.hMenus('ViewMenu') = hFig.uimenu( 'Label', 'View');
    uimenu(obj.hMenus('ViewMenu'), 'Label', 'Show traces', 'Callback', @(hO, hE) obj.showTraces());
    uimenu(obj.hMenus('ViewMenu'), 'Label', 'View all [R]', 'Callback', @(hO, hE) obj.keyPressFigWav([], struct('Key', 'r')));
    uimenu(obj.hMenus('ViewMenu'), 'Label', '[Z]oom selected', 'Callback', @(hO, hE) obj.keyPressFigWav([], struct('Key', 'z')), 'Separator', 'on');
    uimenu(obj.hMenus('ViewMenu'), 'Label', '[W]aveform (toggle)', 'Callback', @(hO, hE) obj.keyPressFigWav([], struct('Key', 'w')));
    uimenu(obj.hMenus('ViewMenu'), 'Label', '[N]umbers (toggle)', 'Callback', @(hO, hE) obj.keyPressFigWav([], struct('Key', 'n')));
    uimenu(obj.hMenus('ViewMenu'), 'Label', 'Show raw waveform', 'Callback', @(hO, hE) obj.toggleRaw(hO), 'Separator', 'on');
    uimenu(obj.hMenus('ViewMenu'), 'Label', 'Display annotated units only', 'Callback', @(hO, hE) obj.updateSubset(obj.hClust.annotatedOnly), 'Separator', 'on');
    uimenu(obj.hMenus('ViewMenu'), 'Label', 'Display non-annotated units only', 'Callback', @(hO, hE) obj.updateSubset(setdiff(1:obj.hClust.nClusters, obj.hClust.annotatedOnly)));
    uimenu(obj.hMenus('ViewMenu'), 'Label', 'Display all units', 'Callback', @(hO, hE) obj.updateSubset(1:obj.hClust.nClusters, 1));
    %uimenu(obj.hMenus('ViewMenu'), 'Label', 'Threshold by sites', 'Callback', @(hO, hE) keyPressFcn_thresh_(hFig, 'n'));
    uimenu(obj.hMenus('ViewMenu'), 'Label', 'Reset window positions', 'Callback', @(hO, hE) obj.resetPositions(), 'Separator', 'on');
    uimenu(obj.hMenus('ViewMenu'), 'Label', 'Show config file', 'Callback', @(hO, hE) obj.hCfg.edit(), 'Separator', 'on');

    obj.hMenus('ProjMenu') = hFig.uimenu( 'Label', 'Projection');
    uimenu(obj.hMenus('ProjMenu'), 'Label', 'vpp', 'Callback', @(hO, hE) obj.updateProjection('vpp'));
    uimenu(obj.hMenus('ProjMenu'), 'Label', 'pca', 'Callback', @(hO, hE) obj.updateProjection('pca'));
    uimenu(obj.hMenus('ProjMenu'), 'Label', 'ppca', 'Callback', @(hO, hE) obj.updateProjection('ppca'));
    % uimenu(obj.hMenus('ProjMenu'), 'Label', 'cov', 'Callback', @(hO, hE) obj.updateProjection('cov'));

    obj.hMenus('PlotMenu') = hFig.uimenu( 'Label', 'Plot');
    uimenu(obj.hMenus('PlotMenu'), 'Label', 'Firing rate vs. aux. input (all units)', 'Callback', @(hO, hE) obj.plotAuxRate(0));
    uimenu(obj.hMenus('PlotMenu'), 'Label', 'Firing rate vs. aux. input (selected unit)', 'Callback', @(hO, hE) obj.plotAuxRate(1));

    obj.hMenus('InfoMenu') = hFig.uimenu( 'Label', '', 'Tag', 'InfoMenu');
    uimenu(obj.hMenus('InfoMenu'), 'Label', 'Annotate unit', 'Callback', @(hO, hE) obj.annotateUnit('', 1));
    uimenu(obj.hMenus('InfoMenu'), 'Label', 'Single unit', 'Callback', @(hO, hE) obj.annotateUnit('single', 0), 'Accelerator', '1');
    uimenu(obj.hMenus('InfoMenu'), 'Label', 'Multi unit', 'Callback', @(hO, hE) obj.annotateUnit('multi', 0), 'Accelerator', '2');
    uimenu(obj.hMenus('InfoMenu'), 'Label', 'Noise', 'Callback', @(hO, hE) obj.annotateUnit('noise', 0));
    uimenu(obj.hMenus('InfoMenu'), 'Label', 'Clear annotation', 'Callback', @(hO, hE) obj.annotateUnit('', 0));
    uimenu(obj.hMenus('InfoMenu'), 'Label', 'Equal to', 'Callback', @(hO, hE) obj.annotateUnit('=', 1));

    obj.hMenus('HistMenu') = hFig.uimenu( 'Label', 'History', 'Tag', 'HistMenu');

    obj.hMenus('HelpMenu') = hFig.uimenu( 'Label', 'Help');
    uimenu(obj.hMenus('HelpMenu'), 'Label', '[H]elp', 'Callback', @(hO, hE) obj.keyPressFigWav([], struct('Key', 'h')));

    drawnow;
    hFig.outerPosition = outerPosition;
end

function autoMergeChecked(obj)
    %AUTOMERGECHECKED Edit > Merge auto, with the cross-correlogram as referee.
    %
    %   Merge auto above thresholds mean-waveform correlation alone, which
    %   cannot separate one neuron split in two from two neurons sitting at
    %   the same distance from the same sites: their mean waveforms are
    %   equally alike either way. That is why its threshold has to be set so
    %   high, and why over-split units survive it.
    %
    %   This asks the spike trains instead. For each candidate pair, count
    %   spike pairs at short lag against pairs at 10-25 ms:
    %
    %     near 0    the two units avoid each other at short lag, which only
    %               one shared refractory period explains -> one cell
    %     near 1    no interaction -> two independent cells
    %     above 5   they fire at the same instant -> duplicate detections
    %
    %   Measured on an 841-unit recording, a correlation cut of 0.95 proposed
    %   200 pairs and only 91 survived this test - so a plain merge at that
    %   setting would have fused 109 pairs of different neurons. A merge
    %   cannot be undone afterwards, which is why the preview below is shown
    %   before anything is written.
    %
    %   It is weakest where the two units differ greatly in size: a few
    %   thousand spikes against a few hundred thousand cannot move the
    %   correlogram, so such pairs are left alone and are better judged by
    %   eye.
    if exist('jrc_automerge', 'file') ~= 2
        errordlg({'This needs the helper scripts on the MATLAB path.', ...
                  '', 'Run this once, then reopen the GUI:', '', ...
                  '    addpath(''D:\JRCLUST_CPP\matlab'')'}, 'Merge auto (checked)');
        return;
    end

    prmFile = obj.hCfg.configFile;
    if isempty(prmFile) || ~isfile(prmFile)
        errordlg('Cannot find this session''s .prm file.', 'Merge auto (checked)');
        return;
    end

    dlgAns = inputdlg({['Waveform correlation above which a pair is a ' ...
                        'CANDIDATE (0-1):'], ...
                       ['Cross-correlogram ratio below which a candidate is ' ...
                        'accepted as one cell (1.0 means no evidence - keep ' ...
                        'this below 1):']}, ...
                      'Merge auto (cross-correlogram checked)', 1, ...
                      {num2str(obj.hCfg.maxUnitSim), '0.50'});
    if isempty(dlgAns)
        return;
    end

    maxSim   = str2double(dlgAns{1});
    maxRatio = str2double(dlgAns{2});
    if ~isfinite(maxSim) || maxSim <= 0 || maxSim > 1 || ...
       ~isfinite(maxRatio) || maxRatio <= 0
        jrclust.utils.qMsgBox('Invalid criteria.');
        return;
    end
    if maxRatio >= 1
        if ~strcmp(questdlg({sprintf(['A ratio cut of %.2f accepts pairs ' ...
                'that show no refractory interaction at all.'], maxRatio), '', ...
                ['1.0 IS the value for two independent cells, so this ' ...
                 'effectively disables the check and merges on waveform ' ...
                 'shape alone.'], '', 'Continue anyway?'}, ...
                'Merge auto (checked)', 'Yes', 'No', 'No'), 'Yes')
            return;
        end
    end

    % The merge reads the _res.mat from disk, so curation held only in the
    % GUI has to be written out first or it would be merged against a stale
    % clustering and then overwritten.
    obj.saveFiles();

    % Always preview. A merge is irreversible once saved, and the count of
    % REJECTED pairs is the number that matters: it says how much damage the
    % unrefereed version would have done at the same correlation cut.
    hBox = jrclust.utils.qMsgBox('Evaluating pairs... (this closes automatically)', 0, 1);
    try
        out = evalc(['jrc_automerge(''' prmFile ''', ''maxUnitSim'', ' ...
                     num2str(maxSim) ', ''maxCrossRatio'', ' num2str(maxRatio) ...
                     ', ''dryRun'', true)']);
    catch ME
        jrclust.utils.tryClose(hBox);
        errordlg(sprintf('Evaluation failed:\n\n%s', ME.message), 'Merge auto (checked)');
        return;
    end
    jrclust.utils.tryClose(hBox);

    cand = regexp(out, 'candidates\s+(\d+)', 'tokens', 'once');
    acc  = regexp(out, 'accepted as one cell\s+(\d+)', 'tokens', 'once');
    rej  = regexp(out, 'rejected \(independent or duplicate\)\s+(\d+)', 'tokens', 'once');
    uu   = regexp(out, '(\d+) units -> (\d+) units', 'tokens', 'once');

    if isempty(acc) || str2double(acc{1}) == 0
        helpdlg({'No pairs passed the cross-correlogram check.', '', ...
                 'Lower the correlation cut to propose more candidates, or ' ...
                 'merge by hand.'}, 'Merge auto (checked)');
        return;
    end

    if isempty(uu), after = '?'; else, after = uu{2}; end
    if ~strcmp(questdlg({sprintf('Candidates (correlation >= %.2f) : %s', maxSim, cand{1}), ...
                         sprintf('Accepted as one cell            : %s', acc{1}), ...
                         sprintf('Rejected as different cells     : %s', rej{1}), '', ...
                         sprintf('%s units  ->  %s units', num2str(obj.hClust.nClusters), after), '', ...
                         'Apply? This cannot be undone once saved.'}, ...
                        'Merge auto (checked)', 'Yes', 'No', 'No'), 'Yes')
        return;
    end

    hBox = jrclust.utils.qMsgBox('Merging... (this closes automatically)', 0, 1);
    try
        jrc_automerge(prmFile, 'maxUnitSim', maxSim, 'maxCrossRatio', maxRatio);
    catch ME
        jrclust.utils.tryClose(hBox);
        errordlg(sprintf('Merge failed:\n\n%s', ME.message), 'Merge auto (checked)');
        return;
    end
    jrclust.utils.tryClose(hBox);

    helpdlg({'Units merged and written to the _res.mat.', '', ...
             'Reopen the session to see them:', '', '    jrc manual <this .prm>', '', ...
             ['Then run Edit > Recompute all mean waveforms - a merged ' ...
              'unit''s stored waveform and quality scores no longer ' ...
              'describe it, so they were set to NaN rather than left ' ...
              'wrong.']}, 'Merge auto (checked)');
end

function relabelUnits(obj)
    %RELABELUNITS Edit > Relabel units.
    %
    %   Re-judges every unit on its spike train and rewrites the annotations
    %   in place, so the verdicts can be refreshed after merging or splitting
    %   without leaving the GUI. Nothing is deleted: only clusterNotes
    %   changes, and the previous 'to_delete' / 'review' labels are cleared
    %   first so that loosening a threshold gives units back rather than only
    %   ever condemning more.
    %
    %   The two thresholds are asked for rather than assumed, because the
    %   right values depend on the recording. Both defaults were calibrated
    %   against units judged by eye rather than derived: the zero-lag cut in
    %   particular was 1.0 on theoretical grounds and deleted a 197,474-spike
    %   unit firing steadily at 21 Hz, because clean units reach 2.05 at the
    %   99th percentile while contaminated ones start at 3.19.
    if exist('jrc_relabel', 'file') ~= 2
        errordlg({'Relabel needs the helper scripts on the MATLAB path.', ...
                  '', 'Run this once, then reopen the GUI:', '', ...
                  '    addpath(''D:\JRCLUST_CPP\matlab'')'}, 'Relabel units');
        return;
    end

    prmFile = obj.hCfg.configFile;
    if isempty(prmFile) || ~isfile(prmFile)
        errordlg('Cannot find this session''s .prm file.', 'Relabel units');
        return;
    end

    dlgAns = inputdlg({['Refractory violation rate above which a unit is ' ...
                        'condemned (fraction of intervals under 1 ms):'], ...
                       ['Autocorrelogram zero-lag ratio above which a unit ' ...
                        'is condemned (1.0 = flat ACG, 3.0 = calibrated):']}, ...
                      'Relabel units', 1, {'0.30', '3.00'});
    if isempty(dlgAns)
        return;
    end

    maxISI = str2double(dlgAns{1});
    maxACG = str2double(dlgAns{2});
    if ~isfinite(maxISI) || maxISI <= 0 || ~isfinite(maxACG) || maxACG <= 0
        jrclust.utils.qMsgBox('Invalid thresholds.');
        return;
    end

    % jrc_relabel reads the _res.mat from disk, so curation held only in the
    % GUI has to be written out first or it would be judged against a stale
    % clustering - and then overwritten by the labels it produces.
    obj.saveFiles();

    hBox = jrclust.utils.qMsgBox('Relabelling... (this closes automatically)', 0, 1);
    try
        jrc_relabel(prmFile, 'maxISIRatio', maxISI, 'maxZeroLag', maxACG);
    catch ME
        jrclust.utils.tryClose(hBox);
        errordlg(sprintf('Relabelling failed:\n\n%s', ME.message), 'Relabel units');
        return;
    end
    jrclust.utils.tryClose(hBox);

    helpdlg({'Labels rewritten to the _res.mat.', '', ...
             'Reopen the session to see them:', '', ...
             '    jrc manual <this .prm>', '', ...
             ['The GUI holds the clustering in memory, so it cannot pick ' ...
              'up the new annotations without reloading.']}, 'Relabel units');
end

function saveAsPruned(obj)
    %SAVEASPRUNED File > Save As Pruned.
    %
    %   Writes a reduced copy of the session with the discarded spikes and
    %   clusters removed from everything - the _res.mat and the _raw, _filt
    %   and _features .jrc binaries - so the pruned dataset is genuinely
    %   smaller rather than merely relabelled.
    %
    %   Edit > Delete annotated cannot do this: it only rewrites the spike
    %   table, leaving every discarded spike in the .jrc files, so nothing
    %   downstream shrinks.
    %
    %   The work is done by jrc_save_as_pruned / jrc_prune, which stream the
    %   .jrc files in fixed-size blocks. Those live outside JRCLUST, so the
    %   menu item checks for them and says what to do rather than throwing an
    %   "Undefined function" dialog at someone who has not added them to the
    %   path.
    if exist('jrc_save_as_pruned', 'file') ~= 2
        errordlg({'Save As Pruned needs the helper scripts on the MATLAB path.', ...
                  '', 'Run this once, then reopen the GUI:', '', ...
                  '    addpath(''D:\JRCLUST_CPP\matlab'')'}, 'Save As Pruned');
        return;
    end

    jrc_save_as_pruned(obj);
end