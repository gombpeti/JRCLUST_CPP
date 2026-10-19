function jrc_save_as_pruned(obj)
%JRC_SAVE_AS_PRUNED File > Save As Pruned, for the curation GUI.
%
%   Restores the menu item JRCLUST used to carry: take the curated result,
%   physically remove the discarded spikes and clusters from EVERYTHING - the
%   _res.mat and the _raw / _filt / _features .jrc binaries - and write the
%   reduced dataset to a 'pruned' folder beside the .prm, with a .prm of its
%   own so it opens like any other session.
%
%   This build of JRCLUST has no such menu item; the only bulk removal is
%   Edit > Delete annotated, which cannot do the job. That path calls
%   deleteUnit once per unit, and deleteUnit snapshots EVERY per-cluster
%   field each time - three mean-waveform arrays and the full nClusters x
%   nClusters similarity matrix among them. On a 1143-unit recording it
%   reached 73 GB and died partway through, leaving the clustering
%   inconsistent. It also only ever edits the spike table: the .jrc files
%   keep every discarded spike, so nothing downstream gets smaller.
%
%   jrc_prune streams the .jrc files in fixed-size blocks instead, so its
%   memory use does not depend on how much is being removed.
%
%   Units annotated 'to_delete' are the ones removed. To include units
%   annotated 'noise', run jrc_mark_noise first.
%
%   Called from the GUI; takes the CurateController. Not meant to be invoked
%   directly - use jrc_prune(prmFile) for that.

    if nargin < 1 || isempty(obj)
        error('jrc_save_as_pruned:noController', ...
              'call jrc_prune(prmFile) directly outside the GUI');
    end

    prmFile = obj.hCfg.configFile;
    if isempty(prmFile) || ~isfile(prmFile)
        errordlg('Cannot find this session''s .prm file.', 'Save As Pruned');
        return;
    end

    % Count what will go, so the question is answerable rather than blind.
    notes = {};
    try
        notes = obj.hClust.clusterNotes(:);
    catch
    end
    nDelete = sum(strcmp(notes, 'to_delete'));
    nNoise  = sum(strcmp(notes, 'noise'));

    msg = sprintf(['Remove %d unit(s) annotated ''to_delete'' from the ' ...
                   'spike table AND the .jrc files, writing a reduced copy ' ...
                   'to\n\n    %s\n\n'], nDelete, fullfile(fileparts(prmFile), 'pruned'));
    if nNoise > 0
        msg = [msg sprintf(['%d unit(s) annotated ''noise'' will be KEPT. ' ...
                            'Run jrc_mark_noise first to include them.\n\n'], nNoise)];
    end
    msg = [msg 'Unsaved curation will be saved first. Continue?'];

    if nDelete == 0
        warndlg('No units are annotated ''to_delete''. Nothing to prune.', ...
                'Save As Pruned');
        return;
    end

    if ~strcmp(questdlg(msg, 'Save As Pruned', 'Yes', 'No', 'No'), 'Yes')
        return;
    end

    % The prune reads the _res.mat from disk, so curation still held only in
    % the GUI has to be written out first or it would be silently dropped.
    obj.saveFiles();

    try
        jrc_prune(prmFile);
    catch ME
        errordlg(sprintf('Pruning failed:\n\n%s', ME.message), 'Save As Pruned');
        rethrow(ME);
    end

    [prmDir, sessionName] = fileparts(prmFile);
    prunedPrm = fullfile(prmDir, 'pruned', [sessionName '_pruned.prm']);
    helpdlg(sprintf(['Pruned dataset written.\n\nOpen it with:\n\n' ...
                     '    jrc manual %s'], prunedPrm), 'Save As Pruned');
end
