function success = saveFiles(obj)
    %SAVEFILES Save the curation to the _res.mat
    %   No question box and no "Saving..." box: both need the GUI to draw, and
    %   when MATLAB's figure drawing stalls they never answer, so the save
    %   never happened. Progress goes to the Command Window instead.
    %   The file is written under a temporary name first and then renamed, so
    %   an interrupted save never leaves a half-written _res.mat behind.
    success = 0;
    resFile = obj.hCfg.resFile;
    tmpFile = regexprep(resFile, '\.mat$', '_saving.mat');
    if strcmp(tmpFile, resFile)
        tmpFile = [resFile, '_saving.mat'];
    end

    fprintf('JRCLUST: saving clustering to %s ...\n', resFile);
    t = tic;
    try
        obj.cRes.curatedOn = now();

        res_ = jrclust.utils.mergeStructs(obj.res, obj.cRes);
        res_ = rmfield(res_, 'hClust');

        hClust = obj.hClust;
        res_.initialClustering = hClust.initialClustering;
        res_.spikeClusters = hClust.spikeClusters;

        % fieldnames contained in dRes or sRes
        dsFields = union(fieldnames(hClust.dRes), fieldnames(hClust.sRes));
        % fieldnames from hClust which are not in dRes or sRes
        hClustOnly = setdiff(fieldnames(hClust), dsFields);
        for i = 1:numel(hClustOnly)
            fn = hClustOnly{i};
            % hClust fields take precedence
            res_.(fn) = hClust.(fn);
        end

        % don't save these fields
        dropFields = {'spikesRaw', 'spikesRawVolt', 'spikesFilt', 'spikesFiltVolt', ...
                      'spikesFilt2', 'spikesFilt3', 'spikeFeatures', 'nClusters', ...
                      'nEdits', 'unitFields', 'hRecs'};
        res_ = rmfield(res_, intersect(fieldnames(res_), dropFields));

        if exist(tmpFile, 'file') == 2
            delete(tmpFile);
        end
        save(tmpFile, '-struct', 'res_', '-v7.3', '-nocompression');
        [ok, msg] = movefile(tmpFile, resFile, 'f');
        if ~ok
            error('could not replace %s: %s (the new result is in %s)', resFile, msg, tmpFile);
        end

        success = 1;
        d = dir(resFile);
        fprintf('JRCLUST: saved %d units in %.0f s (%.0f MB, %s)\n', ...
                hClust.nClusters, toc(t), d.bytes/1e6, datestr(d.datenum));
    catch ME
        fprintf(2, 'JRCLUST: SAVE FAILED: %s\n', ME.message);
        fprintf(2, 'JRCLUST: nothing was lost - your edits are still open. Try File > Save again, or run jrc_save_curation in the Command Window.\n');
    end
end
