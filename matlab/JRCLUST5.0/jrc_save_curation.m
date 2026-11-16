function jrc_save_curation(action)
%JRC_SAVE_CURATION Save the open JRCLUST curation from the Command Window
%   jrc_save_curation          save the open curation to its _res.mat
%   jrc_save_curation close    save, then close all curation windows
%   jrc_save_curation discard  close all curation windows WITHOUT saving
%
%   For when the curation windows stop responding. Press Ctrl+C in the
%   Command Window until the >> prompt is back, then run this. It draws
%   nothing and asks nothing, so it works even when MATLAB's figure drawing
%   is stuck. Edits made in the GUI are all in memory until saved.
    if nargin < 1
        action = 'save';
    end
    action = lower(action);
    if ~ismember(action, {'save', 'close', 'discard'})
        error('Use: jrc_save_curation, jrc_save_curation close, or jrc_save_curation discard');
    end

    hC = findCurateController();
    if isempty(hC)
        error('No open JRCLUST curation found (is the Averaged waveform window still open?).');
    end
    hC.isBusy = 0; % an interrupted action may have left this set

    if ~strcmp(action, 'discard')
        if ~hC.saveFiles()
            error('Save failed - the curation is still open, nothing was lost.');
        end
    end
    if ismember(action, {'close', 'discard'})
        hC.closeFigures();
        fprintf('JRCLUST: curation windows closed.\n');
    end
end

function hC = findCurateController()
    %FINDCURATECONTROLLER The controller is captured by the windows' key handlers
    hC = [];
    figs = findall(groot, 'Type', 'figure');
    props = {'KeyPressFcn', 'CloseRequestFcn', 'WindowButtonDownFcn'};
    for i = 1:numel(figs)
        for j = 1:numel(props)
            try
                fh = get(figs(i), props{j});
            catch
                continue;
            end
            hC = controllerIn(fh);
            if ~isempty(hC)
                return;
            end
        end
    end
end

function hC = controllerIn(fh)
    hC = [];
    if ~isa(fh, 'function_handle')
        return;
    end
    info = functions(fh);
    if isfield(info, 'workspace') && ~isempty(info.workspace)
        ws = info.workspace{1};
        vars = fieldnames(ws);
        for k = 1:numel(vars)
            v = ws.(vars{k});
            if isa(v, 'jrclust.curate.CurateController')
                hC = v;
                return;
            end
        end
    end
end
