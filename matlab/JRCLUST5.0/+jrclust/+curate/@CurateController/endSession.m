function res = endSession(obj)
    %ENDSESSION Finish curating and return results
    % isEnding: 0 = curating, 1 = save question is open, 2 = closing windows
    res = [];
    if obj.isEnding
        return; % a save question is already open (Exit pressed twice), or we are closing
    end

    if obj.isWorking
        dlgAns = questdlg('An operation is in progress. Really quit?', 'Operation in progress', 'No');
        if ~strcmp(dlgAns, 'Yes')
            return;
        end
    end

    obj.isEnding = 1;
    if nargout == 0
        success = obj.saveFiles();
        if ~success
            obj.isEnding = 0; % cancelled
            return;
        end
    end
    obj.isEnding = 2;
    obj.closeFigures();

    obj.currentSite = [];
    obj.selected = [];
    res = obj.cRes;
    obj.cRes = [];
    obj.isEnding = 0; % ended
end
