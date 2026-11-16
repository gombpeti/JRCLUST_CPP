function res = endSession(obj, doSave)
    %ENDSESSION Finish curating and return results
    %   endSession()     - File > Exit / closing the main window: save, then close
    %   endSession(0)    - File > Exit without saving
    %   res = endSession - called by code: return the results, no save
    %   No question boxes: they need the GUI to draw, and when MATLAB's figure
    %   drawing stalls they never answered, so Exit could not finish.
    %   isEnding: 0 = curating, 1 = saving, 2 = closing windows
    res = [];
    if nargin < 2
        doSave = (nargout == 0);
    end
    if obj.isEnding
        return; % already saving or closing
    end

    obj.isEnding = 1;
    if doSave
        if ~obj.saveFiles()
            obj.isEnding = 0; % save failed: stay open so nothing is lost
            fprintf(2, 'JRCLUST: not closing, because the save failed.\n');
            return;
        end
    else
        fprintf('JRCLUST: closing without saving.\n');
    end

    obj.isEnding = 2;
    obj.closeFigures();

    obj.currentSite = [];
    obj.selected = [];
    res = obj.cRes;
    obj.cRes = [];
    obj.isEnding = 0; % ended
end
