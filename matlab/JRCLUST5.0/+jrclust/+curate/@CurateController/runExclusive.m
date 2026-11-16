function runExclusive(obj, fun, varargin)
    %RUNEXCLUSIVE Run a GUI callback, unless another one is still running
    %   MATLAB runs any key press, mouse click or menu click that is waiting
    %   in the queue INSIDE the drawnow/pause/uiwait of the callback that is
    %   still running. A split could then start in the middle of the redraw
    %   after a delete, a Save in the middle of the split, and so on - one
    %   nested in the next. Every key press, click and menu item of the
    %   curation GUI goes through here, so a second one is ignored (with a
    %   note in the Command Window) until the first has finished.
    %   The flag is cleared even on an error or Ctrl+C (onCleanup).
    if obj.isBusy
        fprintf(2, 'JRCLUST is still busy with the previous action - ignored. (Ctrl+C in the Command Window clears it.)\n');
        return;
    end

    obj.isBusy = 1;
    cleanup = onCleanup(@() clearBusy(obj)); %#ok<NASGU>
    fun(varargin{:});
end

function clearBusy(obj)
    try
        obj.isBusy = 0;
    catch
    end
end
