function killFigWav(obj, hFig, ~)
    %KILLFIGWAV Closing the main window: save, then close everything
    if obj.isEnding == 2
        delete(hFig); % endSession is closing the windows
        return;
    end
    if obj.isEnding == 1
        return; % a save is running; it closes the windows when it is done
    end

    if obj.isBusy
        fprintf(2, 'JRCLUST is still busy - close the window again when it has finished.\n');
        return;
    end
    obj.runExclusive(@() obj.endSession(1));
end
