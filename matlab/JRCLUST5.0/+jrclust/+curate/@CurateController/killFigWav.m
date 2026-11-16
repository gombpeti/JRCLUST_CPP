function killFigWav(obj, hFig, ~)
    %KILLFIGWAV Destroy the main figure, close all other figures
    if obj.isEnding == 1
        return; % the save question is open: answer it first, don't kill the window under it
    end

    if obj.isEnding == 0
        obj.endSession(); % asks to save, then closes everything (calls back here with isEnding == 2)
    end

    if obj.isEnding == 2
        delete(hFig);
    end
end
