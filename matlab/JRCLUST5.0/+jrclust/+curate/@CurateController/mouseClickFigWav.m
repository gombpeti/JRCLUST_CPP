function mouseClickFigWav(obj, xyPos, clickType)
    %MOUSECLICKFIGWAV Handle callbacks for mouse clicks in main view
    if obj.isWorking
        return;
    end

    iCluster = obj.showSubset(max(1, min(obj.nShown, round(xyPos(1))))); % floor of x position
    if iCluster < 1 || iCluster > obj.hClust.nClusters
        return;
    end

    if strcmp(clickType, 'normal')
        % plain click starts a fresh selection
        obj.updateSelect(iCluster);
    elseif strcmp(clickType, 'alt') || strcmp(clickType, 'extend')
        % ctrl-click, right click or shift-click ADDS the unit to the
        % selection, and clicking a selected unit again removes it. MATLAB
        % reports ctrl-click and right click both as 'alt', and shift-click
        % as 'extend', so all three accumulate.
        %
        % This used to replace the second slot - [obj.selected(1) iCluster] -
        % which made a selection of three impossible. Accumulating lets
        % FigCorr show every pairwise correlogram at once, which is the
        % point: deciding whether four units are one cell needs all six
        % cross-correlograms, not three separate pairwise comparisons.
        sel = obj.selected;
        if any(sel == iCluster)
            sel(sel == iCluster) = [];      % toggle off
            if isempty(sel)
                sel = iCluster;             % never leave nothing selected
            end
        else
            sel = [sel, iCluster];
        end
        obj.updateSelect(sel);
    else                            % middle click, ignore
        disp(clickType);
        return;
    end

    hFig = obj.hFigs('FigWav');
    hFig.wait(1);
    hFig.wait(0);
end