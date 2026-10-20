classdef TracesController < jrclust.interfaces.FigureController
    %TRACESCONTROLLER Summary of this class goes here
    properties (SetAccess=private)
        hCfg;
        hClust;
        hRec;
    end

    properties (SetAccess=private, Transient)
        hFigPSD;
        hFigTraces;

        showFilter;
        showGrid;
        showSpikes;
        showTraces;

        tracesFull;
        tracesRaw;
        tracesFilt;

        shankOrder;     % permutation of unique(hCfg.shankMap): display order of ALL shanks (visible or hidden)
        shankVisible;   % logical, aligned with shankOrder: true if currently shown as a panel
        shankUI;        % struct array (one per visible panel): dropdown/hideBtn uicontrol handles

        durationEdit;   % uicontrol: window-width edit box
        hHiddenPopup;   % uicontrol: popup listing currently-hidden shanks
        hShowBtn;       % uicontrol: "Show" button for the hidden-shank popup
    end

    %% LIFECYCLE
    methods
        function obj = TracesController(hCfg)
            %TRACESCONTROLLER Construct an instance of this class
            obj.hCfg = hCfg;
        end
    end

    %% KEYPRESS/MOUSECLICK METHODS
    methods
        function keyPressFigTraces(obj, ~, hEvent)
            %KEYPRESSFIGTRACES
            factor = 4^double(jrclust.utils.keyMod(hEvent, 'shift')); % 1 or 4

            switch hEvent.Key
                case 'uparrow'
                    obj.hFigTraces.figData.maxAmp = obj.hFigTraces.figData.maxAmp*sqrt(2)^-factor;
                    obj.updateFigTraces(0);

                case 'downarrow'
                    obj.hFigTraces.figData.maxAmp = obj.hFigTraces.figData.maxAmp*sqrt(2)^factor;
                    obj.updateFigTraces(0);

                case {'leftarrow', 'rightarrow', 'j', 'home', 'end'}
                    switch lower(hEvent.Key)
                        case 'leftarrow'
                            if obj.hFigTraces.figData.windowBounds(1) == 1
                                jrclust.utils.qMsgBox('Beginning of file', 1);
                                return;
                            end

                            windowBounds = obj.hFigTraces.figData.windowBounds - (obj.hFigTraces.figData.windowWidth) * factor; %no overlap
                            if windowBounds(1) < 1
                                windowBounds = [1, obj.hFigTraces.figData.windowWidth];
                            end

                        case 'rightarrow'
                            if obj.hFigTraces.figData.windowBounds(2) == obj.hFigTraces.figData.nSamplesTotal
                                jrclust.utils.qMsgBox('End of file', 1);
                                return;
                            end

                            windowBounds = obj.hFigTraces.figData.windowBounds + (obj.hFigTraces.figData.windowWidth + 1) * factor; %no overlap
                            if windowBounds(2) > obj.hFigTraces.figData.nSamplesTotal
                                windowBounds = [-obj.hFigTraces.figData.windowWidth + 1, 0] + obj.hFigTraces.figData.nSamplesTotal;
                            end

                        case 'home' % beginning of file
                            windowBounds = [1, obj.hFigTraces.figData.windowWidth];

                        case 'end' % end of file
                            windowBounds = [-obj.hFigTraces.figData.windowWidth+1, 0] + obj.hFigTraces.figData.nSamplesTotal;

                        case 'j'
                            dlgAns = inputdlg('Go to time (s)', 'Jump to time', 1, {'0'});

                            if isempty(dlgAns)
                                return;
                            end

                            try
                                windowBounds = round(str2double(dlgAns)*obj.hCfg.sampleRate) + [1, obj.hFigTraces.figData.windowWidth];
                            catch
                                return;
                            end
                    end % switch

                    nTimeTraces = obj.hCfg.nSegmentsTraces;
                    multiBounds = jrclust.views.sampleSkip(windowBounds, obj.hFigTraces.figData.nSamplesTotal, nTimeTraces);

                    tracesRaw_ = cellfun(@(lims) obj.hRec.readRawROI(obj.hCfg.siteMap, lims(1):lims(2)), multiBounds, 'UniformOutput', 0);
                    obj.tracesRaw = jrclust.utils.neCell2mat(tracesRaw_);

                    obj.tracesRaw = u2i(obj.tracesRaw);
                    obj.hFigTraces.figData.windowBounds = windowBounds;
                    obj.updateFigTraces(1);

                case 'w' % set window width
                    curWidthSec = obj.hFigTraces.figData.windowWidth / obj.hCfg.sampleRate;
                    dlgAns = inputdlg('Window width (s)', 'Set window width', 1, {sprintf('%0.3g', curWidthSec)});

                    if isempty(dlgAns)
                        return;
                    end

                    obj.setWindowWidthSec(str2double(dlgAns{1}));

                case 'c' % channel query (operates on the first visible shank panel)
                    visibleFullIdx = find(obj.shankVisible);
                    if isempty(visibleFullIdx)
                        return;
                    end

                    hAx1 = obj.hFigTraces.figData.shankAxes(1);
                    shankID1 = obj.shankOrder(visibleFullIdx(1));
                    siteIdx1 = find(obj.hCfg.shankMap == shankID1);
                    plotKey1 = sprintf('hTrace_%d', shankID1);

                    jrclust.utils.qMsgBox('Draw a rectangle', 1);
                    obj.hFigTraces.addPlot('hRect', @imrect, hAx1);
                    ie = obj.hFigTraces.plotApply('hRect', @isempty);
                    if ie
                        return;
                    end

                    rectPos = obj.hFigTraces.plotApply('hRect', @getPosition);
                    UserData = obj.hFigTraces.plotApply(plotKey1, @get, 'UserData');

                    XData = obj.hFigTraces.plotApply(plotKey1, @get, 'XData');
                    YData = obj.hFigTraces.plotApply(plotKey1, @get, 'YData');
                    inBounds = find(XData >= rectPos(1) & XData <= sum(rectPos([1, 3])) & YData >= rectPos(2) & YData <= sum(rectPos([2, 4])));

                    if isempty(inBounds)
                        obj.hFigTraces.rmPlot('hRect');
                        return;
                    end

                    anchorPoint = round(median(inBounds));

                    [~, iLocal] = ind2sub(UserData.shape, anchorPoint);
                    iSite = siteIdx1(iLocal);

                    mrX = reshape(XData, UserData.shape);
                    mrY = reshape(YData, UserData.shape);

                    hold(hAx1, 'on');

                    obj.hFigTraces.addPlot('hPoint', hAx1, XData(anchorPoint), YData(anchorPoint), 'r*');
                    obj.hFigTraces.addPlot('hLine', hAx1, mrX(:, iLocal), mrY(:, iLocal), 'r-');

                    hold(hAx1, 'off');

                    iChan = obj.hCfg.siteMap(iSite);
                    jrclust.utils.qMsgBox(sprintf('Site: %d/ Chan: %d', iSite, iChan), 1);

                    % clean up
                    obj.hFigTraces.rmPlot('hRect');
                    obj.hFigTraces.rmPlot('hLine');
                    obj.hFigTraces.rmPlot('hPoint');

                case 'e' % export current view
                    jrclust.utils.exportToWorkspace(struct('tracesRaw', obj.tracesRaw, ...
                                                           'tracesFilt', obj.tracesFilt), 1);

                case 'f' % toggle filter
                    obj.toggleFilter();

                case 'g' % toggle grid
                    obj.toggleGrid();

                case 'h'
                    jrclust.utils.qMsgBox(obj.hFigTraces.figData.helpText, 1);

                case 'p' % power spectrum
                    iSite = jrclust.utils.inputdlgNum(sprintf('Site# to show (1-%d, 0 for all)', obj.hCfg.nSites), 'Site#', 0);

                    if isnan(iSite)
                        return;
                    end

                    obj.hFigPSD = jrclust.views.Figure('FigPsd', [.5 0 .5 1], obj.hCfg.configFile, 1, 1); % show to the right

                    % ask user which channels to plot
                    if iSite > 0
                        tracesFilt_ = obj.tracesFilt(:, iSite)';
                    else
                        tracesFilt_ = obj.tracesFilt;
                    end

                    jrclust.views.plotFigPSD(obj.hFigPSD, tracesFilt_, obj.hCfg);

                case 'r' % reset view
                    if obj.hFigTraces.figData.maxAmp ~= obj.hCfg.maxAmp
                        obj.hFigTraces.figData.maxAmp = obj.hCfg.maxAmp;
                        obj.updateFigTraces(1);
                    else
                        jrclust.views.resetFigTraces(obj.hFigTraces, obj.tracesRaw, obj.hCfg);
                    end

                case 's' %show/hide spikes
                    obj.toggleSpikes();

                case 't' %show/hide traces
                    obj.toggleTraces();
            end % switch
        end
    end

    %% UTILITY METHODS
    methods (Hidden)
        function closeFigTraces(obj, hFig, hEvent)
            obj.hRec.closeRaw();
            delete(hFig);
        end

        function toggleFilter(obj)
            obj.showFilter = ~obj.showFilter;
            obj.hFigTraces.figData.filter = jrclust.utils.ifEq(obj.showFilter, 'on', 'off');
            obj.updateFigTraces(0);
        end

        function toggleGrid(obj)
            obj.showGrid = ~obj.showGrid;
            obj.hFigTraces.figData.grid = jrclust.utils.ifEq(obj.showGrid, 'on', 'off');
            arrayfun(@(hAx) grid(hAx, obj.hFigTraces.figData.grid), obj.hFigTraces.figData.shankAxes);
        end

        function toggleSpikes(obj)
            obj.showSpikes = ~obj.showSpikes;
            obj.hFigTraces.figData.spikes = jrclust.utils.ifEq(obj.showSpikes, 'on', 'off');
            obj.updateFigTraces(0);
        end

        function toggleTraces(obj)
            obj.showTraces = ~obj.showTraces;
            obj.hFigTraces.figData.traces = jrclust.utils.ifEq(obj.showTraces, 'on', 'off');
            obj.updateFigTraces(0);
        end

        function updateFigTraces(obj, resetAxes)
            obj.tracesFilt = jrclust.views.plotFigTraces(obj.hFigTraces, obj.hCfg, obj.tracesRaw, resetAxes, obj.hClust);
        end

        function setWindowWidthSec(obj, newWidthSec)
            %SETWINDOWWIDTHSEC Change the displayed time-window width (in seconds) and reload/redraw
            if isnan(newWidthSec) || newWidthSec <= 0
                return;
            end

            newWidth = min(round(newWidthSec * obj.hCfg.sampleRate), obj.hFigTraces.figData.nSamplesTotal);
            obj.hFigTraces.figData.windowWidth = newWidth;

            windowBounds = [obj.hFigTraces.figData.windowBounds(1), obj.hFigTraces.figData.windowBounds(1) + newWidth - 1];
            if windowBounds(2) > obj.hFigTraces.figData.nSamplesTotal
                windowBounds = [-newWidth + 1, 0] + obj.hFigTraces.figData.nSamplesTotal;
            end
            if windowBounds(1) < 1
                windowBounds = [1, newWidth];
            end

            nTimeTraces = obj.hCfg.nSegmentsTraces;
            multiBounds = jrclust.views.sampleSkip(windowBounds, obj.hFigTraces.figData.nSamplesTotal, nTimeTraces);

            tracesRaw_ = cellfun(@(lims) obj.hRec.readRawROI(obj.hCfg.siteMap, lims(1):lims(2)), multiBounds, 'UniformOutput', 0);
            obj.tracesRaw = jrclust.utils.neCell2mat(tracesRaw_);
            obj.tracesRaw = u2i(obj.tracesRaw);
            obj.hFigTraces.figData.windowBounds = windowBounds;

            if ~isempty(obj.durationEdit) && isgraphics(obj.durationEdit)
                set(obj.durationEdit, 'String', sprintf('%0.3g', newWidthSec));
            end

            obj.updateFigTraces(1);
        end

        function onDurationApply(obj, ~, ~)
            %ONDURATIONAPPLY Callback for the persistent window-width edit box/Apply button
            newWidthSec = str2double(get(obj.durationEdit, 'String'));
            obj.setWindowWidthSec(newWidthSec);
        end

        function onSlotShankChanged(obj, slotFullIdx, hDropdown)
            %ONSLOTSHANKCHANGED Callback: user picked a different shank for a visible panel (swaps positions)
            newFullIdx = get(hDropdown, 'Value');
            if newFullIdx == slotFullIdx
                return;
            end

            [obj.shankOrder(slotFullIdx), obj.shankOrder(newFullIdx)] = deal(obj.shankOrder(newFullIdx), obj.shankOrder(slotFullIdx));
            [obj.shankVisible(slotFullIdx), obj.shankVisible(newFullIdx)] = deal(obj.shankVisible(newFullIdx), obj.shankVisible(slotFullIdx));

            obj.rebuildShankPanels();
            obj.updateFigTraces(1);
        end

        function onHideSlot(obj, slotFullIdx)
            %ONHIDESLOT Callback: user hid a visible shank panel
            obj.shankVisible(slotFullIdx) = false;
            obj.rebuildShankPanels();
            obj.updateFigTraces(1);
        end

        function onShowHiddenShank(obj, ~, ~)
            %ONSHOWHIDDENSHANK Callback: bring a hidden shank back into view
            hiddenFullIdx = find(~obj.shankVisible);
            if isempty(hiddenFullIdx)
                return;
            end

            chosen = get(obj.hHiddenPopup, 'Value');
            chosen = min(chosen, numel(hiddenFullIdx));
            obj.shankVisible(hiddenFullIdx(chosen)) = true;

            obj.rebuildShankPanels();
            obj.updateFigTraces(1);
        end

        function updateHiddenPopup(obj)
            %UPDATEHIDDENPOPUP Refresh the "hidden shanks" popup + Show button state
            hiddenFullIdx = find(~obj.shankVisible);
            if isempty(hiddenFullIdx)
                set(obj.hHiddenPopup, 'String', {'(none hidden)'}, 'Value', 1, 'Enable', 'off');
                set(obj.hShowBtn, 'Enable', 'off');
            else
                labels = arrayfun(@(s) sprintf('Shank %d', s), obj.shankOrder(hiddenFullIdx), 'UniformOutput', 0);
                set(obj.hHiddenPopup, 'String', labels, 'Value', 1, 'Enable', 'on');
                set(obj.hShowBtn, 'Enable', 'on');
            end
        end

        function rebuildShankPanels(obj)
            %REBUILDSHANKPANELS (Re)create the row of per-shank subplot panels and their per-panel controls
            % let MATLAB finish dispatching the callback that triggered this (which may be a
            % button we're about to delete below) before we start destroying uicontrols
            drawnow;

            % keep figData in sync: plotFigTraces/resetFigTraces read shankOrder/shankVisible
            % from figData (they're plain functions, not methods, so they can't see obj.* directly)
            obj.hFigTraces.figData.shankOrder = obj.shankOrder;
            obj.hFigTraces.figData.shankVisible = obj.shankVisible;

            % clean up old panel axes
            if isKey(obj.hFigTraces.hSubplots, 'hShankPanels')
                oldAxes = obj.hFigTraces.hSubplots('hShankPanels');
                delete(oldAxes(isgraphics(oldAxes)));
                remove(obj.hFigTraces.hSubplots, 'hShankPanels');
            end

            % clean up old per-panel uicontrols
            for i = 1:numel(obj.shankUI)
                if isgraphics(obj.shankUI(i).dropdown)
                    delete(obj.shankUI(i).dropdown);
                end
                if isgraphics(obj.shankUI(i).hideBtn)
                    delete(obj.shankUI(i).hideBtn);
                end
            end
            obj.shankUI = struct('dropdown', {}, 'hideBtn', {});

            visibleFullIdx = find(obj.shankVisible);
            nVisible = numel(visibleFullIdx);

            if nVisible == 0
                obj.hFigTraces.figData.shankAxes = matlab.graphics.axis.Axes.empty;
                obj.updateHiddenPopup();
                return;
            end

            obj.hFigTraces.addSubplot('hShankPanels', 1, nVisible);
            axesArr = obj.hFigTraces.hSubplots('hShankPanels');

            marginL = 0.05; marginR = 0.01; gap = 0.012;
            axTop = 0.88; axBottom = 0.06;
            panelW = (1 - marginL - marginR - gap*(nVisible-1)) / nVisible;

            for v = 1:nVisible
                hAx = axesArr(v);
                xLeft = marginL + (v-1)*(panelW+gap);

                set(hAx, 'Position', [xLeft, axBottom, panelW, axTop-axBottom], ...
                         'Color', 'k', 'XColor', [1 1 1]*.85, 'YColor', [1 1 1]*.85, ...
                         'GridColor', [1 1 1]*.35, 'MinorGridColor', [1 1 1]*.2, ...
                         'GridAlpha', .35, 'MinorGridAlpha', .25, 'Box', 'on', ...
                         'XLimMode', 'manual', 'YLimMode', 'manual');
                hold(hAx, 'on');

                shankLabels = arrayfun(@(s) sprintf('Shank %d', s), obj.shankOrder, 'UniformOutput', 0);
                slotFullIdx = visibleFullIdx(v);

                hDropdown = uicontrol(obj.hFigTraces.hFig, 'Style', 'popupmenu', 'Units', 'normalized', ...
                    'Position', [xLeft, axTop+0.015, panelW*0.65, 0.035], ...
                    'String', shankLabels, 'Value', slotFullIdx, ...
                    'Callback', @(src, ~) obj.onSlotShankChanged(slotFullIdx, src));

                hHideBtn = uicontrol(obj.hFigTraces.hFig, 'Style', 'pushbutton', 'Units', 'normalized', ...
                    'Position', [xLeft+panelW*0.68, axTop+0.015, panelW*0.3, 0.035], ...
                    'String', 'Hide', ...
                    'Callback', @(~, ~) obj.onHideSlot(slotFullIdx));

                obj.shankUI(v).dropdown = hDropdown;
                obj.shankUI(v).hideBtn = hHideBtn;
            end

            obj.hFigTraces.figData.shankAxes = axesArr;
            obj.updateHiddenPopup();
        end
    end

    %% USER METHODS
    methods
        function show(obj, recID, showLFP, hClust)
            if nargin < 4
                hClust = [];
            end
            if nargin < 3
                showLFP = 0;
            end
            if nargin < 2
                recID = 1;
            end

            obj.hClust = hClust;

            % get file to show
            if numel(obj.hCfg.rawRecordings) > 1
                if isempty(recID)
                    arrayfun(@(i) fprintf('%d: %s\n', i, obj.hCfg.rawRecordings{i}), 1:numel(obj.hCfg.rawRecordings), 'UniformOutput', 0);
                    fprintf('---------------------------------------------\n');
                    recID = str2double(input('Please specify File ID from the list above: ', 's'));
                end

                if isnan(recID)
                    return;
                end

                try
                    recFilename = obj.hCfg.rawRecordings{recID};
                catch
                    return;
                end
            else
                recFilename = obj.hCfg.rawRecordings{1}; % if multiple files exist, load first
            end

            fprintf('Opening %s\n', recFilename);
            obj.hRec = jrclust.detect.newRecording(recFilename, obj.hCfg);

        %     [fid_bin, nBytes_bin] = fopen_(vcFile_bin, 'r');
        %     if isempty(fid_bin)
        %         fprintf(2, '.bin file does not exist: %s\n', vcFile_bin);
        %         return;
        %     end
        %     nSamplesTotal = floor(nBytes_bin / bytesPerSample_(obj.hCfg.vcDataType) / obj.hCfg.nChans);
            nSamplesTotal = obj.hRec.nSamples;
            windowWidth = min(floor(diff(obj.hCfg.dispTimeLimits) * obj.hCfg.sampleRate), nSamplesTotal);

            if obj.hCfg.dispTimeLimits(1) > 0
                iSampleOffset = ceil(obj.hCfg.dispTimeLimits(1) * obj.hCfg.sampleRate) + 1; % offset sample number
            else
                iSampleOffset = 1; %sample start location
            end

            windowBounds = [0, windowWidth - 1] + iSampleOffset;
            if windowBounds(1) < 1
                windowBounds = [1, windowWidth];
            end
            if windowBounds(2) > nSamplesTotal
                windowBounds = [-windowWidth + 1, 0] + nSamplesTotal;
            end

            multiBounds = jrclust.views.sampleSkip(windowBounds, nSamplesTotal, obj.hCfg.nSegmentsTraces);

            obj.tracesFull = [];
            tracesRaw_ = cellfun(@(lims) obj.hRec.readRawROI(obj.hCfg.siteMap, lims(1):lims(2)), multiBounds, 'UniformOutput', 0);
            obj.tracesRaw = jrclust.utils.neCell2mat(tracesRaw_);

        %     if obj.hCfg.tallSkinny
        %         obj.tracesFull = [];
        %         fseek_(fid_bin, iSampleOffset, obj.hCfg);
        %         if obj.hCfg.nSegmentsTraces > 1
        %             obj.tracesRaw = load_bin_multi_(fid_bin, cvn_lim_bin, obj.hCfg)';
        %         else
        %             obj.tracesRaw = jrclust.utils.readBin(fid_bin, obj.hCfg.vcDataType, [obj.hCfg.nChans, windowWidth])'; %next keypress: update tlim_show
        %         end
        %         %     @TODO: load from cvn_lim_bin specifiers. check for end or beginning when keyboard command
        %     else %load whole thing
        %         obj.tracesFull = jrclust.utils.readBin(fid_bin, obj.hCfg.vcDataType, [nSamplesTotal, obj.hCfg.nChans]); %next keypress: update tlim_show
        %         fclose(fid_bin);
        %         fid_bin = [];
        %         %obj.tracesRaw = obj.tracesFull((windowBounds(1):windowBounds(2)), :);
        %         obj.tracesRaw = obj.tracesFull(viRange_bin, :);
        %         disp('Entire raw traces are cached to RAM since fTranspose=0.');
        %     end %if
            obj.tracesRaw = u2i(obj.tracesRaw);

            % shanks: all visible by default, in natural (ascending) order
            obj.shankOrder = unique(obj.hCfg.shankMap);
            obj.shankVisible = true(size(obj.shankOrder));
            obj.shankUI = struct('dropdown', {}, 'hideBtn', {});

            obj.hFigTraces = jrclust.views.Figure('FigTraces', [0 0 .5 1], recFilename, 0, 1);

            figData = struct('windowBounds', windowBounds, ...
                             'hRec', obj.hRec, ...
                             'nSamplesTotal', nSamplesTotal, ...
                             'windowWidth', windowWidth, ...
                             'shankOrder', obj.shankOrder, ...
                             'shankVisible', obj.shankVisible, ...
                             'shankAxes', matlab.graphics.axis.Axes.empty);

            figData.maxAmp = obj.hCfg.maxAmp;
            figData.traceColor = [1 1 1]*.88;
            figData.textColor = [1 1 1]*.9;
            figData.title = '[H]elp; (Sft)[Up/Down]:Scale(%0.1f uV); (Sft)[Left/Right]:Time; [F]ilter; [J]ump T; [W]indow width; [C]han. query; [R]eset view; [P]SD; [S]pike; [E]xport; [T]race; [G]rid';
            figData.helpText = { ...
                'Left/Right: change time (Shift: x4)', ...
                '[J]ump T', ...
                '[W]indow width (s), or use the box at the top', ...
                '[Home/End]: go to beginning/end of file', ...
                '---------', ...
                'Up/Down: change scale (Shift: x4)', ...
                'Zoom: Mouse wheel', ...
                '[x/y/ESC]: zoom direction', ...
                'Pan: hold down the wheel and drag', ...
                '[R]eset view', ...
                '---------', ...
                '[F]ilter toggle', ...
                '[S]pike toggle', ...
                'Gri[D] toggle', ...
                '[T]races toggle', ...
                '---------', ...
                '[C]hannel query (first visible shank panel only)', ...
                '[P]ower spectrum', ...
                '[E]xport to workspace', ...
                '---------', ...
                'Each shank panel has a dropdown (pick which shank to show there, swaps with its previous slot) and a Hide button', ...
                'Use the "Show hidden shank" control (top strip) to bring a hidden shank back', ...
            }; % TODO: '[A]ux channel display'

            obj.showFilter = 0;
            figData.filter = 'off';

            obj.showGrid = 1;
            figData.grid = 'on';

            obj.showSpikes = obj.hCfg.getOr('showSpikesTraces', 0);
	    figData.spikes = 'off';
            if obj.showSpikes
                figData.spikes = 'on';
            end

            obj.showTraces = 1;
            figData.traces = 'on';

            obj.hFigTraces.figData = figData;

            obj.hFigTraces.figApply(@set, 'Color', 'k', 'BusyAction', 'cancel', 'CloseRequestFcn', @obj.closeFigTraces);
            obj.hFigTraces.hFunKey = @obj.keyPressFigTraces;
            obj.hFigTraces.setMouseable();

            % persistent control strip: window width + hidden-shank recovery
            uicontrol(obj.hFigTraces.hFig, 'Style', 'text', 'Units', 'normalized', ...
                'Position', [0.01 0.95 0.09 0.035], 'String', 'Window (s):', ...
                'BackgroundColor', 'k', 'ForegroundColor', 'w', 'HorizontalAlignment', 'left');
            obj.durationEdit = uicontrol(obj.hFigTraces.hFig, 'Style', 'edit', 'Units', 'normalized', ...
                'Position', [0.10 0.95 0.06 0.035], 'String', sprintf('%0.3g', windowWidth/obj.hCfg.sampleRate), ...
                'Callback', @obj.onDurationApply);
            uicontrol(obj.hFigTraces.hFig, 'Style', 'pushbutton', 'Units', 'normalized', ...
                'Position', [0.17 0.95 0.05 0.035], 'String', 'Apply', ...
                'Callback', @obj.onDurationApply);

            uicontrol(obj.hFigTraces.hFig, 'Style', 'text', 'Units', 'normalized', ...
                'Position', [0.24 0.95 0.11 0.035], 'String', 'Show hidden shank:', ...
                'BackgroundColor', 'k', 'ForegroundColor', 'w', 'HorizontalAlignment', 'left');
            obj.hHiddenPopup = uicontrol(obj.hFigTraces.hFig, 'Style', 'popupmenu', 'Units', 'normalized', ...
                'Position', [0.36 0.95 0.10 0.035], 'String', {'(none hidden)'});
            obj.hShowBtn = uicontrol(obj.hFigTraces.hFig, 'Style', 'pushbutton', 'Units', 'normalized', ...
                'Position', [0.47 0.95 0.05 0.035], 'String', 'Show', ...
                'Callback', @obj.onShowHiddenShank);

            obj.rebuildShankPanels();

            obj.tracesFilt = jrclust.views.plotFigTraces(obj.hFigTraces, obj.hCfg, obj.tracesRaw, 1, obj.hClust);
        end
    end
end
