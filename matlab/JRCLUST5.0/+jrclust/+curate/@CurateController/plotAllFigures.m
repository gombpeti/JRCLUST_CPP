function plotAllFigures(obj)
    %PLOTALLFIGURES Plot all figures
    if isempty(obj.hFigs)
        obj.spawnFigures();
    elseif ~all(obj.figApply(@(hFig) hFig.isReady)) % clean up from an aborted session
        obj.closeFigures();
        obj.spawnFigures();
    end

    % plot sim score figure
    if obj.hasFig('FigSim')
        % set key and mouse handles
        hFigSim = obj.hFigs('FigSim');
        obj.updateFigSim();
        hFigSim.hFunKey = @(hO, hE) obj.runExclusive(@obj.keyPressFigSim, hO, hE);
        hFigSim.setMouseable(@(varargin) obj.runExclusive(@obj.mouseClickFigSim, varargin{:}));
    end

    % plot feature projection
    if obj.hasFig('FigProj')
        % set key and mouse handles
        hFigProj = obj.hFigs('FigProj');
        hFigProj.hFunKey = @(hO, hE) obj.runExclusive(@obj.keyPressFigProj, hO, hE);
        hFigProj.setMouseable(); % no special mouse function
    end

    % plot amplitude vs. time
    if obj.hasFig('FigTime')
        % set key and mouse handles
        hFigTime = obj.hFigs('FigTime');
        hFigTime.hFunKey = @(hO, hE) obj.runExclusive(@obj.keyPressFigTime, hO, hE);
    end

    % plot main waveform view
    if obj.hasFig('FigWav')
        % set key and mouse handles
        hFigWav = jrclust.views.plotFigWav(obj.hFigs('FigWav'), obj.hClust, obj.maxAmp, obj.channel_idx, obj.showSubset);
        setFigWavXTicks(hFigWav, obj.hClust, 1); % show cluster counts by default

        hFigWav.hFunKey = @(hO, hE) obj.runExclusive(@obj.keyPressFigWav, hO, hE);
        hFigWav.setMouseable(@(varargin) obj.runExclusive(@obj.mouseClickFigWav, varargin{:}));

        % make this guy the key log
        hFigWav.figApply(@set, 'CloseRequestFcn', @obj.killFigWav);
        obj.addMenu(hFigWav);
    end

    % plot rho-delta figure
    obj.updateFigRD();
    obj.updateFigGuide();

    % update help texts
    helpFigs = fieldnames(obj.helpTexts);
    for i = 1:numel(helpFigs)
        figName = helpFigs{i};

        if obj.hasFig(figName)
            hFig = obj.hFigs(figName);
            hFig.figData.helpText = strjoin(obj.helpTexts.(figName), '\n');
        end
    end

    % select first cluster (also plots other figures)
    obj.updateSelect(obj.showSubset(1), 1);
end
