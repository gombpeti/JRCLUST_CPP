function resFile = jrc_res_path(prmFile, outputDir)
%JRC_RES_PATH Locate the _res.mat belonging to a .prm.
%
%   Follows JRCLUST's own rule - outputDir if the .prm sets one, otherwise the
%   directory holding the .prm - with an override for when that path is stale.

    if nargin < 2, outputDir = ''; end

    if ~exist(prmFile, 'file')
        error('jrc_res_path:noPrm', '%s does not exist', prmFile);
    end

    hCfg = jrclust.Config(prmFile);

    [dirName, sessionName] = fileparts(hCfg.configFile);
    if ~isempty(hCfg.outputDir), dirName = hCfg.outputDir; end
    if ~isempty(outputDir), dirName = outputDir; end

    resFile = fullfile(dirName, [sessionName '_res.mat']);

    if ~exist(resFile, 'file')
        error('jrc_res_path:noRes', ...
              ['%s not found.\n' ...
               'Pass ''outputDir'' if the .prm''s outputDir points elsewhere.'], ...
              resFile);
    end
end
