function success = compileCUDA(nvccPath)
    %COMPILECUDA Compile CUDA codes for JRCLUST.
    gpuD = gpuDevice(1);
    if nargin < 1 || isempty(nvccPath)
        nvccPath = localDefaultNvccPath(gpuD);
    end

    basedir = fullfile(jrclust.utils.basedir(), '+jrclust', '+CUDA');
    cudaFiles = dir(fullfile(basedir, '*.cu'));
    cudaFiles = {cudaFiles.name};

    t1 = tic;
    disp('Compiling CUDA codes...');
    success = true;

    [archFlag, computeCapability] = localCudaArchFlag(gpuD);
    fprintf('Detected GPU: %s\n', localGpuName(gpuD));
    fprintf('Using CUDA architecture: %s (compute capability %s)\n', ...
        archFlag, computeCapability);
    nvccVersion = localNvccVersion(nvccPath, gpuD);
    if ~localCanCompile(computeCapability, nvccVersion)
        success = false;
        return;
    end

    nvccPath = localQuotePath(nvccPath);
    hostCompilerFlag = localHostCompilerFlag();
    if ~isempty(hostCompilerFlag)
        fprintf('Using Visual Studio compiler: %s\n', hostCompilerFlag);
    end

    for i = 1:numel(cudaFiles)
        iFileCU = fullfile(basedir, cudaFiles{i});
        iFilePTX = fullfile(basedir, strrep(cudaFiles{i}, '.cu', '.ptx'));
        cmd = sprintf('%s -ptx -m 64 -arch=%s %s "%s" --output-file "%s"', ...
            nvccPath, archFlag, hostCompilerFlag, iFileCU, iFilePTX);

        fprintf('\t%s\n\t', cmd);

        try
            status = system(cmd);
            success = success && (status == 0);
        catch ME
            warning('Could not compile %s: %s\n', iFileCU, ME.message);
            success = false;
        end
    end

    if ~success
        warning(['CUDA could not be compiled but JRCLUST may work fine. ', ...
            'If not, install CUDA toolkit v%0.1f and run "jrc compile".\n'], ...
            gpuD.ToolkitVersion);
    end

    fprintf('Finished compiling, took %0.1fs\n', toc(t1));
end

function nvccPath = localDefaultNvccPath(gpuD)
    if ispc()
        cudaRoot = fullfile('C:\Program Files\NVIDIA GPU Computing Toolkit', 'CUDA');
        toolkitDir = sprintf('v%0.1f', gpuD.ToolkitVersion);
        nvccPath = fullfile(cudaRoot, toolkitDir, 'bin', 'nvcc.exe');

        if ~exist(nvccPath, 'file')
            cudaDirs = dir(fullfile(cudaRoot, 'v*'));
            if ~isempty(cudaDirs)
                [~, iDir] = max([cudaDirs.datenum]);
                nvccPath = fullfile(cudaRoot, cudaDirs(iDir).name, 'bin', 'nvcc.exe');
            end
        end
    else
        [ecode, out] = system('which nvcc');
        if ecode == 0
            nvccPath = strip(out);
        else
            nvccPath = '/usr/local/cuda/bin/nvcc';
        end
    end
end

function [archFlag, computeCapability] = localCudaArchFlag(gpuD)
    gpuName = lower(localGpuName(gpuD));

    if isprop(gpuD, 'ComputeCapability') && ~isempty(gpuD.ComputeCapability)
        computeCapability = char(gpuD.ComputeCapability);
    elseif contains(gpuName, '5090') || contains(gpuName, '5080') || ...
            contains(gpuName, '5070') || contains(gpuName, '5060') || ...
            contains(gpuName, '5050')
        computeCapability = '12.0';
    elseif contains(gpuName, '4090') || contains(gpuName, '4080') || ...
            contains(gpuName, '4070') || contains(gpuName, '4060')
        computeCapability = '8.9';
    elseif contains(gpuName, '3090') || contains(gpuName, '3080') || ...
            contains(gpuName, '3070') || contains(gpuName, '3060')
        computeCapability = '8.6';
    elseif contains(gpuName, '2080') || contains(gpuName, '2070') || ...
            contains(gpuName, '2060')
        computeCapability = '7.5';
    else
        warning('Could not determine GPU compute capability; falling back to sm_86.');
        computeCapability = '8.6';
    end

    archNumber = regexprep(computeCapability, '[^\d]', '');
    archFlag = sprintf('sm_%s', archNumber);
end

function canCompile = localCanCompile(computeCapability, nvccVersion)
    canCompile = true;
    ccValue = str2double(computeCapability);

    if ccValue >= 12
        if localMatlabReleaseNumber() <= 2024.2
            warning(['This GPU has compute capability %s (Blackwell / RTX 50 series). ', ...
                'MATLAB R2024b natively supports NVIDIA GPUs only up to compute capability 9.x. ', ...
                'Use a newer MATLAB release with native 12.x support, or try MATLAB CUDA forward compatibility.'], ...
                computeCapability);
        end

        if ~isnan(nvccVersion) && nvccVersion < 12.8
            warning(['CUDA Toolkit %0.1f cannot compile sm_120 for RTX 50-series GPUs. ', ...
                'Install CUDA Toolkit 12.8 or newer, then run "jrc compile" again.'], ...
                nvccVersion);
            canCompile = false;
        elseif isnan(nvccVersion)
            warning(['Could not determine the CUDA Toolkit version. ', ...
                'RTX 50-series GPUs require CUDA Toolkit 12.8 or newer for sm_120 compilation.']);
        end
    end
end

function nvccVersion = localNvccVersion(nvccPath, gpuD)
    nvccVersion = NaN;
    nvccPath = regexprep(char(nvccPath), '^"|"$', '');

    [ecode, out] = system(sprintf('"%s" --version', nvccPath));
    if ecode == 0
        token = regexp(out, 'release\s+(\d+\.\d+)', 'tokens', 'once');
        if ~isempty(token)
            nvccVersion = str2double(token{1});
            return;
        end
    end

    if isprop(gpuD, 'ToolkitVersion') && ~isempty(gpuD.ToolkitVersion)
        nvccVersion = gpuD.ToolkitVersion;
    end
end

function releaseNumber = localMatlabReleaseNumber()
    releaseName = version('-release');
    token = regexp(releaseName, '^R?(\d{4})([ab])$', 'tokens', 'once');
    if isempty(token)
        releaseNumber = Inf;
        return;
    end

    releaseNumber = str2double(token{1});
    if strcmpi(token{2}, 'b')
        releaseNumber = releaseNumber + 0.2;
    else
        releaseNumber = releaseNumber + 0.1;
    end
end

function hostCompilerFlag = localHostCompilerFlag()
    hostCompilerFlag = '';
    if ~ispc()
        return;
    end

    clDir = localFindClDir();
    if isempty(clDir)
        warning(['Could not find Visual Studio cl.exe automatically. ', ...
            'If nvcc fails, run MATLAB from a Visual Studio Developer Command Prompt.']);
        return;
    end

    hostCompilerFlag = sprintf('-ccbin "%s"', clDir);
end

function clDir = localFindClDir()
    clDir = localFindClWithVswhere();
    if ~isempty(clDir)
        return;
    end

    clDir = localFindClOnPath();
    if ~isempty(clDir)
        return;
    end

    roots = { ...
        fullfile('C:\Program Files (x86)', 'Microsoft Visual Studio', '2019'), ...
        fullfile('C:\Program Files', 'Microsoft Visual Studio', '2019'), ...
        fullfile('C:\Program Files (x86)', 'Microsoft Visual Studio', '2022'), ...
        fullfile('C:\Program Files', 'Microsoft Visual Studio', '2022')};

    for iRoot = 1:numel(roots)
        clFiles = dir(fullfile(roots{iRoot}, '**', 'bin', 'Hostx64', 'x64', 'cl.exe'));
        if ~isempty(clFiles)
            [~, iCl] = max([clFiles.datenum]);
            clDir = clFiles(iCl).folder;
            return;
        end
    end

    clDir = '';
end

function clDir = localFindClOnPath()
    [ecode, out] = system('where cl');
    if ecode ~= 0 || isempty(strtrim(out))
        clDir = '';
        return;
    end

    lines = regexp(strtrim(out), '\r?\n', 'split');
    clDir = '';
    for iLine = 1:numel(lines)
        [candidateDir, ~, ~] = fileparts(lines{iLine});
        if localIsX64ClDir(candidateDir)
            clDir = candidateDir;
            return;
        end
    end
end

function clDir = localFindClWithVswhere()
    vswherePath = fullfile('C:\Program Files (x86)', 'Microsoft Visual Studio', ...
        'Installer', 'vswhere.exe');
    if ~exist(vswherePath, 'file')
        clDir = '';
        return;
    end

    cmd = sprintf(['"%s" -version [16.0,17.0) ', ...
        '-requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 ', ...
        '-property installationPath'], vswherePath);
    [ecode, out] = system(cmd);
    if ecode ~= 0 || isempty(strtrim(out))
        clDir = '';
        return;
    end

    installPaths = regexp(strtrim(out), '\r?\n', 'split');
    for iPath = 1:numel(installPaths)
        clFiles = dir(fullfile(installPaths{iPath}, 'VC', 'Tools', 'MSVC', ...
            '*', 'bin', 'Hostx64', 'x64', 'cl.exe'));
        if ~isempty(clFiles)
            [~, iCl] = max([clFiles.datenum]);
            clDir = clFiles(iCl).folder;
            return;
        end
    end

    clDir = '';
end

function tf = localIsX64ClDir(clDir)
    clDirLower = lower(char(clDir));
    tf = contains(clDirLower, [filesep, 'hostx64', filesep, 'x64']);
end

function name = localGpuName(gpuD)
    if isprop(gpuD, 'Name') && ~isempty(gpuD.Name)
        name = char(gpuD.Name);
    else
        name = 'Unknown GPU';
    end
end

function pathOut = localQuotePath(pathIn)
    pathOut = char(pathIn);
    if isempty(pathOut)
        return;
    end

    if pathOut(1) ~= '"'
        pathOut = sprintf('"%s"', pathOut);
    end
end
