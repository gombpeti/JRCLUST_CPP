===============================================================================
 JRCLUST - native C++/CUDA spike sorter
 with MATLAB helpers for curation
===============================================================================

A compiled spike sorter. No MATLAB and no Python are needed to run it: one
.exe plus the CUDA runtime. MATLAB is only used afterwards, for manual
curation in JRCLUST's existing GUI.

The CUDA kernels are the ORIGINAL JRCLUST kernels, copied byte for byte. Only
the host code that used to be MATLAB has been rewritten.

  executable      D:\JRCLUST_CPP\build\bin\jrclust.exe
  MATLAB helpers  D:\JRCLUST_CPP\matlab\
  reference JRC   D:\JRCLUST_LINUX          (needed by the MATLAB helpers)


===============================================================================
 1. QUICK START
===============================================================================

Sort a session and prepare it for curation:

    D:\JRCLUST_CPP\build\bin\jrclust.exe sort "H:\path\mysession.prm" --save

Then in MATLAB:

    addpath(genpath('D:\JRCLUST_LINUX'));
    addpath('D:\JRCLUST_CPP\matlab');

    jrc_import_cpp('H:\path\mysession.prm', 'autoMerge', false)
    jrc manual H:\path\mysession.prm

Use 'autoMerge', false because the C++ side has already merged. Leaving it on
would merge a second time.


===============================================================================
 2. BUILDING
===============================================================================

    D:\JRCLUST_CPP\build.bat

Requirements (all already installed on this machine):

    CUDA Toolkit 12.8
    MSVC 14.29  (the v142 toolset from VS2019)
    CMake + Ninja  (bundled with Qt under C:\Qt\Tools)

MSVC is not optional. nvcc on Windows only supports MSVC as host compiler, so
MinGW cannot be used.

For a different GPU, override the architecture (default "89;120" = RTX 4060 Ti
+ RTX 5070 Ti; sm_120 needs CUDA 12.8 or newer):

    cmake -S . -B build -DCMAKE_CUDA_ARCHITECTURES=86

Check the build is sane at any time:

    jrclust.exe selftest

This runs both CUDA kernels against a plain C++ reference and compares. Run it
after any change to the CUDA path.


===============================================================================
 3. COMMANDS
===============================================================================

    jrclust gpuinfo                 list CUDA devices
    jrclust info    <file.prm>      parse config, open recordings, report
    jrclust probe   <file.prm>      probe geometry and site neighbours
    jrclust detect  <file.prm>      filter, threshold, detect spikes only
    jrclust sort    <file.prm>      detect + cluster + merge + classify
    jrclust selftest [nSpikes]      verify the CUDA kernels

'sort' is the equivalent of MATLAB's  jrc detect-sort.

-------------------------------------------------------------------------------
 3.1  sort options
-------------------------------------------------------------------------------

  --seconds <n>          process only the first n seconds (per file)
  --recording <n>        sort ONLY entry n of rawRecordings.
                         Default is to concatenate ALL of them (see section 4).
  --recording-path <p>   override the recording path from the .prm
  --output-dir <d>       override outputDir from the .prm
  --save                 write the .jrc files needed for MATLAB curation
  --cpu                  run rho/delta on the CPU instead of the GPU

  --no-merge             skip auto-merge
  --refrac-merge         veto merges whose spike trains fire independently
  --rm-refrac            drop refractory-period violations within a unit
  --rm-outliers          drop spikes far from their own unit in feature space
  --no-split             do not split clusters holding two populations

  --no-classify          skip noise/single/multi labelling
  --mark-for-delete      label noise units 'to_delete' instead of 'noise',
                         which arms JRCLUST's Edit > Delete annotated
  --label-noise-only     label only suspect units, leaving good ones blank
  --max-vpp <uV>         treat units above this amplitude as suspect
  --matlab-lratio        reproduce MATLAB's L-ratio (see section 8)

-------------------------------------------------------------------------------
 3.2  detect options
-------------------------------------------------------------------------------

  --seconds <n>          process only the first n seconds
  --chunk <n>            chunk size in seconds (default 8, performance only)
  --recording <n>        which entry of rawRecordings (default 0)


===============================================================================
 4. CONCATENATED (MULTI-SESSION) SORTING
===============================================================================

When rawRecordings lists several files, they are treated as ONE continuous
recording and sorted together, so a unit can be followed across sessions.
This is the default, and matches JRCLUST.

    jrclust.exe sort "mysession.prm" --save            all files together
    jrclust.exe sort "mysession.prm" --recording 0     only the first file

How it works:

  1. Thresholds are computed PER FILE. Each session has its own noise level,
     so a shared threshold would bias detection toward the quieter file.
  2. Detection runs per file, so a file boundary never falls inside a filter
     window or a spike waveform.
  3. Spike times are shifted onto one timeline by the cumulative sample count.
  4. Clustering, merging and quality scores run once, over all spikes.

The run reports how units distribute across sessions:

    cross-session:
      units in all 3 sessions: 237 / 369
      units present in exactly 1 session: 37

Units appearing in only one session are either genuinely session-specific or
lost to drift, and are worth checking during curation.

NOTE: all files write to the same outputDir under the same session name. If
you sort them separately with --recording, give each its own --output-dir or
the second run overwrites the first.


===============================================================================
 5. UNIT CLASSIFICATION (noise / single / multi)
===============================================================================

After merging, every unit is labelled:

    0  noise    the waveform is not an action potential
    1  single   clean, well isolated
    2  multi    a real spike shape, but contaminated

These become JRCLUST annotations, visible in the curation GUI.

The decisive noise tests are on SHAPE, not amplitude, because shape is
amplitude-invariant: a large unit right against the electrode still has the
stereotyped shape of an action potential, while a digital glitch does not,
however large or small it is. A plain microvolt ceiling would eventually
discard a real close-by unit while missing every low-amplitude artifact.

Measured on each unit's mean RAW waveform:

    trough half-width    a real spike is 0.15-0.5 ms wide; a glitch is 1 sample
    spikiness            trough vs its immediate neighbours; >2.5 = glitch
    recovery             a real spike repolarises after the trough
    spatial spread       real units span a few neighbouring sites
    SNR                  real units are 3-20; >50 is not biological
    ISI ratio            refractory violations

Bimodal split
-------------
Before classifying, clusters holding TWO separable populations are split in
two and each half judged on its own. This automates what curation does by
hand: divide a suspect cluster, see that one half is clean spikes and the
other is artifact, keep the first.

This matters because whole-cluster tests cannot catch that case. A cluster
that is 98% neuron and 2% artifact still has a clean mean waveform, so nothing
flags it - yet those 2% are the wild outliers you see in the amplitude view.

Turn it off with --no-split.


===============================================================================
 6. MATLAB HELPERS
===============================================================================

Always add both paths first:

    addpath(genpath('D:\JRCLUST_LINUX'));
    addpath('D:\JRCLUST_CPP\matlab');

-------------------------------------------------------------------------------
 6.1  jrc_import_cpp  - build the _res.mat for the curation GUI
-------------------------------------------------------------------------------

    jrc_import_cpp('mysession.prm', 'autoMerge', false)

Reads what  sort --save  wrote and produces <session>_res.mat, so that
jrc manual opens it exactly as if produced by jrc detect-sort.

Options:
    'autoMerge'   default true. Set FALSE when the C++ side already merged.
                  NOTE: merging here renumbers clusters, so the sorter's unit
                  labels no longer line up and are dropped.
    'rmOutliers'  Mahalanobis outlier removal (default false)
    'rmRefrac'    drop refractory violations (default false)
    'outputDir'   where the .jrc files are, if the .prm's outputDir is stale

Mean waveforms, waveform similarity, SNR, ISI ratio, isolation distance and
L-ratio are NOT taken from the C++ side here - they are computed by JRCLUST's
own DensityPeakClustering constructor, so they are consistent with the rest of
the toolbox.

-------------------------------------------------------------------------------
 6.2  jrc_mark_noise  - set annotations without the GUI
-------------------------------------------------------------------------------

    jrc_mark_noise('mysession.prm')

Renames annotations: units labelled 'noise' become 'to_delete'. It DELETES
NOTHING and does not touch single, multi or blank units.

Why it exists: JRCLUST's Edit > Delete annotated only recognises the exact
string 'to_delete'. A unit labelled 'noise' is invisible to it. This script is
just a translator between the two vocabularies.

    jrc_mark_noise(prm, 'units', [3 17 42])      mark specific unit IDs
    jrc_mark_noise(prm, 'undo', true)            'to_delete' back to 'noise'
    jrc_mark_noise(prm, 'dryRun', true)          report, write nothing

If you ran the sorter with --mark-for-delete, the labels are already
'to_delete' and this script has nothing to do. Go straight to jrc_prune.

CAUTION:
    jrc_mark_noise(prm, 'from', {'noise','multi'})
marks every multi unit too. Multi units are real spikes that are merely
contaminated, not junk. Do not include 'multi' unless you mean it.

-------------------------------------------------------------------------------
 6.3  jrc_prune  - delete marked units everywhere, write a pruned copy
-------------------------------------------------------------------------------

    jrc_prune('mysession.prm', 'dryRun', true)     preview
    jrc_prune('mysession.prm')                     do it

Drops every spike in a unit annotated 'to_delete', plus anything already at
spikeClusters == 0, renumbers the survivors, and writes a complete session to
<prm dir>\pruned\ :

    <name>_pruned.prm
    <name>_pruned_res.mat
    <name>_pruned_raw.jrc / _filt.jrc / _features.jrc / _hist.jrc

Open it with:

    jrc manual <prm dir>\pruned\<name>_pruned.prm

Options:
    'outputDir'    where the input _res.mat and .jrc files live
    'prunedDir'    where to write (default <prm dir>\pruned)
    'chunkSpikes'  spikes per I/O block (default 50000). Lower if memory is
                   tight, raise for speed.
    'dryRun'       report what would happen, write nothing

WHY THIS EXISTS
The GUI's Edit > Delete annotated runs out of memory on large sessions. It
subsets the waveform arrays with vals(:,:,indices), and spikesRaw for 15
million spikes is about 11 GB, so the subset needs another 11 GB contiguous
before anything is freed. This script never holds those arrays: it streams the
.jrc files a block of spikes at a time, so peak memory is set by chunkSpikes
rather than by the size of the recording.

Run the dry run first. Pruning writes a new folder and is safe, but the
marking step edits your _res.mat in place.


===============================================================================
 7. TYPICAL WORKFLOWS
===============================================================================

A. Sort, curate by hand
-----------------------
    jrclust.exe sort "mysession.prm" --save
    (MATLAB)
    jrc_import_cpp('mysession.prm', 'autoMerge', false)
    jrc manual mysession.prm

B. Sort, auto-label noise, review, bulk delete
----------------------------------------------
    jrclust.exe sort "mysession.prm" --save --mark-for-delete
    (MATLAB)
    jrc_import_cpp('mysession.prm', 'autoMerge', false)
    jrc manual mysession.prm
      -> View > Display annotated units only     to review the suspects
      -> Edit > Delete annotated                 to remove them

C. Same, but the GUI runs out of memory on delete
-------------------------------------------------
    jrclust.exe sort "mysession.prm" --save --mark-for-delete
    (MATLAB)
    jrc_import_cpp('mysession.prm', 'autoMerge', false)
    jrc_prune('mysession.prm', 'dryRun', true)
    jrc_prune('mysession.prm')
    jrc manual <prm dir>\pruned\<name>_pruned.prm

D. Quick test before committing to a long run
---------------------------------------------
    jrclust.exe sort "mysession.prm" --seconds 120
    (no --save, so it reports and writes nothing)


===============================================================================
 8. DIFFERENCES FROM MATLAB JRCLUST
===============================================================================

Deliberate, each for a reason:

1. Threshold estimation is reproducible. MATLAB's estimateRMS subsamples with
   randsample(), so its thresholds - and therefore its spike counts - differ
   between runs on identical input. This port samples evenly instead.

2. Chunk size cannot change results. Threshold windows are defined in absolute
   time, independent of --chunk, so chunking is purely a performance knob.

3. A MATLAB bug is not reproduced. In detectPeaksSite.m the guard that should
   discard only a peak on the final sample instead removes the last candidate
   unconditionally, discarding one real spike per site per chunk.

4. One PCA basis for the whole recording. MATLAB recomputes principal vectors
   per load, so in a concatenated sort the features from different files sit
   in DIFFERENT bases yet are clustered together anyway. Here windows are
   sampled from all files first and every file projects onto one basis.

5. The L-ratio is computed correctly. computeQualityScores.m passes
   size(siteFeatures,2) as the chi-squared degrees of freedom - that is the
   SPIKE COUNT, where the feature dimensionality belongs. With dof in the
   thousands the statistic degenerates into nOutside/nIn. On test data MATLAB
   reports 8.18 where the correct value is 0.50. Use --matlab-lratio to
   reproduce MATLAB's numbers for comparison.

Arithmetic is done in float; MATLAB detects in int16 with saturation. The
float path is more accurate, but exact spike-for-spike agreement should not be
expected.

Validated: on identical input, the C++ post-processing reproduces MATLAB
exactly - 114 clusters from 124 in both, median SNR 6.35 vs 6.33, isolation
distance 9.86 vs 9.86, 57/114 clean units in both.


===============================================================================
 9. PERFORMANCE
===============================================================================

The pipeline is usually limited by DISK, not CPU. Measured on a 512-channel
20 kHz recording:

    USB-attached drive    248 MB/s     ~11x realtime
    NVMe (Samsung 970)   1836 MB/s     ~88x realtime

If CPU sits at 12-30%, it is waiting for the disk, not short of cores - all
threads are already in use. Moving the data to NVMe is worth about 7x and
beats anything more parallelism could give.

A 46-minute 512-channel recording sorts in about 30 seconds from NVMe.

Chunk size (--chunk) affects speed only, never results. 8 seconds is the
measured sweet spot.

Set JRC_THREADS to override the worker count (defaults to all hardware
threads).


===============================================================================
 10. TROUBLESHOOTING
===============================================================================

"Invalid raw file location in param file" (MATLAB)
    The .prm points at a drive that no longer exists. MATLAB's jrclust.Config
    refuses to load it at all. Fix the paths in the .prm, or work on a copy.
    The C++ side is more forgiving: use --recording-path and --output-dir.

"Out of memory" on Edit > Delete annotated
    Use jrc_prune instead. See section 6.3.

"Dimensions of arrays being concatenated are not consistent" when clicking a
cluster
    spikesBySite must be a 1 x nSites ROW cell array. jrc_import_cpp and
    jrc_prune both build it correctly; this error means a _res.mat produced
    some other way.

"All units selected" from View > Display annotated units only
    That message means EVERY unit is annotated, which is the default (every
    unit gets noise/single/multi). It is not an error. To make that filter
    show only the suspects, sort with --label-noise-only.

"Subset is empty; selecting all units"
    The opposite: no unit is annotated. Either classification was skipped, or
    nothing was flagged.

Output files land somewhere unexpected
    The sorter follows JRCLUST's rule: outputDir from the .prm if set,
    otherwise the .prm's own directory, with the .prm's base name as the
    session name. Override with --output-dir.

Sorting seems slow
    Check the stage breakdown printed by 'detect'. If it says
    "limited by I/O", see section 9.


===============================================================================
 11. OUTPUT FILES
===============================================================================

Written by  sort --save  :

    <session>_filt.jrc       int16   [nFiltSamples, nSiteNeighbors, nSpikes]
    <session>_raw.jrc        int16   [nRawSamples,  nSiteNeighbors, nSpikes]
    <session>_features.jrc   single  [nFeatures, nPeaks, nSpikes]
    <session>_cpp.bin        spike table, clustering, labels

Written by  jrc_import_cpp  :

    <session>_res.mat        the file jrc manual opens

All .jrc files are MATLAB column-major, byte-compatible with MATLAB JRCLUST.

Rough sizes for 5.3 million spikes: filt 1.9 GB, raw 3.6 GB, features 0.8 GB.


===============================================================================
 end
===============================================================================
