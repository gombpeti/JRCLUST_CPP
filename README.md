# JRCLUST — native C++/CUDA port

A standalone, compiled spike sorter. No MATLAB, no Python, no interpreter — a
single `.exe` plus the CUDA runtime.

This is a staged port of [JRCLUST](https://github.com/JaneliaSciComp/JRCLUST).
The CUDA kernels are carried over **unmodified**; only the host code that used to
be MATLAB is rewritten.

## Status

**Stage 1 — complete.** Foundation and GPU path, validated.

| Component | State |
|---|---|
| CMake + CUDA + MSVC build | working, builds a standalone `jrclust.exe` |
| `.prm` parameter parsing | working (117 params, incl. 512×2 matrices) |
| Memory-mapped recording reader | working, verified byte-exact on 52 GiB files |
| `jrc_cuda_rho` kernel | working, **bit-exact** vs CPU reference |
| `jrc_cuda_delta` kernel | working, agrees to float rounding (2.4e-7) |

**Stage 2 — detection complete.** Filter → CAR → threshold → peaks → merge.

| Component | State |
|---|---|
| Probe geometry (siteMap/siteLoc/shankMap) | working, 8 shanks × 512 sites |
| `ndiff` filter (orders 1–3) | working |
| Common average reference (mean/median) | working |
| Threshold estimation | working, reproducible (see note) |
| Peak detection + cross-site merge | working |
| Multithreaded, prefetched pipeline | **88× realtime** on NVMe |

**Stage 2b — features and clustering complete.** The full sort runs end to end.

| Component | State |
|---|---|
| Site neighbour table (`findSiteNeighbors`) | working |
| Spike window extraction | working, 21 samples × 9 sites |
| Local CAR + mean subtraction | working |
| PCA (Jacobi eigensolver) + projection | working, PC1 = 66.8% variance |
| Secondary peak sites (`nPeaksFeatures = 2`) | working |
| rho / delta via the original CUDA kernels | working |
| Quadratic detrend + centre selection | working |
| Cluster assignment, small-cluster removal | working |

**Stage 3 — post-processing in C++.** Mean waveforms, waveform similarity,
auto-merge and quality scores all run natively; the MATLAB step is now optional.

| Component | State |
|---|---|
| Mean waveforms (drift and non-drift paths) | working, matches MATLAB |
| Waveform similarity (shifted Pearson) | working, matches MATLAB |
| Auto-merge (`mergeBySim`, multi-pass) | working, matches MATLAB exactly |
| `orderClusters('clusterSites')` | working |
| Quality scores (SNR, Vpp, ISI, IsoDist, L-ratio) | working |
| Feature-weighted spike positions (`spikePos`) | working |
| Refractory-period cleanup | working |

Not yet ported: `bandpass`/`fir1`/`sgdiff`/`user` filters, `fftClean`, motion
rejection (`blankThresh`), the `interpPC` refinement, Mahalanobis outlier
removal, and the Qt curation GUI.

## Validated against MATLAB

Post-processing was checked against JRCLUST in `D:\JRCLUST_LINUX`, run on
byte-identical input (the C++ sorter's own raw clusters, imported into MATLAB
and post-processed there). 60 s of a 512-channel recording:

| | C++ | MATLAB |
|---|---|---|
| raw clusters | 124 | 124 |
| after auto-merge | **114** (10 merged) | **114** (10 merged) |
| median SNR | 6.35 | 6.33 |
| median Vpp | 108.82 µV | 108.54 µV |
| median ISI ratio | 0.09 | 0.09 |
| median isolation distance | **9.86** | **9.86** |
| units SNR≥5 and ISI≤0.1 | **57 / 114** | **57 / 114** |

Getting there surfaced two things worth recording:

- **`driftMerge` changes the answer.** With it on, `getClusterMean` averages
  only the middle third of a cluster by spike y-position, realigning each spike
  onto the cluster's own site neighbourhood. With it off, it averages the
  spikes nearest the middle of the recording. On this data that is 114 vs 120
  clusters. Both paths are implemented; the `.prm` selects.
- **Spike positions must be the real centroids.** Feeding MATLAB a
  site-location stand-in instead of `spikePos`'s feature-weighted centroid
  collapsed the drift bands and changed the merge result (100 vs 114). The
  sorter now computes and exports the real positions.

## End-to-end result

47-minute, 512-channel, 20 kHz recording (358 active sites), USB-attached drive:

```
detected              5,280,729 spikes in 187.0 s
clusters              1,100
assigned spikes       5,217,838 (98.8 %)
cluster size          42 min, 394 median, 207,214 max

detect                187.01 s
rho (kernel)           13.95 s
delta (kernel)         11.57 s
assign                  1.47 s
TOTAL                 240.54 s      (11.5x realtime)
```

**Caveat on cluster quality.** These are raw density-peak clusters. The
post-processing that JRCLUST applies afterwards — auto-merge by waveform
similarity, refractory-period and outlier spike removal — is not ported yet, and
it shows: the largest cluster fires at 74.8 Hz, which is multi-unit rather than
a single neuron. Treat the current cluster count as an upper bound.

## Performance

512-channel, 20 kHz recording on a Ryzen 9 7950X (16 core / 32 thread):

| Version | Throughput |
|---|---|
| First working version (single-threaded) | 2.1× realtime |
| Multithreaded + site-major layout | 12.1× |
| Fused read/convert/transpose | 34.3× |
| Chunk tuning + I/O prefetch | **88.5×** |

A 47-minute recording detects in **32 seconds**. At that point the pipeline is
I/O bound and near-perfectly balanced — 30.1 s of reading overlapped with
29.3 s of compute in 32.0 s of wall clock, at 1836 MB/s off an NVMe drive.

The same recording on a USB-attached drive runs at 15×, limited entirely by
that drive's ~420 MB/s. **Storage, not CPU, is the ceiling** — put the data on
NVMe.

## Deliberate differences from MATLAB

Three places where this port does not reproduce MATLAB's behaviour, each on
purpose:

1. **Threshold estimation is reproducible.** MATLAB's `estimateRMS` subsamples
   with `randsample()`, so its thresholds — and therefore its spike counts —
   differ between runs on identical input. This port takes an evenly spaced
   subsample from fixed time windows instead. Same estimator for stationary
   noise, but the same answer every time.
2. **Chunk size cannot change results.** Threshold windows are defined in
   absolute time, independent of `chunkSeconds`, so chunking is purely a
   performance knob. Verified: chunk sizes of 1/2/4/8/16 s all yield exactly
   225,002 spikes on the same input.
3. **One MATLAB bug is not reproduced.** In `detectPeaksSite.m`, the guard that
   should discard only a peak sitting on the final sample instead falls into an
   `else` branch that removes the last candidate unconditionally, discarding one
   real spike per site per chunk. This port drops it only when it is genuinely
   at the boundary.

4. **One PCA basis for the whole recording.** MATLAB recomputes principal
   vectors per load, so features from different loads sit in different bases and
   are not strictly comparable. This port builds a single basis from spikes
   sampled across the whole recording. On a file that fits in one MATLAB load
   the two are equivalent; on a long file this one is correct and MATLAB's is
   not.
5. **Principal vectors are reproducible too.** `spikePrVecs.m` also subsamples
   with `randsample()`; this port samples evenly.

Arithmetic is done in `float`; MATLAB detects in `int16` with saturation. The
float path is more accurate, but exact spike-for-spike agreement should not be
expected.

6. **The L-ratio is computed correctly here, and differently from MATLAB.**
   `computeQualityScores.m` passes `size(siteFeatures, 2)` as the chi-squared
   degrees of freedom — that is the *spike count*, where the *feature
   dimensionality* belongs. With dof in the thousands the tail probability
   saturates at 1 and the statistic degenerates into `nOutside / nIn`, which is
   not an L-ratio. On the test data MATLAB reports a median of 8.18 where the
   correct value is 0.50. This port uses the feature count; pass
   `--matlab-lratio` to reproduce MATLAB's numbers for comparison.

## Beyond MATLAB: the refractory merge gate

Merging on waveform similarity alone has a known failure mode — two *different*
neurons recorded on neighbouring sites can have near-identical waveform shapes,
and get merged into one contaminated unit. Waveform shape cannot distinguish
these cases, but the spike trains can: one neuron cannot fire twice inside its
own refractory period, so two units that are really one neuron show a dip at
zero lag in their cross-correlogram, while two distinct neurons fire
independently.

`--refrac-merge` adds that test as a veto on top of the waveform decision. For
each proposed merge it counts spike pairs within ±`refracInt` and compares them
against the number expected from independent trains. A ratio near 1 means
independent — the merge is rejected.

The test only fires when it has power: if the expected coincidence count is
below `minExpectedCoincidences` (default 5), the ratio is Poisson noise rather
than evidence and the waveform decision stands. Without that guard, small
units would be vetoed on a coin flip.

Measured on 300 s of the test recording:

| | default (MATLAB parity) | `--refrac-merge` |
|---|---|---|
| clusters after merge | 155 | 162 |
| merges vetoed | — | 22 |
| units SNR≥5 and ISI≤0.1 | 59 | **61** |

More clusters but *more clean units*: the vetoed merges were combining distinct
neurons into contaminated ones. This is off by default so the pipeline
reproduces MATLAB unless you ask otherwise.

GPU and CPU density-peak paths agree closely but not bit-exactly: on a 60 s
slice they gave 277 vs 278 clusters and 89.5% vs 89.9% spikes assigned. The rho
kernel is bit-exact against the CPU reference; delta differs by ~2.4e-7 from
float accumulation order, which is enough to move a spike across the hard
centre-selection threshold. This is inherent to thresholding a float quantity,
not a defect in either path.

## Requirements

Matches the machine this was developed on:

- **CUDA Toolkit 12.8**
- **MSVC 14.29** (v142 toolset, ships with VS2019) — `nvcc` on Windows requires
  MSVC as the host compiler; MinGW is not supported
- **CMake ≥ 3.21** and **Ninja** (both bundled with Qt under `C:\Qt\Tools`)
- Qt 6.7.3 msvc2019_64 — *not needed until Stage 2*

## Build

```
build.bat
```

This pins the v142 toolset, configures with Ninja, and produces
`build\bin\jrclust.exe`.

To build for a different GPU, override the architecture (default is `89`, the
RTX 4060 Ti):

```
cmake -S . -B build -DCMAKE_CUDA_ARCHITECTURES=86
```

Use a semicolon-separated list (`"75;86;89"`) for a portable fat binary.

## Usage

```
jrclust gpuinfo               show CUDA devices
jrclust info <file.prm>       parse config, open recordings, report geometry
jrclust probe <file.prm>      show probe geometry and site neighbours
jrclust detect <file.prm>     filter, threshold and detect spikes
jrclust sort <file.prm>       detect, then cluster (rho/delta on GPU)
jrclust selftest [nSpikes]    verify the CUDA kernels against a CPU reference
```

`sort` options:

```
--seconds <n>         process only the first n seconds
--recording <n>       which entry of rawRecordings (default 0)
--recording-path <p>  override the recording path in the .prm
--output-dir <d>      override outputDir from the .prm
--save                write .jrc files for MATLAB curation
--no-merge            skip auto-merge
--rm-refrac           drop refractory-period violations
--refrac-merge        veto merges whose spike trains fire independently
--matlab-lratio       reproduce MATLAB's (incorrect) L-ratio
--cpu                 run rho/delta on the CPU reference
```

`--recording-path` and `--output-dir` exist because `.prm` files often carry
paths from the machine they were written on. They override without editing the
file.

## Concatenated (multi-session) sorting

When `rawRecordings` lists several files, JRCLUST treats them as **one
continuous recording** and sorts them together, so a unit can be followed across
sessions. That is the default here too — `sort` concatenates every entry unless
`--recording n` selects one.

How it works, following `detect.m`:

1. **Thresholds are per file.** Each session has its own noise level, so a
   shared threshold would bias detection toward whichever file is quieter.
2. **Detection runs per file**, so a file boundary never falls inside a filter
   window or a spike waveform.
3. **Spike times are shifted onto one timeline** by the cumulative sample count,
   exactly as `detect.m` does with `recOffset`.
4. **Clustering, merging and quality scores run once**, over all spikes.

One thing this port does that MATLAB does not: **the principal-vector basis is
pooled across every file**. MATLAB recomputes principal vectors per load, so in
a concatenated sort the features from different files sit in *different bases*
and are not strictly comparable — yet they are clustered together anyway. Here
windows are sampled from all files first, one basis is computed, and every file
projects onto it.

The run reports how units distribute across sessions, e.g.:

```
cross-session:
  units in all 3 sessions: 237 / 369
  units present in exactly 1 session: 37
  units present in exactly 2 sessions: 95
  units present in exactly 3 sessions: 237

  largest 5 units, spikes per session:
    unit  345:     1916     4289     1983
    unit  301:     2429     2959     1724
```

Units appearing in only one session are either genuinely session-specific or
lost to drift, and are worth checking during curation.

## Manual curation in MATLAB

`jrclust.exe sort <prm> --save` is the equivalent of MATLAB's
`jrc detect-sort`. Two commands get you to the curation GUI:

```bat
jrclust.exe sort "H:\data\mysession.prm" --save
```

then in MATLAB:

```matlab
addpath('D:\JRCLUST_CPP\matlab');
jrc_import_cpp('H:\data\mysession.prm')
jrc manual H:\data\mysession.prm
```

### Post-processing and quality scores

`jrc_import_cpp` runs **auto-merge** by default, the same step
`SortController.m` runs after clustering: merge units whose mean waveforms
correlate above `maxUnitSim`, repeated up to `nPassesMerge` times. Optional
extras:

```matlab
jrc_import_cpp(prm)                                  % auto-merge (default)
jrc_import_cpp(prm, 'autoMerge', false)              % raw clusters
jrc_import_cpp(prm, 'rmRefrac', true)                % drop refractory violations
jrc_import_cpp(prm, 'rmOutliers', true)              % Mahalanobis outliers
```

**Quality scores come for free.** Constructing the clustering runs `commit()`,
which computes `unitSNR`, `unitVpp`, `unitISIRatio`, `unitIsoDist`,
`unitLRatio`, `waveformSim` and the mean waveforms — all by JRCLUST's own code,
so they are consistent with the toolbox and shown in the GUI as usual.

These steps are deliberately left in MATLAB rather than reimplemented in C++:
they cost seconds, not minutes, and reusing JRCLUST's tested implementations
means the merged result matches what `jrc detect-sort` would produce from the
same spikes. The expensive half — detection and clustering — is what runs
natively.

`--save` writes JRCLUST's own binaries, byte-compatible with the MATLAB version:

| File | Contents |
|---|---|
| `<session>_filt.jrc` | `int16 [nFiltSamples, nSiteNeighbors, nSpikes]` |
| `<session>_raw.jrc` | `int16 [nRawSamples, nSiteNeighbors, nSpikes]` |
| `<session>_features.jrc` | `single [nFeatures, nPeaks, nSpikes]` |
| `<session>_cpp.bin` | spike table, clustering, thresholds |

`jrc_import_cpp` then builds `<session>_res.mat`.

**Why the MATLAB step exists.** `_res.mat` is a v7.3 (HDF5) file holding ~56
variables, most of which are *derived* — mean waveforms, waveform similarity,
SNR, ISI ratio, isolation distance. Rather than reimplement all of that (and
link HDF5 into the C++ build), the importer hands the spike table to JRCLUST's
own `DensityPeakClustering` constructor and lets `commit()` compute them. Those
fields are therefore produced by JRCLUST itself and are consistent by
construction.

Verified end to end: after import, `JRC.loadFiles()` — the exact call `jrc
manual` makes — returns a populated `DensityPeakClustering` with mean
waveforms, `waveformSim`, and per-unit quality scores.

Note that `--save` writes next to the `.prm` (or to `outputDir` if the `.prm`
sets one), using the `.prm`'s base name as the session name — the same rule
MATLAB JRCLUST follows.

`detect` options:

```
--seconds <n>     process only the first n seconds
--chunk <n>       chunk size in seconds (default 8; performance only)
--recording <n>   which entry of rawRecordings to process (default 0)
```

`selftest` is the correctness gate: it builds a synthetic site, runs both
kernels on the GPU and the equivalent plain-C++ reference, and compares. Run it
after any change to the CUDA path.

Set `JRC_THREADS` to override the worker count (defaults to all hardware
threads).

## Layout

```
src/
  core/
    Config.{h,cpp}            .prm parser  (port of +jrclust/@Config)
    Probe.{h,cpp}             site geometry, neighbour search
    Recording.{h,cpp}         raw reader   (port of +jrclust/+detect/@Recording)
    MemoryMappedFile.{h,cpp}  cross-platform mmap (replaces MATLAB memmapfile)
    Parallel.h                work-stealing parallelFor
    SymEigen.h                Jacobi eigensolver (replaces MATLAB eig)
  detect/
    Filter.{h,cpp}            ndiff + common average reference
    PeakDetector.{h,cpp}      thresholds, per-site peaks, cross-site merge
    DetectController.{h,cpp}  chunked, prefetched pipeline driver
  features/
    FeatureExtractor.{h,cpp}  windows, local CAR, PCA, projection
  sort/
    Clustering.{h,cpp}        rho/delta driver, detrend, cluster assignment
  cuda/
    DensityPeak.{h,cu}        host launches + CPU reference implementations
    CudaUtils.h               error checking, RAII device buffers
    kernels/                  ORIGINAL JRCLUST .cu files — do not edit
  main.cpp                    CLI
```

### Memory layout

Detection buffers are **site-major**: `samples[site * nSamples + t]`, so each
site's trace is contiguous in time. Filtering, threshold estimation and peak
finding all walk one site along time, and this makes those sequential rather
than strided reads — worth roughly 3× on its own.

The file is the opposite (channel-major, sample-interleaved), so
`Recording::readSitesTransposed` fuses read, integer→float conversion and
transpose into one cache-tiled pass. Doing these as separate passes cost 65% of
total runtime; fused, the transpose disappears from the profile entirely.

### On the kernels

`src/cuda/kernels/*.cu` are byte-for-byte copies of `+jrclust/+CUDA/*.cu`. They
are `#include`d by `DensityPeak.cu` rather than compiled as separate translation
units, so they must **not** be added to the CMake source list.

MATLAB loaded these as PTX at runtime via `parallel.gpu.CUDAKernel`. Here they
are compiled directly into the executable by `nvcc`, which removes the JIT step
and the MATLAB marshaling layer. The kernel code itself is untouched, so GPU
results are identical.

## Data layout

Recordings are channel-interleaved: sample `(channel c, sample s)` lives at
`headerOffset + (c + s*nChans) * sizeof(dtype)`. This is what MATLAB's
column-major `memmapfile` with `dshape = [nChans nSamples]` produced, so files
written for the MATLAB version load unchanged.

## Performance note

The GPU sort is **not** expected to get faster — it was already CUDA. The wins
are in file I/O, filtering, detection, feature extraction, and the elimination
of MATLAB startup and marshaling overhead. The guaranteed win is deployment: a
self-contained executable with no runtime or license dependency.
