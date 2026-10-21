# JRCLUST, patched

A copy of JRCLUST (Apache 2.0, see `LICENSE`) with the fixes below applied.
The original is at <https://github.com/JaneliaSciComp/JRCLUST>; `json/jrc.json`
still identifies this as **4.0.1 "Edward"**, which is what the GUI title bar
reports.

Used by the C++/CUDA sorter in this repository for manual curation. Add both
to the MATLAB path, and make sure no other JRCLUST copy is on the path or is
the current folder, since MATLAB prefers the current folder and will silently
load a different one.

## Fixes

### `+curate/@CurateController/deleteClusters.m` — Delete annotated ran out of memory

The loop appended a snapshot of the whole clustering per deleted unit and held
them all at once:

```matlab
res = struct('spikeClusters', {}, 'metadata', {});
for iCluster = 1:numel(deleteMe)
    res_ = obj.hClust.deleteUnit(args{:});
    res(end+1) = res_;        % every snapshot kept
end
```

Each snapshot carries the spike table plus every per-cluster field —
`spikesByCluster`, four mean-waveform arrays, the nClusters × nClusters
similarity matrix. On a 512-site, 1143-unit recording that is **514 MB each**,
so deleting 120 annotated units asked for **60 GB** and MATLAB stalled partway
through. Only `res(end)` was ever read — as input to the next iteration and as
the argument to `commit` — so the other snapshots were pure waste. Keeping one
makes the cost independent of how many units are deleted.

This is why it worked on smaller probes and stopped working on large ones,
rather than failing outright.

### `+interfaces/@Clustering/computeQualityScores.m` — L-ratio used the wrong degrees of freedom

```matlab
unitLRatio_(iCluster) = sum(1 - chi2cdf(mDistsOutside, size(siteFeatures,2)))/nSpikesOnSite;
```

`siteFeatures` is `squeeze(spikeFeatures(:, 1, idx))`, i.e. **nFeatures ×
nSpikes** — which is why `mahal` is handed its transpose and why `isOnSite`
indexes columns. So `size(siteFeatures, 2)` is the **spike count**, while the
chi-squared degrees of freedom must be the dimensionality of the Mahalanobis
distance, the **feature count** (`size(..., 1)`).

Using the spike count made the survival function hugely permissive: every
outside spike counted as "could belong". On a cleanly separated synthetic unit
the correct L-ratio is 0.00 and the buggy one 6.08; on a measured unit, 0.50
against 8.18. Lower L-ratio means better isolation, so this made every unit
look far worse isolated than it is.

### `+utils/subsample.m` — results were not reproducible

`randsample` drew a different subset on every call. `subsample` feeds
`getPVSamples` and `spikePCA`, so the principal vectors were re-estimated from
a fresh random 10,000 spikes each time and the same recording loaded twice gave
two different feature spaces — which for a sort that concatenates sessions is
the difference between features being comparable across sessions and not. It
also feeds `computeRho`, so the clustering itself varied between runs on
identical input.

Replaced with evenly spaced indices. That is not merely deterministic but the
better estimator here: the inputs are ordered by time, so a regular stride
draws from the whole recording in proportion, whereas a random draw is uniform
only in expectation and can over-weight a high-rate stretch in any single run.

### Selecting more than two units

`updateSelect.m` truncated the selection to `iClusters(1:2)`, so the whole GUI
was capped at a pair. It now accepts any number, in click order, so the first
unit stays primary and the pairwise views keep following it.

`mouseClickFigWav.m` previously replaced the second slot on a right click —
`[obj.selected(1) iCluster]` — which made a selection of three impossible.
Ctrl-click, right click and shift-click now ADD a unit, and clicking a
selected unit again removes it.

`updateFigCorr.m` draws an **N x N correlogram matrix** when more than two
units are selected: autocorrelograms down the diagonal (tinted), cross-
correlograms off it. One and two units behave exactly as before.

That matrix is what makes a group of units decidable in one look, because a
cross-correlogram answers a question waveform similarity cannot:

| cross-correlogram | means |
| --- | --- |
| dip at zero | the two never fire together within a millisecond, which only one refractory period explains - **one cell split in two** |
| flat | independent - **different cells**, however alike their waveforms |
| peak at zero | the same event counted twice - duplicate detections |

`mergeSelected.m` had to change with it. It took `min(selected)` and
`max(selected)`, which was correct while the cap was two but would, with five
units selected, have merged the outermost pair and **silently left the middle
three behind** while reporting success. It now merges all of them, highest id
first so the renumbering after each step cannot invalidate the ids still
pending.

Other views (FigProj, FigTime, FigHist) still read `selected(1)` and
`selected(2)`, so they show the first two and degrade gracefully rather than
breaking.

### `+curate/@CurateController/addMenu.m` — three menu items

- **File → Save As Pruned** — writes a reduced copy with the discarded spikes
  and clusters removed from the `_res.mat` *and* the `.jrc` binaries. `Edit >
  Delete annotated` only rewrites the spike table, so the `.jrc` files keep
  every discarded spike and nothing downstream shrinks.
- **Edit → Merge auto (cross-correlogram checked)** — proposes pairs by
  waveform correlation as the built-in does, then accepts them only if the
  cross-correlogram between the two spike trains shows a shared refractory
  dip. Waveform shape alone cannot separate one neuron split in two from two
  neurons at the same distance from the same sites; on an 841-unit recording a
  correlation cut of 0.95 proposed 200 pairs and only 91 survived the check.
- **Edit → Relabel units** — re-judges every unit and rewrites the
  annotations.

All three call scripts in the parent `matlab/` folder and report what to do if
those are not on the path.
