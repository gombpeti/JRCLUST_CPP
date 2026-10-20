# What to run, in what order

The sorter produces clusters; the scripts below turn those into units you can
trust. The order matters, because each step assumes the previous one has run.

Throughout, `prm` is your parameter file and everything is driven from it:

```matlab
addpath(genpath('D:\JRCLUST_CPP\matlab\JRCLUST5.0'));
addpath('D:\JRCLUST_CPP\matlab');
prm = 'F:\...\amplifier_denoised.prm';
```

Make sure no other JRCLUST copy is on the path or is your current folder.
MATLAB prefers the current folder, so sitting inside a different JRCLUST
silently loads that one instead - which shows up as a `.prm` failing to
validate rather than as anything obviously version related.

---

## The short version

```
1.  jrclust.exe sort ...          detect, cluster, merge          (cmd, ~12 min)
2.  jrc_import_cpp(prm, ...)      build the _res.mat               (MATLAB)
3.  jrc_relabel(prm)              annotate suspect units
4.  jrc manual prm                LOOK AT IT
5.  jrc_prune(prm)                remove what you agreed to
```

Steps 2-5 are repeatable without re-sorting. Only step 1 is expensive, and
only step 5 removes anything.

---

## 1. Sort

```
D:\JRCLUST_CPP\build\bin\jrclust.exe sort "F:\...\amplifier_denoised.prm" --save --veto-spikes --rm-outliers
```

Writes `_filt.jrc`, `_raw.jrc`, `_features.jrc` and `_cpp.bin` to the `.prm`'s
`outputDir`. Add `--output-dir <path>` to send them elsewhere and keep an
existing result.

`--veto-spikes` rejects events that cannot be action potentials before they
reach clustering. It removes spikes but *raises* the number of well isolated
units, so it is not discarding neurons.

## 2. Import

```matlab
jrc_import_cpp(prm, 'autoMerge', false)
```

Builds `_res.mat` from `_cpp.bin` using JRCLUST's own clustering object, and
applies the C++ classifier's `noise` / `single` / `multi` labels.

**Pass `'autoMerge', false`.** The C++ already merged to convergence during
the sort, so MATLAB would only repeat it - and re-merging renumbers every
cluster, which scrambles unit ids before you have looked at anything.

This step is how you get back to a clean slate without re-sorting: it rebuilds
the result from the sorter's own output, discarding any labelling you have
done since.

## 3. Label

```matlab
jrc_relabel(prm, 'dryRun', true)    % preview, writes nothing
jrc_relabel(prm)                    % apply
```

Writes two annotations and touches nothing else:

- `to_delete` - fails a test decisively
- `review` - fails within 1.5x of a cut, so it is kept but flagged

It never moves a spike, so run it as often as you like with different
thresholds; each run clears its own previous labels and judges afresh.

The tests, each calibrated against units judged by eye rather than set from
theory:

| test | catches |
| --- | --- |
| `maxISIRatio` 0.30 | fires below the 1 ms refractory floor |
| `maxZeroLag` 3.0 | autocorrelogram peaks at zero lag instead of dipping |
| `minHalfWidthMs` 0.06 | trough one sample wide |
| `maxCrossings` 12 | waveform rings instead of spiking |
| `maxSNR` 50 | amplitude beyond anything biological |
| `loudSNR` / `loudMaxSpikes` | large amplitude with almost no spikes |

A template-match score is also reported but condemns nothing by default - see
*Template matching* below.

## 4. Look at it

```matlab
jrc manual prm
```

`View > Display annotated units only` shows the flagged ones. Start with
`review`: those are the near-threshold cases, deliberately not deleted,
and the ones where your judgement beats the statistics.

**Never use `Edit > Delete annotated`.** Use step 5 instead - see *Things that
will bite you*.

## 5. Prune

Only `to_delete`:

```matlab
jrc_prune(prm)
```

Or `to_delete` **and** the C++ `noise` labels:

```matlab
jrc_mark_noise(prm, 'from', {'noise'})
jrc_prune(prm)
```

Writes a reduced copy to `<prm folder>\pruned\`, leaving the original intact.
`File > Save As Pruned` in the GUI does the same thing.

The second form removes more but costs you some real cells: on a measured
recording it removed 105 extra units, of which 47 were genuine glitches and 58
had perfectly normal waveforms. Prefer the first form and delete the `noise`
units by hand once you have seen them.

---

## Optional steps

These restructure the clustering, so they are **not** free to re-run and your
existing annotations refer to the old unit ids afterwards. Back up the
`_res.mat` first - it is the only file they modify.

### Merge over-split units

```matlab
jrc_automerge(prm, 'dryRun', true)
jrc_automerge(prm)
```

Proposes pairs by mean-waveform correlation, as `Edit > Merge auto` does, then
accepts them only if the cross-correlogram between the two spike trains shows
a shared refractory dip.

Waveform shape alone cannot tell one neuron split in two from two neurons the
same distance from the same sites - their mean waveforms are equally alike
either way. On an 841-unit recording a correlation cut of 0.95 proposed 200
pairs and only **91** passed the check; the other 109 were independent cells
that a plain merge would have fused. A split you can merge later; a merge you
cannot undo.

Also available as `Edit > Merge auto (cross-correlogram checked)`.

### Split contaminating spikes out of a unit

```matlab
jrc_split_noise(prm, 'dryRun', true)
jrc_split_noise(prm)
```

Scores every individual spike against a template and moves the ones that do
not look like action potentials into a **new cluster immediately after its
parent**, annotated `splitted_to_noise`. Unit ids above each split shift up by
one, so parent and sibling sit adjacent in the GUI and you compare them
directly.

Nothing is deleted - you decide per pair whether to delete the sibling, merge
it back, or keep both.

`'madFactor'` is the knob: 6 splits less, 3 splits more.

### Repair derived data

```matlab
jrc_recompute(prm)
```

Rebuilds mean waveforms, waveform similarity, centroids and quality scores
from the current spike table. `jrc_automerge` and `jrc_split_noise` call it
themselves; run it by hand if you have moved spikes some other way and the
similarity matrix looks wrong or the space key stops finding the nearest unit.

### Repair inconsistent fields

```matlab
jrc_fix_fields(prm, 'dryRun', true)
```

If the GUI throws `Index exceeds the number of array elements` while deleting,
a per-cluster field has drifted out of step with the others. This brings them
all back to the same length.

---

## Template matching

`jrc_relabel` reports a score for how spike-shaped each unit is, derived from
the recording itself with **no labels**:

1. Normalise every unit's mean waveform to its own trough depth, leaving only
   shape, and take the per-sample median of the aligned ones. Real units
   outnumber noise by better than ten to one, so the median is theirs.
2. Check that candidate against membrane biophysics rather than against
   examples: one dominant excursion, a duration no membrane can beat,
   repolarisation, a return to baseline. If it fails, the recording is
   majority noise and the score is refused rather than inverted.

Because every recording builds its own template, nothing transfers between
sessions or animals - which is exactly why it generalises. A different probe,
filter or preparation yields a different template, correctly.

**It condemns nothing by default.** One template assumes one waveform family,
and a legitimate but atypical cell scores low for being unusual rather than
for being fake: on a measured recording, units of 85,000 to 220,000 spikes at
SNR 6.7 to 12.9 scored 0.72 to 0.76 while hand-judged noise scored 0.25 to
0.77. Those overlap, so any cut would be guesswork. Enable it only once you
have judged where the boundary sits on your own data:

```matlab
jrc_relabel(prm, 'minTemplateScore', 0.85)
```

Scored blind against twelve hand-labelled units it got 12 of 12, so the idea
is sound - the threshold is what is unresolved.

---

## Things that will bite you

**`Edit > Delete annotated`** holds one complete copy of the clustering per
deleted unit, all at once - 514 MB each on a 1143-unit recording, so 120 units
asks for 60 GB. The patched JRCLUST in `matlab/JRCLUST5.0` fixes this, but
`jrc_prune` is still the better route because it also shrinks the `.jrc`
files; `Delete annotated` only rewrites the spike table, so every discarded
spike stays on disk.

**`jrc sort` or `jrc recluster` on a curated result** calls `reassign()`,
which resets every label and re-clusters from scratch. It undoes all of the
above. To re-sort, go back to step 1.

**Anything that moves spikes** - merging, splitting - makes the stored mean
waveforms, similarity matrix and quality scores describe a clustering that no
longer exists, and the GUI reads them from the file rather than regenerating
them. Run `jrc_recompute` afterwards; the scripts here do it for you.

**A `.prm` that will not load**, complaining about `figList`, was written by a
JRCLUST build whose parameter list differs from the one reading it. Those are
only remembered figure positions - delete both lines:

```matlab
L = readlines(prm);
L(startsWith(strtrim(L), {'figList','figPos'})) = [];
fid = fopen(prm,'w'); for k=1:numel(L), fprintf(fid,'%s\n',L(k)); end; fclose(fid);
```

**Zero-filled recording gaps** read as enormous simultaneous spikes on every
site. The sorter rejects them by default (`rejectGaps`), finding them from the
data rather than from an acquisition log.
