function [vals, vi] = subsample(vals, k, dim)
    %SUBSAMPLE Sample k items from vals, optionally along dim
    if nargin < 3
        dim = 2;
    end

    if isvector(vals)
        doReshape = isrow(vals);
        [vals, vi] = subsampleVec(vals(:), k);

        % reshape sampled values if necessary
        if doReshape
            vals = vals';
        end
    else
        if isempty(intersect(dim, 1:ndims(vals)))
            error('cannot sample along dimension %d', dim);
        end

        [vals, vi] = subsampleNd(vals, k, dim);
    end
end

%% LOCAL FUNCTIONS
function [vals, inds] = subsampleVec(vals, k)
    %SUBSAMPLEVEC Subsample `k` values from a vector
    if isempty(k) || k >= numel(vals)
        inds = 1:numel(vals);
    else
        inds = evenlySpaced(numel(vals), k);
    end

    vals = vals(inds);
end

function [vals, inds] = subsampleNd(vals, k, dim)
    %SUBSAMPLEND Subsample `k` rows (or cols, or ...) from an Nd array, N >= 2
    if isempty(k) || k >= size(vals, dim)
        inds = 1:size(vals, dim);
    else
        inds = evenlySpaced(size(vals, dim), k);
    end

    if dim == 1
        vals = vals(inds, :);
    else
        vals = vals(:, inds);
    end
end

function inds = evenlySpaced(n, k)
    %EVENLYSPACED k indices spread evenly over 1..n, deterministically.
    %
    %   Replaces randsample(), which drew a DIFFERENT subset on every call
    %   and made the whole pipeline irreproducible. subsample() feeds
    %   getPVSamples and spikePCA, so the principal vectors were re-estimated
    %   from a fresh random 10,000 spikes each time: the same recording
    %   loaded twice gave two different feature spaces. For a sort that
    %   concatenates sessions, that is the difference between features being
    %   comparable across sessions and not. It also feeds computeRho, so the
    %   clustering itself varied between runs on identical input.
    %
    %   Even spacing is not merely deterministic, it is the better estimator
    %   here. The inputs are ordered by time, so a regular stride draws from
    %   the whole recording in proportion, whereas a random draw is uniform
    %   only in expectation and can over-weight a high-rate stretch in any
    %   single run.
    %
    %   Only called when k < n, so the stride is at least 1 and the rounded
    %   values are strictly increasing; unique() is a guard, not a necessity.
    inds = unique(round(linspace(1, n, k)));
end
