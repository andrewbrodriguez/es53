%% Lab 5 – Simple PIF / PEF + Peak Volume (LabChart blocks-safe)
clear; clc;

files = { ...
    'Exercise3-Andrew-Lab5.mat', ...
    'Exercise3-Elena-Lab5.mat', ...
    'Exercise3-Evelyn-Lab5.mat' ...
    'Exercise4-Andrew-Lab5.mat', ...
};

% ---- Tunables ----
minPeakPromAbsFlow = 0.8;   % L/s
minPeakDistanceSec = 0.5;   % s

for f = 1:numel(files)
    S = load(files{f});

    % -------- Flow (channel 1) --------
    ch = 1;                                 % Flow channel
    [~, nbl] = size(S.datastart);
    flow = []; fs = [];
    for bl = 1:nbl
        if S.datastart(ch,bl) == -1, continue; end
        i1 = S.datastart(ch,bl); i2 = S.dataend(ch,bl);
        seg = double(S.data(i1:i2));
        if isfield(S,'scaleunits') && isfield(S,'scaleoffset')
            seg = (seg + double(S.scaleoffset(ch,bl))).*double(S.scaleunits(ch,bl));
        end
        flow = [flow; seg]; %#ok<AGROW>
        if isempty(fs), fs = double(S.samplerate(ch,bl)); end
    end
    if isempty(flow) || numel(flow) < 3
        fprintf('\n%s\n  (No usable samples on channel 1 — skipping)\n', files{f});
        continue;
    end
    t  = (0:numel(flow)-1)'/fs;
    md = max(1, round(minPeakDistanceSec * fs));

    % PIF (positive) & PEF (negative)
    [pifVals, pifIdx] = findpeaks(flow, ...
        'MinPeakProminence', minPeakPromAbsFlow, ...
        'MinPeakDistance',   md, ...
        'SortStr','descend','NPeaks',3);

    [pefVals, pefIdx] = findpeaks(-flow, ...
        'MinPeakProminence', minPeakPromAbsFlow, ...
        'MinPeakDistance',   md, ...
        'SortStr','descend','NPeaks',3);
    pefVals = -pefVals;

    fprintf('\n=== %s ===\n', files{f});
    for i = 1:numel(pifVals)
        fprintf('  PIF #%d = %.3f L/s (%.1f L/min) at t = %.2f s\n', ...
            i, pifVals(i), 60*pifVals(i), t(pifIdx(i)));
    end
    for i = 1:numel(pefVals)
        fprintf('  PEF #%d = %.3f L/s (%.1f L/min) at t = %.2f s\n', ...
            i, abs(pefVals(i)), 60*abs(pefVals(i)), t(pefIdx(i)));
    end

    % -------- Volume (channel 2) --------
    chV = 2;                                % Volume channel
    vol = []; fsV = [];
    for bl = 1:nbl
        if S.datastart(chV,bl) == -1, continue; end
        i1 = S.datastart(chV,bl); i2 = S.dataend(chV,bl);
        seg = double(S.data(i1:i2));
        if isfield(S,'scaleunits') && isfield(S,'scaleoffset')
            seg = (seg + double(S.scaleoffset(chV,bl))).*double(S.scaleunits(chV,bl));
        end
        vol = [vol; seg]; %#ok<AGROW>
        if isempty(fsV), fsV = double(S.samplerate(chV,bl)); end
    end
    if ~isempty(vol)
        tV = (0:numel(vol)-1)'/fsV;
        [vMax, iMax] = max(vol);
        fprintf('  Peak Volume = %.3f L at t = %.2f s\n', vMax, tV(iMax));
    else
        fprintf('  (No usable samples on channel 2 for volume)\n');
    end
end

%%

%% Lab 5 — FVC from Volume using 3 forced breaths (anchor on PIF from Flow)

files = { ...
    'Exercise3-Andrew-Lab5.mat', ...
    'Exercise3-Elena-Lab5.mat', ...
    'Exercise3-Evelyn-Lab5.mat' ...
    'Exercise4-Andrew-Lab5.mat', ...
};

% ---- Tunables ----
minPeakPromAbsFlow = 0.8;   % L/s  (flow prominence for forced breaths)
minPeakDistanceSec = 0.5;   % s    (min spacing between PIFs)
windowSec          = 10.0;  % s    (analysis window length around each PIF)

for f = 1:numel(files)
    S = load(files{f});

    % ---------- stitch Flow (ch 1) across blocks ----------
    chF = 1;
    [~, nbl] = size(S.datastart);
    flow = []; fsF = [];
    for bl = 1:nbl
        if S.datastart(chF,bl) == -1, continue; end
        i1 = S.datastart(chF,bl); i2 = S.dataend(chF,bl);
        seg = double(S.data(i1:i2));
        if isfield(S,'scaleunits') && isfield(S,'scaleoffset')
            seg = (seg + double(S.scaleoffset(chF,bl))).*double(S.scaleunits(chF,bl));
        end
        flow = [flow; seg]; %#ok<AGROW>
        if isempty(fsF), fsF = double(S.samplerate(chF,bl)); end
    end
    if numel(flow) < 3
        fprintf('\n%s\n  (No usable flow samples — skipping)\n', files{f});
        continue;
    end
    tF = (0:numel(flow)-1)'/fsF;

    % --- find 3 largest PIF peaks on Flow to anchor forced breaths ---
    minDist = max(1, round(minPeakDistanceSec * fsF));
    [pifVals, pifIdx] = findpeaks(flow, ...
        'MinPeakProminence', minPeakPromAbsFlow, ...
        'MinPeakDistance',   minDist, ...
        'SortStr','descend', 'NPeaks', 3);
    pifIdx = sort(pifIdx); % analyze in time order
    tPIF   = tF(pifIdx);

    % ---------- stitch Volume (ch 2) across blocks ----------
    chV = 2;
    vol = []; fsV = [];
    for bl = 1:nbl
        if S.datastart(chV,bl) == -1, continue; end
        i1 = S.datastart(chV,bl); i2 = S.dataend(chV,bl);
        seg = double(S.data(i1:i2));
        if isfield(S,'scaleunits') && isfield(S,'scaleoffset')
            seg = (seg + double(S.scaleoffset(chV,bl))).*double(S.scaleunits(chV,bl));
        end
        vol = [vol; seg]; %#ok<AGROW>
        if isempty(fsV), fsV = double(S.samplerate(chV,bl)); end
    end
    if numel(vol) < 3
        fprintf('\n%s\n  (No usable volume samples)\n', files{f});
        continue;
    end
    tV = (0:numel(vol)-1)'/fsV;

    % --- compute FVC for each forced breath anchored at tPIF ---
    FVCs = nan(1,3);
    for k = 1:min(3, numel(tPIF))
        tC = tPIF(k);
        tL = max(tV(1),  tC - windowSec/2);
        tR = min(tV(end), tC + windowSec/2);
        idx = (tV >= tL) & (tV <= tR);
        if nnz(idx) < 3
            % minimal fallback: nearest 1 s span around tC
            [~, iC] = min(abs(tV - tC));
            halfW = max(1, round(0.5*fsV));
            iL = max(1, iC-halfW); iR = min(numel(tV), iC+halfW);
            idx = false(size(tV)); idx(iL:iR) = true;
        end

        vseg = vol(idx); tseg = tV(idx);
        [vPeak, iPeak] = max(vseg);            % peak inhalation (volume max)
        if iPeak < numel(vseg)
            vMinAfter = min(vseg(iPeak:end));  % maximal expiration after the peak
        else
            vMinAfter = vseg(end);
        end
        FVCs(k) = vPeak - vMinAfter;
    end

    % ---- print results ----
    fprintf('\n=== %s ===\n', files{f});
    for k = 1:3
        if isfinite(FVCs(k))
            fprintf('  FVC (breath %d) = %.3f L\n', k, FVCs(k));
        else
            fprintf('  FVC (breath %d) = (not found)\n', k);
        end
    end
    fprintf('  FVC mean (of found breaths) = %.3f L\n', mean(FVCs,'omitnan'));
end

%%
%% Lab 5 — FEV1 for three forced breaths per subject (blocks-safe)
clear; clc;

files = { ...
    'Exercise3-Andrew-Lab5.mat', ...
    'Exercise3-Elena-Lab5.mat', ...
    'Exercise3-Evelyn-Lab5.mat' ...
    'Exercise4-Andrew-Lab5.mat', ...
};

% ---- Tunables ----
minPeakPromAbsFlow = 0.8;   % L/s (flow prominence for forced breaths)
minPeakDistanceSec = 0.5;   % s   (min spacing between PIFs)
windowSec          = 10.0;  % s   (volume analysis window centered on PIF)

for f = 1:numel(files)
    S = load(files{f});

    % ---------- stitch Flow (channel 1) across blocks ----------
    chF = 1;
    [~, nbl] = size(S.datastart);
    flow = []; fsF = [];
    for bl = 1:nbl
        if S.datastart(chF,bl) == -1, continue; end
        i1 = S.datastart(chF,bl); i2 = S.dataend(chF,bl);
        seg = double(S.data(i1:i2));
        if isfield(S,'scaleunits') && isfield(S,'scaleoffset')
            seg = (seg + double(S.scaleoffset(chF,bl))).*double(S.scaleunits(chF,bl));
        end
        flow = [flow; seg]; %#ok<AGROW>
        if isempty(fsF), fsF = double(S.samplerate(chF,bl)); end
    end
    if numel(flow) < 3
        fprintf('\n%s\n  (No usable flow samples — skipping)\n', files{f});
        continue;
    end
    tF = (0:numel(flow)-1)'/fsF;

    % --- find up to 3 largest PIF peaks to anchor the three maneuvers ---
    minDist = max(1, round(minPeakDistanceSec * fsF));
    [~, pifIdx] = findpeaks(flow, ...
        'MinPeakProminence', minPeakPromAbsFlow, ...
        'MinPeakDistance',   minDist, ...
        'SortStr','descend', 'NPeaks', 3);
    pifIdx = sort(pifIdx);            % analyze in time order
    tPIF   = tF(pifIdx);

    % ---------- stitch Volume (channel 2) across blocks ----------
    chV = 2;
    vol = []; fsV = [];
    for bl = 1:nbl
        if S.datastart(chV,bl) == -1, continue; end
        i1 = S.datastart(chV,bl); i2 = S.dataend(chV,bl);
        seg = double(S.data(i1:i2));
        if isfield(S,'scaleunits') && isfield(S,'scaleoffset')
            seg = (seg + double(S.scaleoffset(chV,bl))).*double(S.scaleunits(chV,bl));
        end
        vol = [vol; seg]; %#ok<AGROW>
        if isempty(fsV), fsV = double(S.samplerate(chV,bl)); end
    end
    if numel(vol) < 3
        fprintf('\n%s\n  (No usable volume samples)\n', files{f});
        continue;
    end
    tV = (0:numel(vol)-1)'/fsV;

    % --- compute FEV1 for each forced breath ---
    FEV1s = nan(1,3);
    for k = 1:min(3, numel(tPIF))
        tC = tPIF(k);
        % analyze volume within a window centered on PIF
        tL = max(tV(1),  tC - windowSec/2);
        tR = min(tV(end), tC + windowSec/2);
        idx = (tV >= tL) & (tV <= tR);
        if nnz(idx) < 3
            % fallback: 1 s span around anchor if window too small
            [~, iC] = min(abs(tV - tC));
            iL = max(1, iC - round(0.5*fsV));
            iR = min(numel(tV), iC + round(0.5*fsV));
            idx = false(size(tV)); idx(iL:iR) = true;
        end
        vseg = vol(idx); tseg = tV(idx);

        % 1) volume peak (peak inhalation) within window
        [vPeak, iPeak] = max(vseg);
        t0 = tseg(iPeak);

        % 2) volume exactly 1.0 s after the peak (linear interpolation)
        t1 = t0 + 2.5;
        if tseg(1) <= t1 && tseg(end) >= t1 && numel(tseg) >= 2
            v1s = interp1(tseg, vseg, t1, 'linear');
        else
            % outside window: try using the full series if it brackets t1
            if tV(1) <= t1 && tV(end) >= t1 && numel(tV) >= 2
                v1s = interp1(tV, vol, t1, 'linear');
            else
                % cannot evaluate 1s later — skip this breath
                FEV1s(k) = NaN; continue;
            end
        end

        % 3) FEV1 = drop over 1 s after peak
        FEV1s(k) = vPeak - v1s;
    end

    % ---- print results ----
    fprintf('\n=== %s ===\n', files{f});
    for k = 1:3
        if isfinite(FEV1s(k))
            fprintf('  FEV1 (breath %d) = %.3f L\n', k, FEV1s(k));
        else
            fprintf('  FEV1 (breath %d) = (not found)\n', k);
        end
    end
    fprintf('  FEV1 mean (of found breaths) = %.3f L\n', mean(FEV1s,'omitnan'));
end