%% Table 2 — R–R intervals and HR (STANDING, robust + cleaned)
clear; clc; close all;

% Put your STANDING file names here (swap if different)
files    = {'defne_exercise2.mat','jerry_exercise2.mat','andrew_exercise2.mat','izzy_exercise2.mat'};
subjects = {'Defne','Jerry','Andrew','Izzy'};

n = numel(subjects);
meanRR = nan(1,n);  maxRR = nan(1,n);  minRR = nan(1,n);  stdRR = nan(1,n);  meanHR = nan(1,n);
kept   = zeros(1,n); removed = zeros(1,n);

fprintf('===== R–R INTERVAL ANALYSIS (Standing) =====\n\n');

for i = 1:n
    S = subjects{i}; F = files{i};
    fprintf('--- %s (%s) ---\n', S, F);

    % Load data
    D = load(F);

    % --- detect sampling rate (default 1000 Hz) ---
    fs = detect_fs(D, 1000);

    % --- pick an ECG vector from the file ---
    ecg = detect_ecg_vector(D);
    if isempty(ecg)
        warning('No numeric ECG-like vector found in %s; skipping.\n', F);
        continue;
    end

    % --- compute robust RR stats ---
    [stats, n_keep, n_drop] = rr_stats_from_ecg(ecg, fs);
    kept(i)    = n_keep;  removed(i) = n_drop;
    meanRR(i)  = stats.meanRR;
    maxRR(i)   = stats.maxRR;
    minRR(i)   = stats.minRR;
    stdRR(i)   = stats.stdRR;
    meanHR(i)  = stats.meanHR;

    fprintf('RR (ms): mean=%.2f  max=%.2f  min=%.2f  std=%.2f | mean HR=%.2f bpm | kept %d, removed %d\n\n', ...
        meanRR(i), maxRR(i), minRR(i), stdRR(i), meanHR(i), kept(i), removed(i));
end

% ----- build and display Table 2 -----
Table2 = table(meanRR', maxRR', minRR', stdRR', meanHR', ...
    'RowNames', subjects, ...
    'VariableNames', {'Mean_RR_ms','Max_RR_ms','Min_RR_ms','Std_RR_ms','Mean_HR_bpm'});
disp(Table2);
fprintf('\nTable 2 created successfully!\n');

% ================== helpers ==================
function fs = detect_fs(D, default_fs)
% Return sampling rate if present; otherwise default
    names = fieldnames(D);
    fs = default_fs;
    % Common scalar fields for fs
    for k = 1:numel(names)
        v = D.(names{k});
        if isscalar(v) && isnumeric(v) && v>50 && v<5000 && any(strcmpi(names{k}, ...
                {'fs','Fs','FS','sampling_rate','sample_rate','Hz','FsHz'}))
            fs = double(v); return;
        end
    end
    % Try to infer from a time vector
    for k = 1:numel(names)
        v = D.(names{k});
        if isnumeric(v) && isvector(v) && numel(v)>10
            dv = diff(v(:));
            if all(dv>0) && std(dv)/mean(dv) < 0.01
                est = 1/median(dv);
                if est>50 && est<5000, fs = est; return; end
            end
        end
    end
end

function x = detect_ecg_vector(D)
% Pick a reasonable ECG-like vector (long numeric vector or best column)
    names = fieldnames(D);
    x = [];
    bestScore = -inf;
    for k = 1:numel(names)
        v = D.(names{k});
        if isnumeric(v)
            if isvector(v) && numel(v) > 100
                vv = double(v(:));
                score = min(iqr(vv), std(vv))*1e3 / (1 + max(0,kurtosis(vv)-3)); % prefer dynamic but not spiky
                if score > bestScore, bestScore = score; x = vv; end
            elseif size(v,1) > size(v,2) && size(v,2) <= 4 && size(v,1) > 100
                [~,idx] = max(std(double(v),0,1));
                vv = double(v(:,idx));
                score = min(iqr(vv), std(vv))*1e3 / (1 + max(0,kurtosis(vv)-3));
                if score > bestScore, bestScore = score; x = vv(:); end
            end
        end
    end
end

function [stats, n_kept, n_removed] = rr_stats_from_ecg(ecg, fs, opts)
%RR_STATS_FROM_ECG  Robust RR-interval & mean-HR from a raw ECG trace.
%
% [stats, n_kept, n_removed] = rr_stats_from_ecg(ecg, fs, opts)
%   ecg : numeric vector (ECG, arbitrary units)
%   fs  : sampling rate in Hz
%   opts (optional) fields:
%       .maxBPM            (default 150)  % tighten to avoid T/T+noise double counts
%       .minBPM            (default 30)
%       .qrsBandHz         (default [5 18])  % bandpass around QRS if toolbox available
%       .qrsMaxWidthSec    (default 0.12)
%       .qrsMinWidthSec    (default 0.04)
%       .promScaleStd      (default 1.2)     % prominence = max( promScaleStd*std , promFracRange*range )
%       .promFracRange     (default 0.10)
%       .madK              (default 3)       % MAD threshold for outliers
%
% Outputs:
%   stats.meanRR (ms), stats.maxRR (ms), stats.minRR (ms), stats.stdRR (ms), stats.meanHR (bpm)
%   n_kept : # RR intervals retained after cleaning
%   n_removed : # RR intervals discarded (physiology clamp + MAD)

    if nargin < 3, opts = struct(); end
    maxBPM         = getdef(opts,'maxBPM',150);
    minBPM         = getdef(opts,'minBPM',30);
    qrsBandHz      = getdef(opts,'qrsBandHz',[5 18]);
    qrsMaxWidthSec = getdef(opts,'qrsMaxWidthSec',0.12);
    qrsMinWidthSec = getdef(opts,'qrsMinWidthSec',0.04);
    promScaleStd   = getdef(opts,'promScaleStd',1.2);
    promFracRange  = getdef(opts,'promFracRange',0.10);
    madK           = getdef(opts,'madK',3);

    % ---- sanitize input ----
    x = double(ecg(:));
    x = x(isfinite(x));
    if numel(x) < 5
        stats = empty_stats(); n_kept = 0; n_removed = 0; return;
    end

    % ---- baseline removal (works without toolboxes) ----
    win = max(1, round(0.7*fs));              % ~0.7 s moving-average high-pass
    if mod(win,2)==0, win = win+1; end
    x_detr = x - movmean(x, win, 'Endpoints','fill');

    % ---- optional bandpass around QRS (graceful fallback if no toolbox) ----
    f = x_detr;
    try
        [b,a] = butter(2, qrsBandHz/(fs/2), 'bandpass'); %#ok<*BDSCA>
        f = filtfilt(b,a,f);
    catch
        % fallback: light differentiator then smooth
        f = [0; diff(f)];
        f = movmean(f, round(0.04*fs));
    end

    % ---- choose polarity that yields stronger R peaks ----
    promGuess = max(promScaleStd*std(f), promFracRange*range(f));
    minDist   = round(fs * 60 / maxBPM);      % refractory: ≥ 1/maxBPM seconds
    [p1,~] = findpeaks( f, 'MinPeakDistance',minDist, 'MinPeakProminence',promGuess);
    [p2,~] = findpeaks(-f, 'MinPeakDistance',minDist, 'MinPeakProminence',promGuess);
    if sum(p2) > sum(p1), f = -f; end

    % ---- final peak detection (width guard ~QRS) ----
    args = { 'MinPeakDistance', minDist, ...
             'MinPeakProminence', promGuess, ...
             'MaxPeakWidth', round(qrsMaxWidthSec*fs), ...
             'MinPeakWidth', max(1,round(qrsMinWidthSec*fs)) };
    [~, locs] = findpeaks(f, args{:});

    if numel(locs) < 3
        stats = empty_stats(); n_kept = 0; n_removed = 0; return;
    end

    % ---- RR intervals (ms) ----
    RR_ms = diff(locs) * 1000 / fs;

    % ---- physiologic clamp first (30–150 bpm by default) ----
    lo_ms = 60000 / maxBPM;                    % e.g., 400 ms @ 150 bpm
    hi_ms = 60000 / minBPM;                    % e.g., 2000 ms @ 30 bpm
    phys_ok = (RR_ms >= lo_ms) & (RR_ms <= hi_ms);
    RR_ms = RR_ms(phys_ok);

    % ---- robust MAD outlier removal ----
    if isempty(RR_ms)
        stats = empty_stats(); n_kept = 0; n_removed = 0; return;
    end
    medRR = median(RR_ms);
    madv  = median(abs(RR_ms - medRR));
    if madv == 0, madv = 1e-6; end
    keep  = abs(RR_ms - medRR) <= madK*madv;
    n_removed = sum(~phys_ok) + sum(~keep);
    RR_ms = RR_ms(keep);
    n_kept = numel(RR_ms);

    if n_kept < 1
        stats = empty_stats(); return;
    end

    % ---- stats ----
    stats.meanRR = mean(RR_ms);
    stats.maxRR  = max(RR_ms);
    stats.minRR  = min(RR_ms);
    stats.stdRR  = std(RR_ms);
    stats.meanHR = mean(60000 ./ RR_ms);       % mean of instantaneous HR

end

% --------- helpers ----------
function v = getdef(s, field, default)
    if isfield(s, field) && ~isempty(s.(field))
        v = s.(field);
    else
        v = default;
    end
end

function s = empty_stats()
    s = struct('meanRR',nan,'maxRR',nan,'minRR',nan,'stdRR',nan,'meanHR',nan);
end