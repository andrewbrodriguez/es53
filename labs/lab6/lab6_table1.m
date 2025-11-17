%% Table 1 numbers from Exercise 1 files (stricter breath detection)
subjectFiles = { ...
    'Exercise1_Andrew_lab6.mat', ...
    'Exercise1_Elena_lab6.mat', ...
    'Exercise1_Evelyn_lab6.mat' ...
    };

subjectNames = {'Andrew','Elena','Evelyn'};

nSubj = numel(subjectFiles);

avgRate_bpm  = nan(nSubj,1);
inhaleHold_s = nan(nSubj,1);
exhaleHold_s = nan(nSubj,1);

for s = 1:nSubj
    S = load(subjectFiles{s});

    fs   = S.samplerate(1,1);
    N    = numel(S.data);
    tAll = (0:N-1)/fs;

    % rebuild resp channel (ch 1)
    resp = nan(1,N);
    for b = 1:size(S.datastart,2)
        i1 = S.datastart(1,b);
        i2 = S.dataend(1,b);
        resp(i1:i2) = S.data(i1:i2);
    end

    % comments (col 3), fix zeros
    comSamples = round(S.com(:,3));
    comSamples(comSamples < 1) = 1;
    comSamples(comSamples > N) = N;

    if numel(comSamples) < 5
        error('Need 5 comments for %s', subjectFiles{s});
    end

    comTimes = tAll(comSamples);

    % segments
    t_base_start   = comTimes(1);
    t_base_end     = comTimes(2);
    t_inhold_start = comTimes(2);
    t_inhold_end   = comTimes(3);
    t_exhold_start = comTimes(4);
    t_exhold_end   = comTimes(5);

    %% 1) average breathing rate during baseline
    baseMask = tAll >= t_base_start & tAll <= t_base_end;
    tBase    = tAll(baseMask) - t_base_start;
    yBase    = resp(baseMask);

    % smooth to kill tiny wiggles
    yBase_s  = smoothdata(yBase, 'movmean', round(0.4*fs));  % 0.4 s window

    % find peaks conservatively:
    % - at least 2.5 s apart  -> max ~24 bpm
    % - prominence = 20% of range
    yRange = max(yBase_s) - min(yBase_s);
    minProm = 0.1 * yRange;
    if minProm == 0
        minProm = 0.05;  % tiny fallback
    end

    [pks, locs] = findpeaks(yBase_s, fs, ...
        'MinPeakDistance', 2.5, ...   % seconds
        'MinPeakProminence', minProm);

    nBreaths = numel(locs);
    baseDur  = t_base_end - t_base_start;   % seconds

    avgRate_bpm(s) = (nBreaths / baseDur) * 60;

    %% 2) inhale breath-hold duration (2 -> 3)
    inhaleHold_s(s) = t_inhold_end - t_inhold_start;

    %% 3) exhale breath-hold duration (4 -> 5)
    exhaleHold_s(s) = t_exhold_end - t_exhold_start;
end

T = table(subjectNames(:), avgRate_bpm, inhaleHold_s, exhaleHold_s, ...
    'VariableNames', {'Subject','AvgBreathingRate_bpm','InhaleHold_s','ExhaleHold_s'});

disp(T);