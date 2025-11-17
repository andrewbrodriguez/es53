%% Figure 1: Exercise 1, 3 subjects, 3 segments
subjectFiles = { ...
    'Exercise1_Andrew_lab6.mat', ...
    'Exercise1_Elena_lab6.mat', ...
    'Exercise1_Evelyn_lab6.mat' ...
    };

subjectLabels = {'Andrew','Elena','Evelyn'};

figure(1); clf;
set(gcf,'Name','Figure 1 - Respiration','Color','w');

for s = 1:numel(subjectFiles)
    S = load(subjectFiles{s});

    % --- basic time setup ---
    fs   = S.samplerate(1,1);          % 100 Hz in these labs
    N    = numel(S.data);
    tAll = (0:N-1)/fs;

    % --- rebuild respiration channel (channel 1) ---
    resp = nan(1,N);
    nBlocks = size(S.datastart, 2);
    for b = 1:nBlocks
        i1 = S.datastart(1,b);
        i2 = S.dataend(1,b);
        resp(i1:i2) = S.data(i1:i2);
    end

    % --- comments: in your files, sample index is column 3 ---
    comSamples = round(S.com(:,3));    % 5×1
    % fix the "0" case (Andrew): time 0 → index 1
    comSamples(comSamples < 1) = 1;
    % just in case anything goes past the end
    comSamples(comSamples > N) = N;

    if numel(comSamples) < 5
        error('File %s still doesn''t have 5 usable comment indices.', subjectFiles{s});
    end

    % convert to times
    comTimes = tAll(comSamples);

    % your scheme:
    % 1→2 = baseline
    seg1_t1 = comTimes(1); seg1_t2 = comTimes(2);
    % 2→3 = breath hold
    seg2_t1 = comTimes(2); seg2_t2 = comTimes(3);
    % 4→5 = exhale
    seg3_t1 = comTimes(4); seg3_t2 = comTimes(5);

    %% ---------- row 1: baseline ----------
    row = 1;
    mask = tAll >= seg1_t1 & tAll <= seg1_t2 - 10;
    subplot(3,3,(row-1)*3 + s);
    plot(tAll(mask) - seg1_t1, resp(mask), 'k', 'LineWidth', 1);
    if row == 1
        title(subjectLabels{s});
    end
    ylabel('Voltage (mV)');
    title('Baseline', subjectLabels{s});
    box off;

    %% ---------- row 2: breath hold (2→3) ----------
    row = 2;
    mask = tAll >= seg2_t1 & tAll <= seg2_t2;
    subplot(3,3,(row-1)*3 + s);
    plot(tAll(mask) - seg2_t1, resp(mask), 'k', 'LineWidth', 1);
    ylabel('Voltage (mV)');
    hold on;
    title('Inhale and Hold', subjectLabels{s});
    hold off;
    box off;

    %% ---------- row 3: exhale (4→5) ----------
    row = 3;
    mask = tAll >= seg3_t1 & tAll <= seg3_t2;
    subplot(3,3,(row-1)*3 + s);
    plot(tAll(mask) - seg3_t1, resp(mask), 'k', 'LineWidth', 1);
    ylabel('Voltage (mV)');
    xlabel('Time (s)');
    hold on;
    title('Exhale and Hold', subjectLabels{s});
    hold off;
    box off;
end
