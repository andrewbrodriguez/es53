%% Lab 6 – Table 2 from channel 3 (shifted time), all 3 subjects

subjectFiles = { ...
    'Exercise1_Andrew_lab6.mat', ...
    'Exercise1_Elena_lab6.mat', ...
    'Exercise1_Evelyn_lab6.mat' ...
    };

subjectNames = {'Andrew','Elena','Evelyn'};

nSubj = numel(subjectFiles);

avg_norm = nan(nSubj,1);
avg_inh  = nan(nSubj,1);
avg_exh  = nan(nSubj,1);

for s = 1:nSubj
    S = load(subjectFiles{s});

    fs   = S.samplerate(1,1);
    N    = numel(S.data);
    tAll = (0:N-1)/fs;

    % ----- rebuild channel 3 for THIS subject -----
    ch3 = nan(1, N);
    for b = 1:size(S.datastart,2)
        i1 = S.datastart(3,b);
        i2 = S.dataend(3,b);
        if i1 > 0 && i2 > 0
            ch3(i1:i2) = S.data(i1:i2);
        end
    end

    % ----- when does ch3 actually start? -----
    iFirst3 = find(~isnan(ch3), 1, 'first');
    if isempty(iFirst3)
        % no ch3 data at all
        continue
    end
    tStart3 = tAll(iFirst3);
    % shifted time so ch3 "starts" at 0
    t3_shift = tAll - tStart3;

    % ----- get comment times (the desired windows) -----
    comSamples = round(S.com(:,3));
    comSamples(comSamples < 1) = 1;
    comSamples(comSamples > N) = N;
    comTimes = tAll(comSamples);

    % your scheme:
    % 1→2 = normal
    t_norm1 = comTimes(1);
    t_norm2 = comTimes(2);
    % 2→3 = inhale hold
    t_inh1  = comTimes(2);
    t_inh2  = comTimes(3);
    % 4→5 = exhale hold
    t_exh1  = comTimes(4);
    t_exh2  = comTimes(5);

    % ----- masks on *shifted* time -----
    mask_norm = t3_shift >= t_norm1 & t3_shift <= t_norm2;
    mask_inh  = t3_shift >= t_inh1  & t3_shift <= t_inh2;
    mask_exh  = t3_shift >= t_exh1  & t3_shift <= t_exh2;

    % ----- averages for this subject -----
    avg_norm(s) = mean(ch3(mask_norm), 'omitnan');
    avg_inh(s)  = mean(ch3(mask_inh),  'omitnan');
    avg_exh(s)  = mean(ch3(mask_exh),  'omitnan');
end

% put it in a table
T = table(subjectNames(:), avg_norm, avg_inh, avg_exh, ...
    'VariableNames', {'Subject','HR_Normal_bpm','HR_InhaleHold_bpm','HR_ExhaleHold_bpm'});

disp(T);