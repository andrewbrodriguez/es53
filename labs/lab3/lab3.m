% File Name: lab2_arodriguez.m
% Author: Andrew Rodriguez
% Created: October 9th 2025
% Description: Lab 3 Figure 1 / 5 Table 3

clear; clc;


%% Figure 1

matFile = 'Lab3_Ex1_Defne.mat';
S = load(matFile);

% Normalize strings
if ischar(S.comtext),   S.comtext = cellstr(S.comtext); end
if ischar(S.titles),    S.titles  = cellstr(S.titles);  end

data = S.data(:);
datastart = S.datastart(:);
dataend = S.dataend(:); 
samplerate = S.samplerate(:);
tickrate = S.tickrate;
comtext = S.comtext(:);
com = S.com;

reqChan = 2;
chanIdx = reqChan;
hasData = (chanIdx <= numel(datastart)) && datastart(chanIdx) > 0 && dataend(chanIdx) > 0;
hasFs = chanIdx <= numel(samplerate) && samplerate(chanIdx) > 0;

if ~hasData || ~hasFs
    ecgIdx = find(contains(upper(S.titles),'ECG'), 1, 'first');
    if ~isempty(ecgIdx)
        chanIdx = ecgIdx;
    else
        valid = find(datastart>0 & dataend>0, 1, 'first');
        if ~isempty(valid), chanIdx = valid; end
    end
end

i0 = datastart(chanIdx);
i1 = dataend(chanIdx);
sig = data(i0:i1);
fs  = samplerate(chanIdx); if fs<=0, fs = tickrate; end

getCommentTime = @(label) local_comment_time(label, comtext, com, chanIdx, tickrate);

t_hand = getCommentTime('Hand movement');
t_jump = getCommentTime('Jumping');

lenA = max(0, floor(t_hand * fs));
segA = sig(1 : min(lenA, numel(sig)));

iJump = max(1, 1 + floor(t_jump * fs));
availB = numel(sig) - (iJump - 1);
lenB   = min(lenA, availB);
segB   = sig(iJump : iJump + lenB - 1);

tA = (0:numel(segA)-1) / fs;
tB = (0:numel(segB)-1) / fs;

figure('Color','w'); clf
plot(tA, segA, 'DisplayName','Still'); hold on
plot(tB, segB, 'DisplayName','Jumping Jacks');
grid on
xlabel('Time (s)');
ylabel('Voltage (V)');
legend('Location','best');


saveas(gcf, 'Figure1.png');

function t = local_comment_time(label, comtext, com, chanIdx, tickrate)
    % label match is case/space-insensitive
    want = find(strcmpi(strtrim(comtext), strtrim(label)), 1, 'first');
    if isempty(want)
        error('Comment label "%s" not found in comtext.', label);
    end
    % com columns (LabChart): [channel, ?, tick_index, ?, text_index]
    rows = com(com(:,1)==chanIdx & com(:,5)==want, :);
    if isempty(rows)
        error('No "%s" comment found for channel %d.', label, chanIdx);
    end
    ticks = rows(1,3);    % use the first occurrence
    t = ticks / tickrate;
end

%% Table 3

names = ["andrew","jerry","defne","izzy"];
files = names + "_exercise2.mat";

rows = struct('name',[],'Amp_relaxed_mV',[],'Amp_standing_mV',[],'Amp_exercise_mV',[]);

for k = 1:numel(files)
    F = files(k);
    S = load(F);
    if ischar(S.titles),  S.titles  = cellstr(S.titles);  end
    titles = string(S.titles(:));

    data = S.data(:);
    datastart = S.datastart;
    dataend = S.dataend;
    samplerate = S.samplerate(:);
    tickrate = S.tickrate;

    who = lower(string(names(k)));
    % last valid channel for Defne; ECG else
    if any(strcmpi(who, ["defne","defnes"]))
        valid = find(all([datastart>0, dataend>0],2));
        if isempty(valid)
            error('%s: no valid channels.', F);
        end
        chanIdx = valid(end);
    else
        chanIdx = find(contains(upper(titles),'ECG'), 1, 'first');
        if isempty(chanIdx)
            v = [];
            if isvector(datastart), v = var(double(data(datastart(1):dataend(1)))); 
            else
                v = zeros(size(datastart,1),1);
                for c = 1:size(datastart,1)
                    mask = datastart(c,:)>0 & dataend(c,:)>0;
                    if any(mask)
                        i0 = datastart(c, find(mask,1,'first'));
                        i1 = dataend(c,   find(mask,1,'first'));
                        v(c) = var(double(data(i0:i1)));
                    end
                end
            end
            if ~isempty(v), [~, chanIdx] = max(v); else, chanIdx = 1; end
        end
    end


    if isvector(datastart)
        ds_row = datastart(:).';   de_row = dataend(:).';
    else
        ds_row = datastart(chanIdx, :);
        de_row = dataend(chanIdx, :);
    end
    mask = ds_row>0 & de_row>0;
    ds_row = ds_row(mask);  de_row = de_row(mask);

    fs = samplerate(min(chanIdx, numel(samplerate)));
    if fs<=0
        if isscalar(tickrate), fs = double(tickrate);
        else
            if chanIdx<=numel(tickrate) && tickrate(chanIdx)>0, fs = double(tickrate(chanIdx));
            else, tz = tickrate(tickrate>0); fs = double(~isempty(tz)*tz(1) + isempty(tz)*1000);
            end
        end
    end

    if numel(ds_row) >= 3
        seg_relax = double(data(ds_row(1):de_row(1)));
        seg_stand = double(data(ds_row(2):de_row(2)));
        seg_post  = double(data(ds_row(3):de_row(3)));
    else
        [seg_relax, seg_stand, seg_post] = segments_from_comments(S, chanIdx, fs);
    end

    % amplitudes (mV): average dominant positive wave
    rows(k).name            = names(k);
    rows(k).Amp_relaxed_mV  = mean_peak_amplitude_mV(seg_relax, fs);
    rows(k).Amp_standing_mV = mean_peak_amplitude_mV(seg_stand, fs);
    rows(k).Amp_exercise_mV = mean_peak_amplitude_mV(seg_post,  fs);
end

T3 = struct2table(rows)
writetable(T3, 'Table3.csv');
fprintf('Saved Table3.csv\n');


function amp_mV = mean_peak_amplitude_mV(x, fs)
    if isempty(x) || numel(x) < round(0.5*fs), amp_mV = NaN; return; end
    x = x(:) - median(x,'omitnan');

    prom0 = 0.4 / 1e3;
    mpd   = round(0.5*fs); 

    % Peaks for amplitude
    [pks,~] = findpeaks(x, 'MinPeakDistance', mpd, 'MinPeakProminence', prom0);
    if numel(pks) < 2
        [pks,~] = findpeaks(x, 'MinPeakDistance', mpd, 'MinPeakProminence', 0.2*std(x,'omitnan'));
    end
    if isempty(pks), amp_mV = NaN; return; end

    % Trim 10% tails; convert V → mV
    pks = sort(pks);
    n = numel(pks); cut = floor(0.1*n);
    core = pks(max(1,1+cut):min(n,n-cut));
    amp_mV = mean(core,'omitnan') * 1e3;
end

% Rare fallback for non-fragmented files: split via comments
function [segA, segB, segC] = segments_from_comments(S, chanIdx, fs)
    data = S.data(:);
    if isvector(S.datastart), i0 = S.datastart(1); i1 = S.dataend(1);
    else, i0 = S.datastart(chanIdx,1); i1 = S.dataend(chanIdx,1); end
    sig = double(data(i0:i1));
    T   = (numel(sig)-1)/fs;

    C = double(S.com);
    rows_chan = C(C(:,1)==chanIdx, :);
    if isempty(rows_chan), rows_chan = C; end

    tr = S.tickrate;
    if isscalar(tr), tick = double(tr);
    else
        if chanIdx<=numel(tr) && tr(chanIdx)>0, tick = double(tr(chanIdx));
        else, tz = tr(tr>0); tick = double(~isempty(tz)*tz(1) + isempty(tz)*fs);
        end
    end

    t = rows_chan(:,3)./tick;
    textCol  = min(5,size(rows_chan,2));
    idx = rows_chan(:,textCol);
    idx(idx<1)=1; idx(idx>numel(S.comtext))=numel(S.comtext);
    L = lower(strtrim(string(S.comtext(idx))));

    t_relax = 0;
    t_stand = first_or_nan(t(contains(L,"stand")));
    t_post  = first_or_nan(t(contains(L,["post","exerc"])));

    if isnan(t_stand), t_stand = T/3; end
    if isnan(t_post),  t_post  = 2*T/3; end

    segA = slice_by_time(sig, fs, t_relax, t_stand);
    segB = slice_by_time(sig, fs, t_stand, t_post);
    segC = slice_by_time(sig, fs, t_post,  T);
end

function y = slice_by_time(x, fs, t0, t1)
    i0 = max(1, 1 + floor(t0 * fs));
    i1 = max(i0, 1 + floor(t1 * fs));
    i1 = min(i1, numel(x));
    y  = x(i0:i1);
end
function out = first_or_nan(v), if isempty(v), out = NaN; else, out = v(1); end, end

%% Figure 5
names = ["andrew","jerry","defne","izzy"];
files = names + "_exercise2.mat";

positions = ["Relaxed","Standing","Post-Exercise"];
nS = numel(names); nP = numel(positions);
MEAN = nan(nP, nS);     % rows = positions, cols = subjects
SD   = nan(nP, nS);

for s = 1:nS
    S = load(files(s));

    % --- normalize ---
    if ischar(S.titles), S.titles = cellstr(S.titles); end
    titles = string(S.titles(:));
    data = S.data(:);
    datastart = S.datastart;  dataend = S.dataend;
    samplerate = S.samplerate(:);  tickrate = S.tickrate;

    % --- choose channel ---
    who = lower(string(names(s)));
    if any(strcmpi(who, ["defne","defnes"]))
        % bottom row = last valid channel
        if isvector(datastart)
            chanIdx = 1;
        else
            valid = find(any(datastart>0 & dataend>0,2));
            if isempty(valid), error('%s: no valid channels', files(s)); end
            chanIdx = valid(end);
        end
    else
        % prefer ECG else fallback
        chanIdx = find(contains(upper(titles),'ECG'), 1, 'first');
        if isempty(chanIdx), chanIdx = 1; end
    end

    % --- fs ---
    fs = samplerate(min(chanIdx, numel(samplerate)));
    if fs<=0
        if isscalar(tickrate), fs = double(tickrate);
        else
            if chanIdx<=numel(tickrate) && tickrate(chanIdx)>0, fs = double(tickrate(chanIdx));
            else, tz = tickrate(tickrate>0); fs = double(~isempty(tz)*tz(1) + isempty(tz)*1000);
            end
        end
    end

    % --- get the 3 fragments for this channel ---
    if isvector(datastart)
        ds_row = datastart(:).'; de_row = dataend(:).';
    else
        ds_row = datastart(chanIdx,:); de_row = dataend(chanIdx,:);
    end
    mask = ds_row>0 & de_row>0;
    ds_row = ds_row(mask); de_row = de_row(mask);

    if numel(ds_row) < 3
        error('%s: expected 3 stop/start segments for Exercise 2.', files(s));
    end

    seg{1} = double(data(ds_row(1):de_row(1)));   % Relaxed
    seg{2} = double(data(ds_row(2):de_row(2)));   % Standing
    seg{3} = double(data(ds_row(3):de_row(3)));   % Post-Exercise

    % --- compute per-beat peak amplitudes (mV), then mean & SD ---
    for p = 1:3
        amps = peak_amps_mV(seg{p}, fs);          % vector (mV)
        if ~isempty(amps)
            % robust: trim 10% tails for mean/SD
            amps = sort(amps); n = numel(amps); cut = floor(0.1*n);
            core = amps(max(1,1+cut):min(n,n-cut));
            MEAN(p,s) = mean(core,'omitnan');
            SD(p,s)   = std(core, 'omitnan');
        end
    end
end

%% ----- plot (grouped by position, color by subject) -----
figure('Color','w'); clf
cats = categorical(positions); cats = reordercats(cats, positions);
hb = bar(cats, MEAN, 'grouped'); hold on
% error bars
for s = 1:nS
    x = hb(s).XEndPoints;
    errorbar(x, MEAN(:,s), SD(:,s), 'k','linestyle','none','LineWidth',1);
end
ylabel('Amplitude (mV)');
legend(cellstr(names), 'Location','northeastoutside');
grid on
% save (optional)
% print(gcf,'Figure5_bar_mean_sd.png','-dpng','-r300');
% print(gcf,'Figure5_bar_mean_sd.pdf','-dpdf');

%% ===== helper: per-beat peak amplitudes in mV (R, or P if larger) =====
function amps_mV = peak_amps_mV(x, fs)
    amps_mV = [];
    if isempty(x) || numel(x) < round(0.5*fs), return; end
    x = x(:) - median(x,'omitnan');               % baseline remove


    prom0 = 0.4 / 1e3;
    mpd   = round(0.5*fs); 


    % peaks = per-beat amplitudes (from baseline), convert to mV
    [pks,~] = findpeaks(x, 'MinPeakDistance', mpd, 'MinPeakProminence', prom0);
    if numel(pks) < 2
        [pks,~] = findpeaks(x, 'MinPeakDistance', mpd, 'MinPeakProminence', 0.2*std(x,'omitnan'));
    end
    amps_mV = pks * 1e3;                           % V -> mV
end