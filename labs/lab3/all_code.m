% File Name: lab3_arodriguez.m
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


%% FIX FOR TABLE 2:

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


% Lab 3
% Isabelle Agarwal code
s = load("Lab3_Ex2_Defne.mat");
s_chunk = s.data(s.datastart(3,1):s.dataend(3,1));
startIndex = s.com(1,1);
endIndex = startIndex + 10*s.samplerate(3,1);
snip = s_chunk(startIndex:endIndex);
sTime = [0:length(snip)-1]/1000;
%Jerry
t = load("jerry_exercise2.mat");
t_chunk = t.data(t.datastart(1,1):t.dataend(1,1));
startIndex1 = t.samplerate(1,1);
endIndex1 = startIndex1 + 10*t.samplerate(1,1);
snip1 = t_chunk(startIndex1:endIndex1);
sTime1 = [0:length(snip1)-1]/1000;
%Izzy
x = load("izzy_exercise2.mat");
x_chunk = x.data(x.datastart(1,1):x.dataend(1,1));
startIndex2 = 10*x.samplerate(1,1);
endIndex2 = startIndex2 + 10*x.samplerate(1,1);
snip2 = x_chunk(startIndex2:endIndex2);
sTime2 = [0:length(snip2)-1]/1000;
%Andrew
z= load("andrew_exercise2.mat");
z_chunk = z.data(z.datastart(1,1):z.dataend(1,1));
startIndex3 = 10*z.samplerate(1,1);
endIndex3 = startIndex3 + 10*z.samplerate(1,1);
snip3 = z_chunk(startIndex3:endIndex3);
sTime3 = [0:length(snip3)-1]/1000;
%Full figure
figure;
hold on;
plot(sTime, snip);
plot(sTime1, snip1);
plot(sTime2, snip2);
plot(sTime3, snip3);
hold off;
legend("Defne", "Jerry", "Isabelle", "Andrew");
xlabel("Time(s)");
ylabel("Voltage (V)");
%% Figure 3
s = load("Lab3_Ex2_Defne.mat");
s_chunk = s.data(s.datastart(3,2):s.dataend(3,2));
startIndex = s.com(1,1);
endIndex = startIndex + 10*s.samplerate(3,2);
snip = s_chunk(startIndex:endIndex);
sTime = [0:length(snip)-1]/1000;
%Jerry
t = load("jerry_exercise2.mat");
t_chunk = t.data(t.datastart(1,2):t.dataend(1,2));
startIndex1 = t.samplerate(1,2);
endIndex1 = startIndex1 + 10*t.samplerate(1,2);
snip1 = t_chunk(startIndex1:endIndex1);
sTime1 = [0:length(snip1)-1]/1000;
%Izzy
x = load("izzy exercise2.mat");
x_chunk = x.data(x.datastart(1,2):x.dataend(1,2));
startIndex2 = 10*x.samplerate(1,2);
endIndex2 = startIndex2 + 10*x.samplerate(1,2);
snip2 = x_chunk(startIndex2:endIndex2);
sTime2 = [0:length(snip2)-1]/1000;
%Andrew
z= load("andrew_exercise2.mat");
z_chunk = z.data(z.datastart(1,2):z.dataend(1,2));
startIndex3 = 7*z.samplerate(1,2);
endIndex3 = startIndex3 + 10*z.samplerate(1,2);
snip3 = z_chunk(startIndex3:endIndex3);
sTime3 = [0:length(snip3)-1]/1000;
%Full figure
figure;
hold on;
plot(sTime, snip);
plot(sTime1, snip1);
plot(sTime2, snip2);
plot(sTime3, snip3);
hold off;
legend("Defne", "Jerry", "Isabelle", "Andrew");
xlabel("Time(s)");
ylabel("Voltage (V)");
%% Figure 6
s = load("Lab3_Ex3_Defne.mat");
s_chunk = s.data(s.datastart(3,1):s.dataend(3,1));
startIndex = 0.4*s.samplerate(3,1);
endIndex = startIndex + 15*s.samplerate(3,1);
snip = s_chunk(startIndex:endIndex);
sTime = [0:length(snip)-1]/1000;
%Jerry
t = load("jerry_exercise3.mat");
t_chunk = t.data(t.datastart(3,1):t.dataend(3,1));
startIndex1 = 0.4*t.samplerate(3,1);
endIndex1 = startIndex1 + 15*t.samplerate(3,1);
snip1 = t_chunk(startIndex1:endIndex1);
sTime1 = [0:length(snip1)-1]/1000;
%Izzy
x = load("izzy_exercise3.mat");
x_chunk = x.data(x.datastart(3,1):x.dataend(3,1));
startIndex2 = 0.4*x.samplerate(3,1);
endIndex2 = startIndex2 + 15*x.samplerate(3,1);
snip2 = x_chunk(startIndex2:endIndex2);
sTime2 = [0:length(snip2)-1]/1000;
%Andrew
z= load("andrew_exercise3.mat");
z_chunk = z.data(z.datastart(3,1):z.dataend(3,1));
startIndex3 = 0.4*z.samplerate(3,1);
endIndex3 = startIndex3 + 15*z.samplerate(3,1);
snip3 = z_chunk(startIndex3:endIndex3);
sTime3 = [0:length(snip3)-1]/1000;
%Full figure
s1 = load("Lab3_Ex3_Defne.mat");
s_chunk1 = s1.data(s1.datastart(2,1):s1.dataend(2,1));
startIndexs = 0.4*s1.samplerate(2,1);
endIndexs = startIndexs + 15*s1.samplerate(2,1);
snips = s_chunk1(startIndexs:endIndexs);
sTimes = [0:length(snips)-1]/1000;
%Jerry
t1 = load("jerry_exercise3.mat");
t_chunk1 = t1.data(t1.datastart(2,1):t1.dataend(2,1));
startIndext = 0.4*t1.samplerate(2,1);
endIndext = startIndext + 15*t1.samplerate(2,1);
snipt = t_chunk1(startIndext:endIndext);
sTimet = [0:length(snipt)-1]/1000;
%Izzy
x1 = load("izzy_exercise3.mat");
x_chunkx = x1.data(x1.datastart(2,1):x1.dataend(2,1));
startIndexx = 0.4*x1.samplerate(2,1);
endIndexx = startIndexx + 15*x1.samplerate(2,1);
snipx = x_chunkx(startIndexx:endIndexx);
sTimex = [0:length(snipx)-1]/1000;
%Andrew
z1= load("andrew_exercise3.mat");
z_chunkz = z1.data(z1.datastart(2,1):z1.dataend(2,1));
startIndexz = 0.4*z1.samplerate(2,1);
endIndexz = startIndexz + 15*z1.samplerate(2,1);
snipz = z_chunkz(startIndexz:endIndexz);
sTimez = [0:length(snipz)-1]/1000;
%Full figure
figure;
subplot(2,1,1)
hold on;
plot(sTime, snip);
plot(sTime1, snip1);
plot(sTime2, snip2);
plot(sTime3, snip3);
hold off;
legend("Defne", "Jerry", "Isabelle", "Andrew");
xlabel("Time(s)");
ylabel("EKG Amplitude (V)");
title("EKG Plot: All Group Members");
subplot(2,1,2)
hold on;
plot(sTimes, snips);
plot(sTimet, snipt);
plot(sTimex, snipx);
plot(sTimez, snipz);
hold off;
legend("Defne", "Jerry", "Isabelle", "Andrew");
xlabel("Time(s)");
ylabel("Transducer Output (V)");
title("Blood Volume Pulse: All Group Members");

% File Name: lab3_fig4_fig7_arodriguez.m
% Author: Defne Tiftik
% Created: October 9, 2025
% Description: Lab 3  –  Figure 4 (Defne, Exercise 2) and Figure 7 (All Subjects, Exercise 4)
clear; clc;

%% ------------------------------------------------------------------------
%  FIGURE 4 : Defne – ECG under different positions (Exercise 2)

matFile = 'defne_exercise2.mat';
S = load(matFile);

%normalize cell strings (LabChart quirk)
if ischar(S.titles), S.titles = cellstr(S.titles); end
titles = string(S.titles(:));
data  = S.data(:);
datastart = S.datastart;
dataend   = S.dataend;
samplerate = S.samplerate(:);
tickrate   = S.tickrate;

%choose ECG channel (prefer title containing "ECG")
chanIdx = find(contains(upper(titles),'ECG'),1,'first');
if isempty(chanIdx)
    valid = find(datastart>0 & dataend>0,1,'first');
    if isempty(valid), chanIdx = 1; else, chanIdx = valid; end
end

%sampling rate fallback logic
fs = samplerate(min(chanIdx,numel(samplerate)));
if fs <= 0
    if isscalar(tickrate)
        fs = double(tickrate);
    elseif chanIdx<=numel(tickrate) && tickrate(chanIdx)>0
        fs = double(tickrate(chanIdx));
    else
        tz = tickrate(tickrate>0);
        fs = double(~isempty(tz)*tz(1) + isempty(tz)*1000);
    end
end

%determine 3 posture segments (Lying / Sitting / Standing)
if isvector(datastart)
    ds_row = datastart(:).'; de_row = dataend(:).';
else
    ds_row = datastart(chanIdx,:); de_row = dataend(chanIdx,:);
end
mask = ds_row>0 & de_row>0;
ds_row = ds_row(mask); de_row = de_row(mask);

if numel(ds_row) < 3
    error('Expected 3 segments for lying/sitting/standing in %s', matFile);
end

seg{1} = double(data(ds_row(1):de_row(1)));   % Lying
seg{2} = double(data(ds_row(2):de_row(2)));   % Sitting
seg{3} = double(data(ds_row(3):de_row(3)));   % Standing

labels = {'Lying','Sitting','Standing'};
cols   = {'b','r','g'};

%convert time to milliseconds and trim to 10 000 ms (10 s)
for k = 1:3
    t{k} = (0:numel(seg{k})-1)/fs * 1000;     % ms
    idx = t{k} <= 10000;                      % keep first 10 s
    t{k} = t{k}(idx);
    seg{k} = seg{k}(idx);
end

% plot (wide figure)
figure('Color','w','Position',[100 100 1200 400]); clf
for k = 1:3
    plot(t{k}, seg{k}, 'Color', cols{k}, 'LineWidth',1.3, ...
         'DisplayName', labels{k}); hold on
end
xlabel('Time (ms)');
ylabel('ECG Amplitude (V)');
title('Figure 4: ECG of one subject (Defne – Lying, Sitting, Standing, first 10 s)');
legend('Location','northeastoutside');
grid on
xlim([0 10000]);

%% ------------------------------------------------------------------------
%  FIGURE 7  :  Blood Volume Pulse  –  All Subjects (Exercise 4)
names = ["andrew","izzy","jerry","defne"];
files = names + "_exercise4.mat";
cols  = {'b','r','g','m'};

figure('Color','w'); clf

% Above-head (first 15 s)
subplot(2,1,1); hold on
for k = 1:numel(files)
    F = files(k);
    if ~isfile(F), warning('Missing %s',F); continue; end
    S = load(F);
    if ischar(S.titles), S.titles = cellstr(S.titles); end
    data = S.data(:);
    datastart = S.datastart; dataend = S.dataend;
    samplerate = S.samplerate(:); tickrate = S.tickrate;

    % pick first valid pulse/ECG channel (same logic)
    chanIdx = find(contains(upper(string(S.titles)),'PULSE'),1,'first');
    if isempty(chanIdx), chanIdx = find(datastart>0 & dataend>0,1,'first'); end

    fs = samplerate(min(chanIdx,numel(samplerate)));
    if fs<=0
        if isscalar(tickrate), fs=double(tickrate);
        elseif chanIdx<=numel(tickrate) && tickrate(chanIdx)>0
            fs=double(tickrate(chanIdx));
        else
            tz = tickrate(tickrate>0);
            fs = double(~isempty(tz)*tz(1) + isempty(tz)*1000);
        end
    end

    i0 = datastart(chanIdx); i1 = dataend(chanIdx);
    sig = double(data(i0:i1));
    N = numel(sig); t = (0:N-1)/fs;
    segN = round(15*fs);
    idx = 1:min(segN,N);
    plot(t(idx), sig(idx), 'Color',cols{k}, 'LineWidth',1.2);
end
title('Figure 7A: Blood Volume Pulse (Hand Above Head)');
xlabel('Time (s)'); ylabel('Amplitude (V)');
legend(cellstr(names),'Location','northeastoutside');
grid on

% Below-hip (next 15 s or last 15 s)
subplot(2,1,2); hold on
for k = 1:numel(files)
    F = files(k);
    if ~isfile(F), continue; end
    S = load(F);
    data = S.data(:);
    datastart = S.datastart; dataend = S.dataend;
    samplerate = S.samplerate(:); tickrate = S.tickrate;

    chanIdx = find(contains(upper(string(S.titles)),'PULSE'),1,'first');
    if isempty(chanIdx), chanIdx = find(datastart>0 & dataend>0,1,'first'); end

    fs = samplerate(min(chanIdx,numel(samplerate)));
    if fs<=0
        if isscalar(tickrate), fs=double(tickrate);
        elseif chanIdx<=numel(tickrate) && tickrate(chanIdx)>0
            fs=double(tickrate(chanIdx));
        else
            tz = tickrate(tickrate>0);
            fs = double(~isempty(tz)*tz(1) + isempty(tz)*1000);
        end
    end

    i0 = datastart(chanIdx); i1 = dataend(chanIdx);
    sig = double(data(i0:i1));
    N = numel(sig); t = (0:N-1)/fs;
    segN = round(15*fs);
    if N >= 2*segN
        idx = (segN+1):(2*segN);
    else
        idx = max(1,N-segN+1):N;
    end
    plot(t(idx)-t(idx(1)), sig(idx), 'Color',cols{k}, 'LineWidth',1.2);
end
title('Figure 7B: Blood Volume Pulse (Hand Below Hip)');
xlabel('Time (s)'); ylabel('Amplitude (V)');
legend(cellstr(names),'Location','northeastoutside');
grid on

% File: lab3_table4.m
% Author: Defne Tiftik
clear; clc;

names = ["andrew","izzy","jerry","defne"];
files = ["andrew_exercise4.mat","izzy_exercise4.mat", ...
         "jerry_exercise4.mat","defne_exercise4.mat"];

rows = struct('Name',[],'AmpAbove_mV',[],'AmpBelow_mV',[], ...
              'HRAbove_bpm',[],'HRBelow_bpm',[]);

for k = 1:numel(files)
    F = files(k);
    if ~isfile(F)
        warning('File %s not found, skipping.', F);
        continue;
    end
    S = load(F);

    % pick a likely pulse channel 
    if ischar(S.titles), S.titles = cellstr(S.titles); end
    titles = string(S.titles(:));
    data = S.data(:);
    datastart = S.datastart;
    dataend   = S.dataend;
    samplerate = S.samplerate(:);
    tickrate   = S.tickrate;

    chanIdx = find(contains(upper(titles),'PULSE'),1,'first');
    if isempty(chanIdx)
        valid = find(datastart>0 & dataend>0,1,'first');
        if isempty(valid), chanIdx = 1; else, chanIdx = valid; end
    end

    % sampling rate (same fallback logic as your earlier code) 
    fs = samplerate(min(chanIdx,numel(samplerate)));
    if fs<=0
        if isscalar(tickrate), fs=double(tickrate);
        elseif chanIdx<=numel(tickrate) && tickrate(chanIdx)>0
            fs=double(tickrate(chanIdx));
        else
            tz=tickrate(tickrate>0);
            fs=double(~isempty(tz)*tz(1)+isempty(tz)*1000);
        end
    end

    % extract full pulse vector
    i0 = datastart(chanIdx); i1 = dataend(chanIdx);
    sig = double(data(i0:i1));

    % split into 2 segments: 0–15 s, 15–30 s
    segSamps = round(15*fs);
    above = sig(1:min(segSamps,end));
    below = sig(min(segSamps+1,end):min(2*segSamps,end));

    % compute amplitude (mV)
    ampAbove = mean_peak2peak_mV(above);
    ampBelow = mean_peak2peak_mV(below);

    % compute HR (bpm)
    hrAbove = heart_rate_bpm(above,fs);
    hrBelow = heart_rate_bpm(below,fs);

    rows(k).Name = names(k);
    rows(k).AmpAbove_mV = ampAbove;
    rows(k).AmpBelow_mV = ampBelow;
    rows(k).HRAbove_bpm = hrAbove;
    rows(k).HRBelow_bpm = hrBelow;
end

T4 = struct2table(rows);
disp(T4);
writetable(T4,'Table4.csv');
fprintf('Saved Table4.csv\n');

%% helper functions
function amp_mV = mean_peak2peak_mV(x)
    if isempty(x), amp_mV = NaN; return; end
    x = x(:) - median(x);
    [pks, locs] = findpeaks(x, 'MinPeakProminence', 0.05*std(x));
    if numel(pks)<2, amp_mV=NaN; return; end
    [trs,~] = findpeaks(-x, 'MinPeakProminence', 0.05*std(x));
    if isempty(trs), trs = -min(x); end
    meanpk = mean(pks,'omitnan'); meantr = mean(-trs,'omitnan');
    amp_mV = (meanpk - meantr)*1e3;  % V → mV
end

function hr = heart_rate_bpm(x,fs)
    if isempty(x), hr = NaN; return; end
    x = x(:) - median(x);
    [pks,locs] = findpeaks(x,'MinPeakProminence',0.05*std(x),'MinPeakDistance',round(0.3*fs));
    if numel(locs) < 2, hr = NaN; return; end
    ibi = diff(locs)/fs;             % inter-beat intervals (s)
    hr = 60/mean(ibi,'omitnan');     % beats per minute
end

%% Figure 8 and Table 5
clear; clc; close all;
files = { ...
   'defne_exercise4.mat', ...
   'jerry_exercise4.mat', ...
   'andrew_exercise4.mat', ...
   'izzy_exercise4.mat'};
subjects = {'Defne', 'Jerry', 'Andrew', 'Izzy'};
meanVals = zeros(length(subjects), 2);
stdVals = zeros(length(subjects), 2);
fprintf('===== DATA ANALYSIS =====\n\n');
for i = 1:length(files)
   fprintf('--- %s ---\n', subjects{i});
  
   data = load(files{i});
   varNames = fieldnames(data);
   vals = [];
  
   for v = 1:length(varNames)
       val = data.(varNames{v});
       if isnumeric(val) && size(val,2) >= 2
           vals = val;
           break;
       end
   end
  
   if isempty(vals)
       error('No valid 2-column numeric data found in file: %s', files{i})
   end
  
   hipData = vals(:,1);
   aboveData = vals(:,2);
  
   fprintf('Column 1 (hip): mean=%.6f, std=%.6f\n', mean(hipData), std(hipData));
   fprintf('Column 2 (above): mean=%.6f, std=%.6f\n', mean(aboveData), std(aboveData));
  
   % Check if identical
   if isequal(hipData, aboveData)
       fprintf('*** WARNING: BOTH COLUMNS ARE IDENTICAL! ***\n');
   elseif abs(mean(hipData) - mean(aboveData)) < 1e-10
       fprintf('*** WARNING: Means are essentially the same! ***\n');
   else
       fprintf('OK: Data differs between columns\n');
   end
  
   % Apply absolute value for Defne AND Andrew
   if strcmpi(subjects{i},'defne') || strcmpi(subjects{i},'andrew')
       hipData = abs(hipData);
       aboveData = abs(aboveData);
   end
  
   meanVals(i,1) = mean(aboveData);
   stdVals(i,1) = std(aboveData);
   meanVals(i,2) = mean(hipData);
   stdVals(i,2) = std(hipData);
  
   fprintf('After processing: Above=%.6f±%.6f, Hip=%.6f±%.6f\n\n', ...
       meanVals(i,1), stdVals(i,1), meanVals(i,2), stdVals(i,2));
end
fprintf('\n===== SUMMARY =====\n');
for i = 1:length(subjects)
   diffPercent = abs(meanVals(i,1) - meanVals(i,2)) / max(meanVals(i,:)) * 100;
   fprintf('%s: Difference = %.2f%%\n', subjects{i}, diffPercent);
end
% Transpose to group by activity
meanGrouped = meanVals';
stdGrouped = stdVals';
figure('Color', 'w');
b = bar(meanGrouped, 'grouped');
% Set specific colors for ALL 4 subjects
colors = [0 0.4470 0.7410;      % Blue for Defne
         0.8500 0.3250 0.0980;  % Orange for Jerry
         0.4660 0.6740 0.1880;  % Green for Andrew
         0.9290 0.6940 0.1250]; % Yellow for Izzy
for j = 1:length(b)
   b(j).FaceColor = colors(j,:);
   b(j).EdgeColor = 'k';
   b(j).LineWidth = 1;
end
hold on;
% Add error bars with enhanced visibility
numGroups = size(meanGrouped, 1);
numBars = size(meanGrouped, 2);
groupWidth = min(0.8, numBars / (numBars + 1.5));
for j = 1:numBars
   x = (1:numGroups) - groupWidth/2 + (2*j-1) * groupWidth / (2*numBars);
   y = meanGrouped(:,j);
   err = stdGrouped(:,j);
  
   % Make error bars thicker and with larger caps
   errorbar(x, y, err, 'k', 'linestyle', 'none', 'LineWidth', 2, 'CapSize', 10);
end
set(gca, 'XTick', 1:2, 'XTickLabel', ...
   {'Avg Pulse Amplitude Above Head', 'Avg Pulse Amplitude Below Hip'}, ...
   'FontSize', 11);
xlabel('Activity', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Voltage (V)', 'FontSize', 12, 'FontWeight', 'bold');
% Auto-scale y-axis
maxVal = max(meanGrouped(:) + stdGrouped(:));
ylim([0, maxVal * 1.2]);
legend(subjects, 'Location', 'northeast', 'FontSize', 10);
grid on;
hold off;
% Create a second figure zoomed in to show error bars
figure('Color', 'w');
b2 = bar(meanGrouped, 'grouped');
for j = 1:length(b2)
   b2(j).FaceColor = colors(j,:);
   b2(j).EdgeColor = 'k';
   b2(j).LineWidth = 1;
end
hold on;
for j = 1:numBars
   x = (1:numGroups) - groupWidth/2 + (2*j-1) * groupWidth / (2*numBars);
   y = meanGrouped(:,j);
   err = stdGrouped(:,j);
   errorbar(x, y, err, 'r', 'linestyle', 'none', 'LineWidth', 3, 'CapSize', 15);
end
set(gca, 'XTick', 1:2, 'XTickLabel', ...
   {'Above Head', 'Below Hip'}, 'FontSize', 11);
xlabel('Activity', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Voltage (V)', 'FontSize', 12, 'FontWeight', 'bold');
% Zoom in on the data range
minVal = min(meanGrouped(:) - stdGrouped(:));
maxVal = max(meanGrouped(:) + stdGrouped(:));
ylim([max(0, minVal - 0.0005), maxVal + 0.001]);
title('Zoomed View - Error Bars More Visible', 'FontSize', 14, 'FontWeight', 'bold');
legend(subjects, 'Location', 'best', 'FontSize', 10);
grid on;
hold off;
%% Table 1
clear; clc; close all;
% File names and subject names
files = { ...
   'defne_exercise2.mat', ...
   'jerry_exercise2.mat', ...
   'andrew_exercise2.mat', ...
   'izzy_exercise2.mat'};
subjects = {'Defne', 'Jerry', 'Andrew', 'Izzy'};
% Sampling frequency (adjust this based on your data!)
fs = 1000; % Hz - CHANGE THIS to match your actual sampling rate
% Storage for results
numSubjects = length(subjects);
meanRR = zeros(1, numSubjects);
maxRR = zeros(1, numSubjects);
minRR = zeros(1, numSubjects);
stdRR = zeros(1, numSubjects);
meanHR = zeros(1, numSubjects);
fprintf('===== R-R INTERVAL ANALYSIS =====\n\n');
% Process each subject
for i = 1:length(files)
   fprintf('--- Processing %s ---\n', subjects{i});
  
   % Load the data
   try
       data = load(files{i});
   catch
       fprintf('ERROR: Could not load file %s\n\n', files{i});
       continue;
   end
  
   % Extract ECG signal (adjust variable name as needed)
   varNames = fieldnames(data);
   ecg = [];
  
   for v = 1:length(varNames)
       val = data.(varNames{v});
       if isnumeric(val) && min(size(val)) == 1
           ecg = val(:); % Make it a column vector
           fprintf('Using variable: %s (length: %d samples)\n', varNames{v}, length(ecg));
           break;
       end
   end
  
   if isempty(ecg)
       fprintf('ERROR: No suitable ECG signal found\n\n');
       continue;
   end
  
   % Apply absolute value if needed (for subjects with negative voltages)
   if strcmpi(subjects{i},'defne') || strcmpi(subjects{i},'andrew')
       ecg = abs(ecg);
   end
  
   % Find R-peaks using findpeaks
   % You may need to adjust MinPeakHeight and MinPeakDistance
   [peaks, locs] = findpeaks(ecg, ...
       'MinPeakHeight', 0.6*max(ecg), ...  % Peaks must be at least 60% of max
       'MinPeakDistance', 0.5*fs);          % Minimum 0.5 seconds between peaks
  
   fprintf('Found %d R-peaks\n', length(locs));
  
   % Calculate R-R intervals (in samples, then convert to seconds)
   RR_samples = diff(locs);           % Difference between consecutive peak locations
   RR_seconds = RR_samples / fs;      % Convert to seconds
   RR_ms = RR_seconds * 1000;         % Convert to milliseconds
  
   % Calculate statistics
   meanRR(i) = mean(RR_ms);
   maxRR(i) = max(RR_ms);
   minRR(i) = min(RR_ms);
   stdRR(i) = std(RR_ms);
  
   % Convert R-R intervals to heart rate (bpm)
   % HR (bpm) = 60 / RR_interval (seconds)
   HR_values = 60 ./ RR_seconds;
   meanHR(i) = mean(HR_values);
  
   fprintf('R-R Interval Stats (ms):\n');
   fprintf('  Mean: %.2f ms\n', meanRR(i));
   fprintf('  Max:  %.2f ms\n', maxRR(i));
   fprintf('  Min:  %.2f ms\n', minRR(i));
   fprintf('  Std:  %.2f ms\n', stdRR(i));
   fprintf('Mean Heart Rate: %.2f bpm\n\n', meanHR(i));
end
% Create Table 1
fprintf('\n===== TABLE 1: R-R INTERVAL AND HEART RATE STATISTICS =====\n\n');
% Create the table
Table1 = table(meanRR', maxRR', minRR', stdRR', meanHR', ...
   'RowNames', subjects, ...
   'VariableNames', {'Mean_RR_ms', 'Max_RR_ms', 'Min_RR_ms', 'Std_RR_ms', 'Mean_HR_bpm'});
disp(Table1);
% Save the table if desired
% writetable(Table1, 'Table1_RR_Intervals.csv', 'WriteRowNames', true);
fprintf('\nTable created successfully!\n');
%%Table 2
clear; clc; close all;
% Subject names
subjects = {'Defne', 'Jerry', 'Andrew', 'Izzy'};
% Sampling frequency (adjust this based on your data!)
fs = 1000; % Hz - CHANGE THIS to match your actual sampling rate
% Storage for results
numSubjects = length(subjects);
meanRR = zeros(1, numSubjects);
maxRR = zeros(1, numSubjects);
minRR = zeros(1, numSubjects);
stdRR = zeros(1, numSubjects);
meanHR = zeros(1, numSubjects);
fprintf('========================================\n');
fprintf('TABLE 2: SITTING POSITION (Exercise 2)\n');
fprintf('========================================\n\n');
files = { ...
   'defne_exercise2.mat', ...
   'jerry_exercise2.mat', ...
   'andrew_exercise2.mat', ...
   'izzy_exercise2.mat'};
for i = 1:length(files)
   fprintf('--- Processing %s ---\n', subjects{i});
  
   % Load the data
   try
       data = load(files{i});
   catch
       fprintf('ERROR: Could not load file %s\n\n', files{i});
       continue;
   end
  
   % Extract ECG signal
   varNames = fieldnames(data);
   ecg = [];
  
   for v = 1:length(varNames)
       val = data.(varNames{v});
       if isnumeric(val) && min(size(val)) == 1
           ecg = val(:);
           fprintf('Using variable: %s (length: %d samples)\n', varNames{v}, length(ecg));
           break;
       end
   end
  
   if isempty(ecg)
       fprintf('ERROR: No suitable ECG signal found\n\n');
       continue;
   end
  
   % Apply absolute value if needed
   if strcmpi(subjects{i},'defne') || strcmpi(subjects{i},'andrew')
       ecg = abs(ecg);
   end
  
   % Find R-peaks
   [peaks, locs] = findpeaks(ecg, ...
       'MinPeakHeight', 0.6*max(ecg), ...
       'MinPeakDistance', 0.5*fs);
  
   fprintf('Found %d R-peaks\n', length(locs));
  
   % Calculate R-R intervals
   RR_samples = diff(locs);
   RR_seconds = RR_samples / fs;
   RR_ms = RR_seconds * 1000;
  
   % Calculate statistics
   meanRR(i) = mean(RR_ms);
   maxRR(i) = max(RR_ms);
   minRR(i) = min(RR_ms);
   stdRR(i) = std(RR_ms);
  
   % Convert to heart rate
   HR_values = 60 ./ RR_seconds;
   meanHR(i) = mean(HR_values);
  
   fprintf('R-R Interval Stats (ms): Mean=%.2f, Max=%.2f, Min=%.2f, Std=%.2f\n', ...
       meanRR(i), maxRR(i), minRR(i), stdRR(i));
   fprintf('Mean Heart Rate: %.2f bpm\n\n', meanHR(i));
end
% CREATE AND DISPLAY TABLE 2
fprintf('\n====================================================\n');
fprintf('           TABLE 2: FINAL RESULTS\n');
fprintf('====================================================\n\n');
Table2 = table(meanRR', maxRR', minRR', stdRR', meanHR', ...
   'RowNames', subjects, ...
   'VariableNames', {'Mean_RR_ms', 'Max_RR_ms', 'Min_RR_ms', 'Std_RR_ms', 'Mean_HR_bpm'});
disp(Table2);
fprintf('\nTable 2 created successfully!\n');


