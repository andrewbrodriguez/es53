%% ===== Biceps EMG RMS vs Weight (incl. 0-lb baseline; x-axis = lbs) =====
clear; clc;
% ------------ SETTINGS ------------
filePattern = 'andrew_data_exercise1.mat';   % adjust if your filenames differ
weight_lbs  = [0 2.2 4.4 6.6 8.8];       % <-- put your actual weights (lbs) here
weights     = 0:numel(weight_lbs)-1;     % weight labels expected in notes (0,1,2,3,4)
trimStart   = 1.0;                       % s after segment start to skip (settling)
trimEnd     = 0.5;                       % s before segment end to skip (movement)
scale_to_mV = 1000;                      % set to 1 if your data are already in mV
% ------------ FILES ------------
files = dir(filePattern);
if isempty(files), error('No files matching "%s".', filePattern); end
% Safe indexing helpers
getRC = @(M,r,c) M(r, min(c, size(M,2)));   % matrix [rows x cols]
getC  = @(M,c)   M(1, min(c, size(M,2)));   % row vector (e.g., tickrate)
RMS_by_subj = nan(numel(files), numel(weights));   % rows=subj, cols=weights
subjNames   = strings(numel(files),1);
for s = 1:numel(files)
   % ---- Load LabChart struct ----
   S = load(files(s).name);
   fn = fieldnames(S);
   if numel(fn)==1 && isstruct(S.(fn{1}))
       L = S.(fn{1});
   else
       L = S;
   end
   subjNames(s) = erase(files(s).name, '.mat');
   % ---- Find channels by name ----
   chTitles = strtrim(cellstr(L.titles));    % e.g., {'RMS Biceps','RMS Triceps','Biceps','Triceps'}
   chn_rms_biceps = find(strcmpi(chTitles,'RMS Biceps'),1); if isempty(chn_rms_biceps), chn_rms_biceps = 1; end
   chn_biceps_raw = find(strcmpi(chTitles,'Biceps'),1);     if isempty(chn_biceps_raw), chn_biceps_raw = 3; end
   % ---- Auto-pick longest block on raw biceps channel ----
   nBlocks = size(L.datastart,2);
   durations = zeros(1,nBlocks);
   for b = 1:nBlocks
       i0b = getRC(L.datastart , chn_biceps_raw, b);
       i1b = getRC(L.dataend   , chn_biceps_raw, b);
       fsb = getRC(L.samplerate, chn_biceps_raw, b);
       if i0b>0 && i1b>=i0b && fsb>0
           durations(b) = (i1b - i0b + 1)/fsb;
       end
   end
   [~, bn] = max(durations);
   if durations(bn) <= 0
       warning('No valid data block for %s. Skipping.', files(s).name);
       continue
   end
   % ---- Raw biceps data from chosen block ----
   fs_raw = getRC(L.samplerate, chn_biceps_raw, bn);
   i0_raw = getRC(L.datastart , chn_biceps_raw, bn);
   i1_raw = getRC(L.dataend   , chn_biceps_raw, bn);
   x_raw  = double(L.data(i0_raw:i1_raw));
   tEnd   = (i1_raw - i0_raw + 1)/fs_raw;
   % ---- Notes: global ticks -> seconds relative to this block start ----
   tick     = getC(L.tickrate, bn);                                 % ticks/sec
   t0_block = getRC(L.firstsampleoffset, chn_biceps_raw, bn)/tick;  % block start (s, absolute)
   note_tick = L.com(:,3);                                          % global tick index
   note_txt  = strtrim(cellstr(L.comtext));
   t_note    = note_tick/tick - t0_block;                           % seconds rel. to block start
   % Keep notes inside this block's duration
   keep     = (t_note >= 0) & (t_note <= tEnd);
   t_note   = t_note(keep);
   note_txt = note_txt(keep);
   [t_note, ord] = sort(t_note);
   note_txt = note_txt(ord);
   % ---- Parse weight number from note text; no number => treat as 0 (rest)
   W_notes = zeros(size(t_note));
   for k = 1:numel(t_note)
       m = regexp(lower(note_txt{k}), '\d+', 'match');
       if ~isempty(m), W_notes(k) = str2double(m{1}); end
   end
   % ---- Build segments: include baseline from 0 s (weight 0) ----
   startTimes = [0; t_note(:)];
   labels     = [0; W_notes(:)];
   endTimes   = [t_note(:); tEnd];
   % ---- Choose RMS source ----
   use_rms_channel = getRC(L.datastart, chn_rms_biceps, bn) > 0;
   if use_rms_channel
       fs_rms = getRC(L.samplerate, chn_rms_biceps, bn);
       j0     = getRC(L.datastart , chn_rms_biceps, bn);
       j1     = getRC(L.dataend   , chn_rms_biceps, bn);
       x_rms  = double(L.data(j0:j1));       % precomputed RMS trace (likely V)
   else
       fs_rms = fs_raw;
       x_rms  = x_raw;                        % we'll compute RMS from raw windows
   end
   % ---- Average RMS per weight (stable window in each segment) ----
   for k = 1:numel(labels)
       wlbl = labels(k);
       if ~ismember(wlbl, weights), continue; end
       t1 = startTimes(k) + trimStart;
       t2 = endTimes(k)   - trimEnd;
       if t2 <= t1, continue; end
       if use_rms_channel
           a = max(1, round(t1*fs_rms) + 1);
           b = min(numel(x_rms), round(t2*fs_rms));
           seg = x_rms(a:b);
           avgrms = mean(seg);
       else
           a = max(1, round(t1*fs_raw) + 1);
           b = min(numel(x_raw), round(t2*fs_raw));
           seg = x_raw(a:b);
           seg = seg - mean(seg);
           avgrms = sqrt(mean(seg.^2));
       end
       col = find(weights == wlbl, 1, 'first');       % 1..numel(weights)
       if ~isempty(col)
           RMS_by_subj(s, col) = avgrms * scale_to_mV;   % store in mV
       end
   end
end
%% ---------- CLEAN SUBJECT NAMES ----------
names = regexprep(string(subjNames), '\.mat$','');
names = regexprep(names, 'Data\s*Exercise\s*1', '', 'ignorecase');  % remove "Data Exercise 1"
names = regexprep(names, '_', ' ');
%% ---------- TABLE 1 (rows=weights incl. rest; cols=subjects + Group Mean/SD) ----------
% Pretty row labels
words = ["One","Two","Three","Four","Five","Six","Seven","Eight","Nine","Ten"];
rowLabels = strings(numel(weight_lbs),1);
rowLabels(1) = "Rest (" + string(weight_lbs(1)) + " lbs)";
for i = 2:numel(weight_lbs)
   plural = ""; if i >= 3, plural = "s"; end   % "weights" for 3,4,...
   rowLabels(i) = words(i-1) + " weight" + plural + " (" + string(weight_lbs(i)) + " lbs)";
end
% Group stats across subjects (handles NaNs)
groupMean = mean(RMS_by_subj, 1, 'omitnan').';   % (weights x 1)
groupSD   = std( RMS_by_subj, 0, 1, 'omitnan').';
% Build table in desired order
Table1 = table(rowLabels, weight_lbs(:), 'VariableNames', {'Condition','Weight_lbs'});
for s = 1:size(RMS_by_subj,1)
   colName = matlab.lang.makeValidName(string(names(s)));
   Table1.(colName) = round(RMS_by_subj(s,:).', 3);
end
Table1.GroupMean_mV = round(groupMean, 3);
Table1.GroupSD_mV   = round(groupSD,   3);
disp('Table 1 (RMS Biceps EMG, mV):');
disp(Table1);
% Save CSV for your report
writetable(Table1, 'Table1_RMS_Biceps_by_Weight_AllSubjects.csv');
% UI table (nice for screenshots)
fig = uifigure('Name','Table 1: RMS Biceps EMG (mV)');
uit = uitable(fig, 'Data', Table1);
uit.Position = [20 20 1000 280];
%% ---------- GROUP BAR (mean ± SD; x = lbs) ----------
figure;
bar(weight_lbs, groupMean); hold on;
errorbar(weight_lbs, groupMean, groupSD, 'k.', 'LineWidth', 1.4);
hold off; grid on;
xlim([min(weight_lbs) max(weight_lbs)]);
xticks(weight_lbs);
xlabel('Weight (lbs)');
ylabel('Biceps EMG Average RMS (mV)');
title('Biceps EMG RMS vs Weight (Group mean \pm SD)');
%% ---------- PER-SUBJECT LINES (incl. 0 lb; x = lbs) ----------
figure; hold on
for s = 1:size(RMS_by_subj,1)
   plot(weight_lbs, RMS_by_subj(s,:), '--*', 'LineWidth', 1.6, 'MarkerSize', 8);
end
hold off; grid on
xlim([min(weight_lbs) max(weight_lbs)]);
xticks(weight_lbs);
xlabel('Weight (lbs)');
ylabel('RMS Biceps Amplitude (mV)');
legend(names, 'Location', 'northwest');
title('RMS Biceps EMG vs Weight (per subject, incl. 0 lb)');

