
%% Lab 6
% Title: Lab6
% Author: Elena Luciano
% Date: 2025-12-11
% Description: Lab5
%% Figure 3
function fig3_rebreathing_pulseHR_180
files.Andrew   = 'Exercise3_Andrew_lab6.mat';
files.Elena    = 'Exercise3_Elena_lab6.mat';
files.Evelyn   = 'Exercise3_Evelyn_lab6.mat';   % Subject C = Evelyn
SPO2.Andrew = [98 98 98 97 97 96 97 98 98 97 97 98 97 98 97 98 98 97];
SPO2.Elena  = [99 99 99 99 99 98 98 98 98 98 98 98 98 97 85 85 82 82];
SPO2.Evelyn = [98 98 98 97 98 98 97 97 98 98 98 98 98 98 98 98 98 98];
rc.Andrew = [1 3; 2 3];   % Breath row1 col3, Pulse row2 col3
rc.Elena  = [1 3; 2 3];   % same
% Evelyn handled explicitly
subjects = {'Andrew','Elena','Evelyn'};
titles   = {'Andrew','Elena','Evelyn'};
maxT = 180;   % seconds
figure('Position',[80 80 1120 840]);
for k = 1:3
   nm   = subjects{k};
   spo2v= SPO2.(nm);
   subplot(3,1,k);
  
   if strcmp(nm,'Evelyn')
      
       S  = load(files.Evelyn);
       fs = double(S.samplerate(1));
       raw = double(S.data(:));
       ds = double(S.datastart(:));   % 3x1
       de = double(S.dataend(:));     % 3x1
       breathC = raw(ds(1):de(1));    % row 1 = Breath
       pulseC  = raw(ds(2):de(2));    % row 2 = Pulse
       N = min(numel(breathC), numel(pulseC));
       breath = breathC(1:N);
       pulse  = pulseC(1:N);
       t      = (0:N-1)'/fs;
   else
      
       S  = load(files.(nm));
       fs = double(S.samplerate(1));
       raw = double(S.data(:));
       ds = double(S.datastart);
       de = double(S.dataend);
       if isvector(ds), ds = ds(:); end
       if isvector(de), de = de(:); end
       idx = rc.(nm);
       breath = slice_seg(raw, ds, de, idx(1,1), idx(1,2));  % Breath
       pulse  = slice_seg(raw, ds, de, idx(2,1), idx(2,2));  % Pulse
       N = min(numel(breath), numel(pulse));
       breath = breath(1:N);
       pulse  = pulse(1:N);
       t      = (0:N-1)'/fs;
   end
  
   keep = t <= maxT;
   t      = t(keep);
   breath = breath(keep);
   pulse  = pulse(keep);
  
   HR = hr_from_pulse(pulse, fs, t);   % bpm, same length as t
   % A) Breathing (left axis)
   yyaxis left
   plot(t, breath, 'b', 'LineWidth', 1.1); hold on
   ylabel('voltage (mV)'); box on
   % B) Heart Rate + SpO2 (right axis)
   yyaxis right
   plot(t, HR, 'r', 'LineWidth', 1.1); hold on
   if ~isempty(spo2v)
       t_sp = (0:numel(spo2v)-1)' * 10;  % every 10 s
       keep_sp = t_sp <= maxT;           % cut SPO2 at 180 s too
       t_sp = t_sp(keep_sp);
       sp   = spo2v(keep_sp);
       plot(t_sp, sp, 'm', 'LineWidth', 1.6);
   end
   % Nice y-limits for HR + SpO2
   allHR  = HR(isfinite(HR));
   allSpO = spo2v(:);
   yr_low  = min([60; min(allHR,[],'omitnan'); min(allSpO,[],'omitnan')]);
   yr_high = max([140; max(allHR,[],'omitnan'); max(allSpO,[],'omitnan')]);
   ylim([yr_low yr_high]);
   ylabel('heart rate (bpm) / SpO2 (%)');
   % Legend & titles
   hb  = plot(nan,nan,'b','LineWidth',1.1);
   hhr = plot(nan,nan,'r','LineWidth',1.1);
   hsp = plot(nan,nan,'m','LineWidth',1.6);
   legend([hb hhr hsp], {'Breathing','Heart Rate (from pulse)','SpO2'}, 'Location','northeast');
   title(titles{k});
   if k == 3
       xlabel('Time (s)');
   end
end
sgtitle('Exercise 3 (Rebreathing, first 180 s): Breathing (mV), Heart Rate from Pulse (bpm), SpO_2 (%)');
end
function x = slice_seg(raw, ds, de, r, c)
if isvector(ds), ds = ds(:); end
if isvector(de), de = de(:); end
if r<=size(ds,1) && c<=size(ds,2) && r<=size(de,1) && c<=size(de,2)
   s = ds(r,c); e = de(r,c);
   if s>0 && e>0 && e>=s && e<=numel(raw)
       x = raw(s:e);
       return
   end
end
x = [];
end
function HR = hr_from_pulse(pulse, fs, t)
sig = double(pulse(:)) - median(double(pulse(:)));
[~, locs] = findpeaks(sig, ...
   'MinPeakDistance', round(0.35*fs), ...        
   'MinPeakProminence', max(0.1*std(sig), 0.02*range(sig)));
if numel(locs) < 2
   HR = nan(size(t));
   return;
end
% Inter-beat intervals
ibi = diff(locs) / fs;      % seconds
hr_inst = 60 ./ ibi;        % bpm
% Put HR at midpoints between beats
tmid = (locs(1:end-1) + locs(2:end)) / 2 / fs;
% Interpolate to full time vector and smooth
HR = interp1(tmid, hr_inst, t, 'linear', 'extrap');
HR = clean_hr(HR);
end
function HR = clean_hr(hr)
HR = double(hr(:));
HR(~isfinite(HR) | HR<=30 | HR>=220) = NaN;    
HR = fillmissing(HR,'linear','EndValues','nearest');
HR = movmedian(HR,5,'omitnan');
HR = movmean(HR,5,'omitnan');
end
%% Table 4
files.Andrew = 'Exercise3_Andrew_lab6.mat';
files.Elena  = 'Exercise3_Elena_lab6.mat';
files.Evelyn   = 'Exercise3_Evelyn_lab6.mat';
tStart = 60;
tEnd   = 120;
% ANDREW
S = load(files.Andrew);
fsA  = double(S.samplerate(1));
rawA = double(S.data(:));
dsA  = double(S.datastart);
deA  = double(S.dataend);
% Breath: row 1 col 3
breathA = rawA(dsA(1,3):deA(1,3));
tA = (0:length(breathA)-1)'/fsA;
% Window 60–120 s
idxA = (tA >= tStart) & (tA <= tEnd);
[pksA, locsA] = findpeaks(breathA(idxA), tA(idxA));
fprintf('Andrew: %d breaths in 60–120 s\n', numel(locsA));
% ELENA
S = load(files.Elena);
fsE  = double(S.samplerate(1));
rawE = double(S.data(:));
dsE  = double(S.datastart);
deE  = double(S.dataend);
% Breath: row 1 col 3
breathE = rawE(dsE(1,3):deE(1,3));
tE = (0:length(breathE)-1)'/fsE;
% Window 60–120 s
idxE = (tE >= tStart) & (tE <= tEnd);
[pksE, locsE] = findpeaks(breathE(idxE), tE(idxE));
fprintf('Elena:  %d breaths in 60–120 s\n', numel(locsE));
% EVELYN
S = load(files.Evelyn);
fsC  = double(S.samplerate(1));
rawC = double(S.data(:));
dsC  = double(S.datastart);
deC  = double(S.dataend);
% Breath: row 3 col 1
breathC = rawC(dsC(3,1):deC(3,1));
tC = (0:length(breathC)-1)'/fsC;
% Window 60–120 s
idxC = (tC >= tStart) & (tC <= tEnd);
[pksC, locsC] = findpeaks(breathC(idxC), tC(idxC));
fprintf('Evelyn: %d breaths in 60–120 s\n', numel(locsC));
