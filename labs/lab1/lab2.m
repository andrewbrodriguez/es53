% File Name: lab2_arodriguez.m
% Author: Andrew Rodriguez
% Created: September 23rd 2025
% Description: Lab 2, Coactivation — Table 1 and Figure 2 (Nina uses block 2)

clear; clc;

scale_to_mV = 1000;                                   % set to 1 if already mV
subjects = ["andrew","daniel","evelyn","nina"];
files    = subjects + "_data_exercise2.mat";

% Snipped data time points for each activation 
W = struct();
W.andrew.biceps_act   = [20 22; 32 34; 42 44; 50 52; 59 60];
W.andrew.triceps_act  = [25 27; 36 39; 46 48; 55 57; 62 63];

W.daniel.biceps_act   = [20 22; 31 33; 40 42; 48 50; 57 58];
W.daniel.triceps_act  = [26 28; 36 38; 44 46; 52 54; 60 62];

W.evelyn.biceps_act   = [29 30; 37 39; 42 45; 49 52; 57 59];
W.evelyn.triceps_act  = [34 36; 40 42; 46 49; 54 56; 61 63];

W.nina.biceps_act     = [ 5  7; 13 16; 22 24; 29 31; 37 39];
W.nina.triceps_act    = [ 9 11; 18 20; 26 28; 33 35; 42 44];

rowNames = ["Biceps Contracting"; "Triceps Contracting"];
nSubj = numel(files);
RMS_Biceps  = nan(2, nSubj);
RMS_Triceps = nan(2, nSubj);
ratio_bicepsContr  = nan(nSubj,1); % (biceps/triceps) during biceps-contracting
ratio_tricepsContr = nan(nSubj,1); % (triceps/biceps) during triceps-contracting

for s = 1:nSubj
    subj = subjects(s);
    fp   = files(s);
    fprintf('>> %s\n', fp);

    S = load(fp); fn = fieldnames(S);
    if numel(fn)==1 && isstruct(S.(fn{1})), L = S.(fn{1}); else, L = S; end

    chRMS_B = 1;
    chRMS_T = 2;

    % Nina-only: force block 2; others: only one block
    if strcmpi(subj, "nina")
        block_override = 2;
    else
        block_override = [];
    end

    win_bi = W.(char(subj)).biceps_act;
    win_tr = W.(char(subj)).triceps_act;

    bi_B_vals = mean_over_windows(L, chRMS_B, win_bi, scale_to_mV, block_override);
    bi_T_vals = mean_over_windows(L, chRMS_T, win_bi, scale_to_mV, block_override);

    if ~isempty(bi_B_vals) && ~isempty(bi_T_vals)
        RMS_Biceps(1, s)  = mean(bi_B_vals, 'omitnan');
        RMS_Triceps(1, s) = mean(bi_T_vals, 'omitnan');
        if RMS_Triceps(1,s) > 0
            ratio_bicepsContr(s) = RMS_Biceps(1,s) / RMS_Triceps(1,s);
        end
    else
        warning('No usable biceps-contracting windows in %s', fp);
    end

    % ----- Triceps contracting -----
    tr_T_vals = mean_over_windows(L, chRMS_T, win_tr, scale_to_mV, block_override);
    tr_B_vals = mean_over_windows(L, chRMS_B, win_tr, scale_to_mV, block_override);

    if ~isempty(tr_T_vals) && ~isempty(tr_B_vals)
        RMS_Biceps(2, s)  = mean(tr_B_vals, 'omitnan');
        RMS_Triceps(2, s) = mean(tr_T_vals, 'omitnan');
        if RMS_Biceps(2,s) > 0
            ratio_tricepsContr(s) = RMS_Triceps(2,s) / RMS_Biceps(2,s);
        end
    else
        warning('No usable triceps-contracting windows in %s', fp);
    end
end

% Table 1
T = table(rowNames, 'VariableNames', {'Condition'});
for s = 1:nSubj
    name_str = char(subjects(s)); name_str(1) = upper(name_str(1));
    base = matlab.lang.makeValidName(name_str);
    T.(base + "_RMS_Biceps_mV")  = round(RMS_Biceps(:,s), 3);
    T.(base + "_RMS_Triceps_mV") = round(RMS_Triceps(:,s), 3);
end
disp('----- Table 1. Coactivation RMS (mV) -----'); disp(T);
writetable(T, 'Table1_Coactivation_RMS_Biceps_Triceps.csv');

% UI
try
    fT = uifigure('Name','Table 1: Coactivation RMS (mV)');
    uitable(fT, 'Data', T, 'Position', [20 20 1000 280]);
catch, end

% Figure 2
grpMeans = [mean(ratio_bicepsContr,  'omitnan'), mean(ratio_tricepsContr, 'omitnan')];
grpSDs   = [std( ratio_bicepsContr, 0,'omitnan'), std( ratio_tricepsContr,0,'omitnan')];

x = categorical({'Biceps contracting','Triceps contracting'});
x = reordercats(x, cellstr(x));

figure;
bar(x, grpMeans); hold on;
errorbar(1:2, grpMeans, grpSDs, 'k.', 'LineWidth', 1.25);
ylabel('RMS ratio (contracting / non-contracting)'); grid on;
saveas(gcf, 'Figure2_Coactivation_RatioBar.png');

fprintf('\nGroup ratio (Biceps contracting):  mean = %.3f, SD = %.3f, n = %d\n', ...
    grpMeans(1), grpSDs(1), sum(~isnan(ratio_bicepsContr)));
fprintf('Group ratio (Triceps contracting): mean = %.3f, SD = %.3f, n = %d\n', ...
    grpMeans(2), grpSDs(2), sum(~isnan(ratio_tricepsContr)));

% =================== Helpers ===================
function vals = mean_over_windows(L, chIdx, abs_windows, scale_to_mV, block_override)
    % If block_override provided, use exactly that block; otherwise search all blocks.
    if isempty(abs_windows), vals = []; return; end
    if nargin < 5, block_override = []; end
    if size(abs_windows,2) ~= 2, abs_windows = reshape(abs_windows, [], 2); end
    nW = size(abs_windows,1);
    vals = nan(nW,1);

    if isempty(block_override)
        % Search all blocks; each window assigned to the first block that contains it
        nBlocks = size(L.datastart, 2);
        for b = 1:nBlocks
            [ok, fs, t0, x] = pull_block(L, chIdx, b, scale_to_mV);
            if ~ok, continue; end
            len = numel(x);
            for k = 1:nW
                if ~isnan(vals(k)), continue; end
                t1 = abs_windows(k,1) - t0;  t2 = abs_windows(k,2) - t0;
                a = max(1, floor(t1*fs)+1);  b2 = min(len, floor(t2*fs));
                if b2 > a, vals(k) = mean(x(a:b2)); end
            end
        end
    else
        % Use only the specified block
        [ok, fs, t0, x] = pull_block(L, chIdx, block_override, scale_to_mV);
        if ok
            len = numel(x);
            for k = 1:nW
                t1 = abs_windows(k,1) - t0;  t2 = abs_windows(k,2) - t0;
                a = max(1, floor(t1*fs)+1);  b2 = min(len, floor(t2*fs));
                if b2 > a, vals(k) = mean(x(a:b2)); end
            end
        end
    end
    vals = vals(~isnan(vals)); % drop windows that didn't map to data
end

function [ok, fs, t0, x] = pull_block(L, chIdx, b, scale_to_mV)
    ok = false; fs = []; t0 = []; x = [];
    if chIdx > size(L.datastart,1) || b > size(L.datastart,2), return; end
    j0 = L.datastart(chIdx, b); j1 = L.dataend(chIdx, b);
    if j0<=0 || j1<=j0, return; end
    if chIdx > size(L.samplerate,1) || b > size(L.samplerate,2), return; end
    fs = L.samplerate(chIdx, b); if fs<=0, return; end

    % tickrate handling
    if isvector(L.tickrate), tick = L.tickrate(min(numel(L.tickrate), b));
    else, tick = L.tickrate(1, min(size(L.tickrate,2), b));
    end
    if tick<=0, tick = 1; end

    if size(L.firstsampleoffset,1) >= chIdx && size(L.firstsampleoffset,2) >= b
        t0 = L.firstsampleoffset(chIdx, b)/tick;
    else
        t0 = 0;
    end

    x = double(L.data(j0:j1)) * scale_to_mV;  % RMS channel (mV)
    ok = true;
end