S = load('ECGdata.mat'); 
t = S.ptime(:);
p = S.pdata(:);
p10 = smooth(p, 10,   'moving'); 
p100 = smooth(p, 100,  'moving'); 
p1000 = smooth(p, 1000, 'moving'); 

figure('Name','ECG smoothing (original vs spans of 10, 100, 1000 pts)');
tiledlayout(4,1,'Padding','compact','TileSpacing','compact');

nexttile;  plot(t, p, 'LineWidth', 1); % Original
title('Original ECG'); ylabel('mV'); grid on;

nexttile;  plot(t, p10, 'LineWidth', 1); % 10 point
title('Smoothed (10 points)'); ylabel('mV'); grid on;

nexttile;  plot(t, p100, 'LineWidth', 1); % 100 point
title('Smoothed (100 points)'); ylabel('mV'); grid on;

nexttile;  plot(t, p1000, 'LineWidth', 1); % 1000 point
title('Smoothed (1000 points)'); xlabel('Time (s)'); ylabel('mV'); grid on;
saveas(gcf, "q6a.png")
% Smoothing lets us average areas to reduce high frequency noise.
% But for smoothing to really work we need a lot of data, and it needs to
% be high frequnecy. What smoothing does is it takes a set of points
% and it averages them to give us a single point for that set of points. 
% So for 10 points we get light smoothing but the QRS peaks are preserved
% for 100 points we stronger smoothing with noticeable rounding of QRS
% then at 1000 points the shape is largely lost and theres a lag 

%% 6b
fs = 1/median(diff(t)); % sampling rate
prom = 0.25*std(p10); % prominence threshold (how high or low)
mindist = 0.50*fs; % peaks are 0.5 fs apart

% Peaks (max) and troughs (min)
[pks, ipks] = findpeaks(p10,  'MinPeakProminence', prom, 'MinPeakDistance', mindist);
[trs, itrs] = findpeaks(-p10, 'MinPeakProminence', prom, 'MinPeakDistance', mindist);
trs = -trs;

figure('Name','ECG peaks & troughs (span=10)'); hold on; grid on; box on;
plot(t, p10, 'LineWidth', 1);
scatter(t(ipks), pks, 36, 'ys', 'filled', 'MarkerEdgeColor','k');  % peaks
scatter(t(itrs), trs, 36, 'ys', 'filled', 'MarkerEdgeColor','k');  % troughs
xlabel('Time (s)'); ylabel('mV'); title('Smoothed ECG (smoothed span of 10 pts): peaks & troughs');
legend('Smoothed ECG','Peaks','Troughs','Location','best');
saveas(gcf, "q6b.png")