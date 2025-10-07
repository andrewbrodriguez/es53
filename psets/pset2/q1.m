T_final = 10; % so we run for 8 ms
pulse_I = 100; % use a pulse heigh of 100
pulse_w_ms = 0.1; % and then the pulse width is 0.1 ms

[t,y] = run_hh_model(T_final, pulse_I, pulse_w_ms); % run model
Vm = y(:,1); m = y(:,2); h = y(:,3); n = y(:,4);
% first element is membrane potential, second is m, third is h, fourth is n
plot(t, Vm, 'k','LineWidth',1.5); grid on;

%%
% Conductances from gates 
gNa = 120*(m.^3).*h;
gK  = 36*(n.^4);

figure('Color','w');

% So then we have the the membrane potential
subplot(3,1,1);
plot(t, Vm, 'k','LineWidth',1.5); grid on;
xlabel('Time (ms)'); ylabel('V_m (mV)'); title('(i) Membrane potential');

% Then the conductance
subplot(3,1,2);
plot(t, gNa, 'b', 'LineWidth',1.2); hold on;
plot(t, gK , 'r', 'LineWidth',1.2); grid on;
xlabel('Time (ms)'); ylabel('Conductance (mS/cm^2)');
title('(ii) g_{Na} and g_{K}'); legend({'g_{Na}','g_{K}'},'Location','northeast');

% Then the open probabilities
subplot(3,1,3);
plot(t, m, 'g', 'LineWidth',1.2); hold on;
plot(t, h, 'm', 'LineWidth',1.2);
plot(t, n, 'Color',[0.3 0.3 0.3], 'LineWidth',1.2);
grid on;
xlabel('Time (ms)'); ylabel('Gate value');
title('(iii) Gates m, h, n');
legend({'m','h','n'},'Location','southeast');


saveas(gcf, 'q1a.png');

%% q1b, now we are doing 100 loops

clear; close all; clc; % clean slate

T_final = 8; % still running 8 ms
pulse_w_ms = 0.1; % still 0.1 ms pulse
amps = 0:1:100;  % now we are ranging from 0 to 100 as our intensity

peakVm = nan(size(amps));
Vm_traces = cell(size(amps)); % saves these for later usage

% now we wanna plot all 100 on the same graph
f1 = figure('Name','Q1(b) Vm sweeps','Color','w'); hold on; box on; grid on;
for k = 1:numel(amps)
    [t, y] = run_hh_model(T_final, amps(k), pulse_w_ms); % run a model
    v = y(:,1); % pull out just our membrane potential
    Vm_traces{k} = v; % append to our list for later
    peakVm(k) = max(v); % check max

    % and plot
    plot(t, v, 'Color', [1.0 0.0 0.0], 'LineWidth', 0.5);
end

xlabel('Time (ms)'); ylabel('V_m (mV)');
title('V_m(t) for 0:100 \muA/cm^2');

saveas(gcf, 'q1b.png');
%% q1c, now we are figuring out the relationship between input and peak

f2 = figure('Name','q1c Peak V vs I','Color','w'); box on; grid on; hold on;
plot(amps, peakVm, 'o-', 'LineWidth', 1.2, 'MarkerSize', 4);
yline(0,'k:','0 mV');

% We can use find to find the first current amount to cause an AP
thrIdx = find(peakVm >= 0, 1, 'first'); 
thrI = amps(thrIdx);
xline(thrI, 'r--', sprintf('~%d \\muA/cm^2', thrI));
saveas(gcf, 'q1c.png');

%% q1d, now lets compare just above threshold to a strong signal

I_thr = 18; % 18 is just above the threshold 
I_str = 100; % strong pulse

% Now lets run the two cases
[t_thr, y_thr] = run_hh_model(T_final, I_thr, pulse_w_ms);
[t_str, y_str] = run_hh_model(T_final, I_str, pulse_w_ms);
% and pull out the vmem curves
vm_thr = y_thr(:,1);
vm_str = y_str(:,1);

figure('Name','q1d threshold vs strong','Color','w'); hold on; box on; grid on;
plot(t_thr, vm_thr, 'b', 'LineWidth', 1.6);
plot(t_str, vm_str, 'r', 'LineWidth', 1.6);
xlabel('Time (ms)'); ylabel('V_m (mV)');
legend({sprintf('18 uA/cm2 (just above threshold)', I_thr), '100 uA/cm2 (strong)'}, 'Location','northeast');
saveas(gcf, 'q1d.png');