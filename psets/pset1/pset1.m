% File Name: pset1_arodriguez.m
% Author: Andrew Rodriguez
% Created: September 13th 2025
% Decription: work for PSET 1 problems

% HH gating variables and conductances for a clamp: -65 mV -> +23 mV

%% 4a
% constants from the handout table
n0 = 0.3177; n_inf = 0.9494; tau_n = 1.2028; 
m0 = 0.0529; m_inf = 0.9953; tau_m = 0.1577;  
h0 = 0.5961; h_inf = 0.0009; tau_h = 1.0022;  

% time (ms)
t = linspace(0,8,201);

% Actual equations
n_t = n_inf - (n_inf - n0) .* exp(-t ./ tau_n);
m_t = m_inf - (m_inf - m0) .* exp(-t ./ tau_m);
h_t = h_inf - (h_inf - h0) .* exp(-t ./ tau_h);

%% 4b
% maximal conductances
gNa_bar = 120; 
gK_bar  = 36;

% conductances vs time
gNa_t = gNa_bar .* (m_t.^3) .* h_t;
gK_t  = gK_bar  .* (n_t.^4);


%% 4c
% plot
figure; hold on;
plot(t, gK_t,  'r', 'LineWidth', 1.8);   % gK red
plot(t, gNa_t, 'y', 'LineWidth', 1.8);   % gNa yellow
xlabel('time (ms)');
ylabel('conductance (mS/cm^2)');
legend('g_K(t)','g_{Na}(t)','Location','northeast');
title('Voltage clamp: -65 \rightarrow +23 mV');
grid on; box on;

% save
exportgraphics(gcf,'gK_gNa_0to8ms.png','Resolution',300);

%% 4e
E_Na = 66; % Assuming we use numbers from class
E_K  = -94; 

% Chord conductance
Vmem_t = (gK_t*E_K + gNa_t*E_Na) ./ (gK_t + gNa_t);

% plot
figure; plot(t, Vmem_t, 'LineWidth', 1.8);
xlabel('time (ms)'); ylabel('V_{mem} (mV)');
title('V_{mem}(t) from chord conductance using g_K(t) and g_{Na}(t)');
grid on; box on;

% save
exportgraphics(gcf,'Vmem_from_chord.png','Resolution',300);



%% 5a

% time
dt = 0.001;
t  = -30:dt:30;
% sin c
y = sinc(t);

figure('Name','Figure 5a'); 
plot(t, y, 'k', 'LineWidth', 1.5);
xlabel('t'); ylabel('sinc(t)'); title('Figure 5a: sinc(t)'); grid on; box on;
exportgraphics(gcf,'5a.png','Resolution',300);

%% 5b

% Derivative
dy = gradient(y, dt); 

figure('Name','Figure 5b');
plot(t, dy, 'k', 'LineWidth', 1.5);
xlabel('t'); ylabel('d/dt sinc(t)'); title('Figure 5b: derivative of sinc'); grid on; box on;
exportgraphics(gcf,'5b.png','Resolution',300);


%% 5c

% Zero-crossings of derivative
sgn = sign(dy);
ix  = find(diff(sgn) ~= 0);
t_zc = t(ix) - dy(ix) .* (t(ix+1)-t(ix)) ./ (dy(ix+1)-dy(ix)); 

figure('Name','Figure 5c'); hold on;
plot(t(dy>0),  dy(dy>0),  'g.', 'MarkerSize', 5);   % positive portions
plot(t(dy<0),  dy(dy<0),  'r.', 'MarkerSize', 5);   % negative portions
scatter(t_zc, zeros(size(t_zc)), 36, 'ys', 'filled', 'MarkerEdgeColor','k');
xlabel('t'); ylabel('d/dt sinc(t)');
legend({'deriv > 0','deriv < 0','zero-crossings'}, 'Location','northeast');
title('Figure 5c: derivative sign + zero-crossings'); grid on; box on;

exportgraphics(gcf,'5c.png','Resolution',300);


%% Figure 5d: original sinc colored by derivative sign; peaks/troughs marked

y_zc = interp1(t, y, t_zc, 'linear');

figure('Name','Figure 5d'); hold on;
plot(t(dy>0), y(dy>0), 'g.', 'MarkerSize', 5);     % where derivative is positive
plot(t(dy<0), y(dy<0), 'r.', 'MarkerSize', 5);     % where derivative is negative
scatter(t_zc, y_zc, 36, 'ys', 'filled', 'MarkerEdgeColor','k'); % peaks/troughs of sinc
xlabel('t'); ylabel('sinc(t)');
legend({'deriv > 0','deriv < 0','peaks & troughs'}, 'Location','northeast');
title('Figure 5d: sinc colored by derivative sign, extrema marked');
grid on; box on;
exportgraphics(gcf,'5d.png','Resolution',300);