% PV loop plot
S = load('PVloop.mat');       
M = S.PVloop;
V = M(1,:); % volumes in mL
P = M(2,:); % pressures in mmHg

figure; hold on; box on; grid on;
plot(V, P, 'k-', 'LineWidth', 2);   
plot([V(end) V(1)], [P(end) P(1)], 'k-', 'LineWidth', 2); 
xlabel('Volume (mL)');
ylabel('Pressure (mmHg)');
title('Pressure–Volume Loop');
xlim([0 150]); ylim([0 160]);
saveas(gcf, "q2a.png")

%% Q2b

V = S.PVloop(1,:);  P = S.PVloop(2,:);
EDV = max(V); ESV = min(V); SV = EDV - ESV; 
MAP = 2*(87.785/3) + max(P)/3; % diastolic pressure is aortic valve opens
% systolic pressure is max pressure
SW_J = MAP * SV * 133.322e-6; % then convert to joules

fprintf('EDV=%.1f mL, ESV=%.1f mL, SV=%.1f mL\n',EDV,ESV,SV);
fprintf('Stroke work ≈ %.3f J per beat\n', SW_J);

%% Q2c

% Area in mmHg·mL and convert to Joules
SW_mmHg_mL = abs(trapz(V, P));    
SW_J = SW_mmHg_mL * 133.322e-6; 

% Plot and shade the loop
figure; hold on; box on;
fill(V, P, [0.9 0.9 1], 'EdgeColor', 'none'); % shaded area = work
plot(V, P, 'k-', 'LineWidth', 2);
xlabel('Volume (mL)'); ylabel('Pressure (mmHg)');
title(sprintf('PV Loop — Stroke Work = %.2f J', SW_J));
fprintf('Stroke Work = %.2f J per beat (area method)\n', SW_J);
saveas(gcf, "q2c.png")