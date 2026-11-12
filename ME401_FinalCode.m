%% ME 401 - 1970 Plymouth Road Runner 383 V8
% Testing different cycles to optimize engine efficiency
% Units: pressure [psi], volume [in^3], temperature [R], mass [lbm]
% Work/heat [in-lbf]
%
% Assumptions (from mind map): working fluid is air (ideal gas),
% ideal cycles, no friction/heat losses, closed system, uniform properties.

clear; clc; close all;

%% ---------------- Colors (no black) ----------------
% Dark navy background so the white axes are readable.
cBg      = [0.10 0.14 0.25];   % navy background (figure and axes)
cWhite   = [1 1 1];            % axes, tick labels, grid, text
cOtto    = [0.30 0.75 1.00];   % Fig 1: cyan-blue
cDiesel  = [1.00 0.55 0.20];   % Fig 2: orange
cBrayton = [0.45 0.90 0.45];   % Fig 3: green
cCarnot  = [0.90 0.50 1.00];   % Fig 4: magenta
cBars    = [0.70 0.72 0.76];   % Fig 5: all bars the same light grey
set(groot,'defaultFigureColor',cBg,'defaultAxesColor',cBg, ...
          'defaultFigureInvertHardcopy','off', ...
          'defaultAxesXColor',cWhite,'defaultAxesYColor',cWhite, ...
          'defaultAxesGridColor',cWhite,'defaultTextColor',cWhite);

%% ---------------- Inputs ----------------
Pmin = 14.7;     % psi  (state 1 pressure, atmospheric)
Pmax = 150;      % psi  (Brayton high pressure; edit as needed)
Tmin = 540;      % R    (state 1 temperature, ~80 F)
Tmax = 3000;     % R    (peak cycle temperature)
Vmax = 383;      % in^3 (engine displacement)
r    = 9.4;      % compression ratio (9.4:1)
k    = 1.4;      % specific heat ratio of air
R    = 53.35*12; % gas constant of air, in-lbf/(lbm-R)  (53.35 ft-lbf/lbm-R)

cv = R/(k-1);    % specific heats, in-lbf/(lbm-R)
cp = k*R/(k-1);

%% ---------------- Calculations ----------------
Vmin = Vmax/r;                 % in^3
% Mass of air: evaluated at state 1 (bottom dead center: Pmin, Vmax, Tmin).
% NOTE: your pseudocode used Pmax, Vmax, Tmax; those are not a single state
% of the cycle, so state 1 is used here. See the note at the bottom.
m = Pmin*Vmax/(R*Tmin);        % lbm

%% ======== OTTO CYCLE ========
% 1->2 isentropic compression, 2->3 constant-volume heat addition,
% 3->4 isentropic expansion, 4->1 constant-volume heat rejection
V1 = Vmax;  P1 = Pmin;  T1 = Tmin;
V2 = Vmin;
P2 = P1*(V1/V2)^k;                 % find P2
T2 = T1*(V1/V2)^(k-1);
V3 = V2;    T3 = Tmax;
P3 = P2*T3/T2;                     % constant volume, ideal gas
V4 = Vmax;
P4 = P3*(V3/V4)^k;                 % find P4

eta_otto = 1 - 1/r^(k-1);          % efficiency
Qin_otto  = m*cv*(T3 - T2);
Wout_otto = eta_otto*Qin_otto;

figure(1); hold on; grid on;
Va = linspace(V1,V2,100);  plot(Va, P1*(V1./Va).^k, 'Color', cOtto, 'LineWidth', 2);
plot([V2 V3],[P2 P3], 'Color', cOtto, 'LineWidth', 2);
Vb = linspace(V3,V4,100);  plot(Vb, P3*(V3./Vb).^k, 'Color', cOtto, 'LineWidth', 2);
plot([V4 V1],[P4 P1], 'Color', cOtto, 'LineWidth', 2);
text([V1 V2 V3 V4],[P1 P2 P3 P4],{' 1',' 2',' 3',' 4'}, 'Color', cOtto, 'FontWeight','bold');
xlabel('Volume (in^3)'); ylabel('Pressure (psi)');
title(sprintf('Otto Cycle P-V Diagram (\\eta = %.1f%%)', 100*eta_otto), 'Color', cOtto);

%% ======== DIESEL CYCLE ========
% 1->2 isentropic compression, 2->3 constant-pressure heat addition,
% 3->4 isentropic expansion, 4->1 constant-volume heat rejection
D1_V = Vmax; D1_P = Pmin; D1_T = Tmin;
D2_V = Vmin;
D2_P = D1_P*(D1_V/D2_V)^k;
D2_T = D1_T*(D1_V/D2_V)^(k-1);
D3_T = Tmax;
D3_P = D2_P;                       % constant pressure (see note on Pmax)
D3_V = m*R*D3_T/D3_P;              % V3 = mRT3/P3
rc   = D3_V/Vmin;                  % cutoff ratio
D4_V = Vmax;
D4_P = D3_P*(D3_V/D4_V)^k;         % find P4

if D3_V >= Vmax || rc <= 1
    warning('Diesel: cutoff volume V3 = %.1f in^3 is not between Vmin and Vmax. Check Tmax.', D3_V);
end

eta_diesel = 1 - (1/r^(k-1)) * (rc^k - 1)/(k*(rc - 1));
Qin_diesel  = m*cp*(D3_T - D2_T);
Wout_diesel = eta_diesel*Qin_diesel;

figure(2); hold on; grid on;
Va = linspace(D1_V,D2_V,100); plot(Va, D1_P*(D1_V./Va).^k, 'Color', cDiesel, 'LineWidth', 2);
plot([D2_V D3_V],[D2_P D3_P], 'Color', cDiesel, 'LineWidth', 2);
Vb = linspace(D3_V,D4_V,100); plot(Vb, D3_P*(D3_V./Vb).^k, 'Color', cDiesel, 'LineWidth', 2);
plot([D4_V D1_V],[D4_P D1_P], 'Color', cDiesel, 'LineWidth', 2);
text([D1_V D2_V D3_V D4_V],[D1_P D2_P D3_P D4_P],{' 1',' 2',' 3',' 4'}, 'Color', cDiesel, 'FontWeight','bold');
xlabel('Volume (in^3)'); ylabel('Pressure (psi)');
title(sprintf('Diesel Cycle P-V Diagram (\\eta = %.1f%%)', 100*eta_diesel), 'Color', cDiesel);

%% ======== BRAYTON CYCLE ========
% 1->2 isentropic compression, 2->3 constant-pressure heat addition,
% 3->4 isentropic expansion, 4->1 constant-pressure heat rejection
rp   = Pmax/Pmin;                          % pressure ratio
B1_P = Pmin; B2_P = Pmax; B3_P = Pmax; B4_P = Pmin;
B2_V = Vmin;                               % V2 = Vmin
B1_V = B2_V*(B2_P/B1_P)^(1/k);             % calculate V1
B4_V = Vmax;                               % V4 = Vmax
B3_V = B4_V*(B4_P/B3_P)^(1/k);             % calculate V3

if B3_V <= B2_V
    warning('Brayton: V3 <= V2, so no heat is added. Lower Pmax (rp must be < r^(k/2)).');
end

eta_brayton = 1 - 1/rp^((k-1)/k);
B2_T = Tmin*rp^((k-1)/k);
Qin_brayton  = m*cp*(Tmax - B2_T);
Wout_brayton = eta_brayton*Qin_brayton;

figure(3); hold on; grid on;
Va = linspace(B1_V,B2_V,100); plot(Va, B1_P*(B1_V./Va).^k, 'Color', cBrayton, 'LineWidth', 2);
plot([B2_V B3_V],[B2_P B3_P], 'Color', cBrayton, 'LineWidth', 2);
Vb = linspace(B3_V,B4_V,100); plot(Vb, B3_P*(B3_V./Vb).^k, 'Color', cBrayton, 'LineWidth', 2);
plot([B4_V B1_V],[B4_P B1_P], 'Color', cBrayton, 'LineWidth', 2);
text([B1_V B2_V B3_V B4_V],[B1_P B2_P B3_P B4_P],{' 1',' 2',' 3',' 4'}, 'Color', cBrayton, 'FontWeight','bold');
xlabel('Volume (in^3)'); ylabel('Pressure (psi)');
title(sprintf('Brayton Cycle P-V Diagram (\\eta = %.1f%%)', 100*eta_brayton), 'Color', cBrayton);

%% ======== CARNOT CYCLE ========
% Between TH = Tmax and TL = Tmin.
% 2->1 isothermal expansion (TH), 1->4 isentropic expansion,
% 4->3 isothermal compression (TL), 3->2 isentropic compression
TH = Tmax;  TL = Tmin;
eta_carnot = 1 - TL/TH;

C4_V = Vmax;  C4_P = Pmin;                         % bottom dead center, at TL
C1_V = C4_V*(TL/TH)^(1/(k-1));                     % isentropic 1->4 (solved for V1)
C1_P = C4_P/((TL/TH)^(k/(k-1)));                   % P1
C3_V = Vmax/r;                                     % isothermal compression ratio = r
C3_P = C4_P*C4_V/C3_V;
C2_V = C3_V*(TL/TH)^(1/(k-1));                     % find V2
C2_P = C3_P/((TL/TH)^(k/(k-1)));                   % find P2

Qin_carnot  = m*R*TH*log(C1_V/C2_V);               % isothermal expansion at TH
Wout_carnot = eta_carnot*Qin_carnot;

figure(4); hold on; grid on;
Vi = linspace(C2_V,C1_V,100); plot(Vi, C2_P*C2_V./Vi, 'Color', cCarnot, 'LineWidth', 2);      % 2->1
Vi = linspace(C1_V,C4_V,100); plot(Vi, C1_P*(C1_V./Vi).^k, 'Color', cCarnot, 'LineWidth', 2); % 1->4
Vi = linspace(C4_V,C3_V,100); plot(Vi, C4_P*C4_V./Vi, 'Color', cCarnot, 'LineWidth', 2);      % 4->3
Vi = linspace(C3_V,C2_V,100); plot(Vi, C3_P*(C3_V./Vi).^k, 'Color', cCarnot, 'LineWidth', 2); % 3->2
text([C1_V C2_V C3_V C4_V],[C1_P C2_P C3_P C4_P],{' 1',' 2',' 3',' 4'}, 'Color', cCarnot, 'FontWeight','bold');
% Between Tmin and Tmax the isentropes span a huge range of V and P, so on linear
% axes the cycle is squashed against the axes. On log-log axes isotherms
% (PV = const) and isentropes (PV^k = const) are straight lines, so the
% cycle shows up as a box-like parallelogram.
set(gca,'XScale','log','YScale','log'); grid on;
xlabel('Volume (in^3), log scale'); ylabel('Pressure (psi), log scale');
title(sprintf('Carnot Cycle P-V Diagram (\\eta = %.1f%%)', 100*eta_carnot), 'Color', cCarnot);

%% ======== EFFICIENCY COMPARISON ========
names = {'Otto','Diesel','Brayton','Carnot'};
etas  = 100*[eta_otto eta_diesel eta_brayton eta_carnot];

figure(5);
hb = bar(etas); grid on;
set(hb,'FaceColor',cBars);
set(gca,'XTickLabel',names);
ylabel('Thermal efficiency (%)');
title('Cycle Efficiency Comparison', 'Color', cWhite);
labels = arrayfun(@(e) sprintf('%.1f%%', e), etas, 'UniformOutput', false);  % works in all MATLAB versions
text(1:4, etas, labels, 'HorizontalAlignment','center','VerticalAlignment','bottom');

%% ======== OUTPUTS ========
Qin  = [Qin_otto  Qin_diesel  Qin_brayton  Qin_carnot];
Wout = [Wout_otto Wout_diesel Wout_brayton Wout_carnot];
fprintf('\nMass of air m = %.4f lbm, Vmin = %.1f in^3\n\n', m, Vmin);
fprintf('%-9s %12s %14s %14s\n','Cycle','Efficiency','Qin (in-lbf)','Wout (in-lbf)');
for i = 1:4
    fprintf('%-9s %11.1f%% %14.1f %14.1f\n', names{i}, etas(i), Qin(i), Wout(i));
end

%% ---------------- Notes on changes from the pseudocode ----------------
% 1) m = Pmax*Vmax/(R*Tmax) mixes states. Mass is found at state 1 (Pmin, Vmax, Tmin).
% 2) Diesel P4: pseudocode used (V3/Vmin); isentropic expansion goes V3 -> V4 = Vmax,
%    so it is (V3/V4).
% 3) Diesel P3: constant-pressure heating means P3 = P2 (from the isentropic compression),
%    so Pmax is only used for Brayton.
% 4) Otto P3 was not in the pseudocode; it comes from P3 = P2*T3/T2 (constant volume).
% 5) Carnot state points are built from V4 = Vmax and an isothermal ratio equal to r.
% 6) The Carnot P-V diagram uses log-log axes so it is readable (see comment in that section).
