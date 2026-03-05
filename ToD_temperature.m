%% === 95%-Konfidenzintervall der geschätzten "time since death" (T0 = 33.3 °C) ===
% nutzt T_measured, TU und die Fit-Parameter mit Kovarianzmatrix

T0 = 33.3;  % fixer Startwert (lebende Stirn)
alpha = 0.4728;  beta = 0.1923;
gamma = -0.5789;  delta = 0.0776;

% --- 1) Nominale Todeszeit t* bestimmen ---
f_T = @(tt) ((alpha*exp(-beta*tt) - gamma*exp(-delta*tt))*(T0 - TU) + TU) - T_measured;
t_star = fzero(f_T, [0, max(tFit)]);   % h

% --- 2) Ableitungen des Modells am Punkt t* ---
e_bt = exp(-beta*t_star);
e_dt = exp(-delta*t_star);

% partielle Ableitungen nach Parametern
dT_da =  e_bt        * (T0 - TU);
dT_db = -alpha*t_star*e_bt * (T0 - TU);
dT_dg = -e_dt        * (T0 - TU);
dT_dd =  gamma*t_star*e_dt * (T0 - TU);

% partielle Ableitung nach der Zeit t
dT_dt = (-alpha*beta*e_bt + gamma*delta*e_dt) * (T0 - TU);

% Gradient nach Parametern
g = [dT_da, dT_db, dT_dg, dT_dd];  % Zeilenvektor

% --- 3) Varianzpropagation (Delta-Methode) ---
var_t = (g * covMatrix * g.') / (dT_dt^2);
se_t  = sqrt(max(var_t,0));

% --- 4) 95%-Konfidenzintervall für t* ---
n   = length(yData);
k   = length(params);
tval = tinv(0.975, max(n-k,1));
CI_t = [t_star - tval*se_t, t_star + tval*se_t];

% --- 5) Ausgabe ---
fprintf('\n=== 95%%-Konfidenzintervall für "time since death" (T0=33.3°C) ===\n');
fprintf('Nominale Schätzung: %.3f h\n', t_star);
fprintf('Standardfehler: %.3f h\n', se_t);
fprintf('95%%-CI: [%.3f, %.3f] h\n', CI_t(1), CI_t(2));
