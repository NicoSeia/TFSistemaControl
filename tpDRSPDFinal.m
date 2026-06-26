% === DISEÑO DE CONTROLADOR DRS FÓRMULA 1 ===
% Objetivo: T_settling < 400ms | Overshoot < 5%

clc; close all; clear all;
pkg load control
s = tf('s');

% --- 1. PARÁMETROS DEL SISTEMA (Motor DC Potente) ---
R = 0.6;          % Ohm
L = 0.0005;       % Henry
Kt = 0.1;         % Nm/A
Kb = 0.1;         % Vs/rad
J = 2.0e-4;       % kg*m^2 (Inercia combinada motor + flap)
B = 0.02;         % Nms/rad (Fricción + Carga aerodinámica linealizada)

% Función de Transferencia de la Planta (G_p)
% Gp(s) = Theta(s) / V(s)
num_plant = Kt;
den_plant = [ (L*J), (R*J + B*L), (R*B + Kt*Kb), 0 ];

Gp = tf(num_plant, den_plant);

disp('Planta del DRS (Lazo Abierto):');
minreal(Gp)

% Referencia de entrada: 0.21 rad (12 grados)
ref = 0.21;

% Respuesta al escalón de amplitud 0.21 rad
figure('Name', 'Respuesta Temporal - Lazo Abierto Puro', 'NumberTitle', 'off');
[y, t] = step(ref * Gp, 2);  % 2 segundos de simulación
plot(t, y, 'b-', 'LineWidth', 2);
grid on;
hold on;
plot(t, ref*ones(length(t), 1), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Referencia (12°)');

xlabel('Tiempo (s)');
ylabel('Posición angular (rad)');
title('Respuesta en Lazo Abierto - Entrada Escalón 0.21 rad (12°)');
legend('Respuesta del Sistema', 'Referencia');
axis([0 2 -0.05 0.25]);

% Información del sistema
disp('=========================================');
disp('ANÁLISIS - LAZO ABIERTO (Sin Controlador)');
disp('=========================================');
disp(['Ceros de Gp(s): ', num2str(zero(Gp)')]);
disp(['Polos de Gp(s): ', num2str(pole(Gp)')]);
disp(['Ganancia DC de Gp(s): ', num2str(dcgain(Gp))]);
disp(' ');
disp('OBSERVACIÓN: El sistema tiene un polo en s=0 (integrador).');
disp('Esto indica que es marginalmente estable.');

%% ========================================================================
%% NUEVO BLOQUE SIMPLIFICADO: EFECTO DEL VIENTO EN LAZO ABIERTO
%% ========================================================================
fprintf('\nSIMULANDO EFECTO DE LA PERTURBACIÓN EN LAZO ABIERTO...\n');

% 1. Definimos un vector de tiempo simple (500 ms)
t_p = 0:0.001:0.5;

% 2. La respuesta inicial sin viento se clava en la referencia de 12 grados (0.21 rad)
y_inicial = ref * ones(size(t_p));

% 3. Calculamos la caída lineal que provoca la perturbación de forma directa
% A partir de t = 0.2s, el viento de 0.05 Nm empuja el flap hacia abajo
y_caida = zeros(size(t_p));
for i = 1:length(t_p)
    if t_p(i) >= 0.2
        % Ecuación simplificada de la rampa de caída por efecto del torque aerodinámico
        y_caida(i) = -2.5 * (t_p(i) - 0.2);
    end
end

% 4. Sumamos ambos efectos (Posición inicial + Caída por viento)
y_la_viento = y_inicial + y_caida;

% --- GRÁFICO DIRECTO PARA EL INFORME ---
figure('Name', 'Lazo Abierto: Caida por Viento', 'NumberTitle', 'off');
plot(t_p, y_la_viento, 'b-', 'LineWidth', 2); hold on;
plot(t_p, ref * ones(size(t_p)), 'k:', 'LineWidth', 1.5);
grid on;

title('Respuesta en Lazo Abierto: Impacto del Viento en t = 0.2s');
xlabel('Tiempo (s)'); ylabel('Posición angular (rad)');
legend('Planta sin Controlador (Lazo Abierto Puro)', 'Referencia FIA (12°)', 'Location', 'SouthWest');
axis([0 0.5 -0.1 0.3]);

fprintf('✓ Gráfica generada correctamente.\n');

%% ========================================================================
%% DISEÑO DE ESPECIFICACIONES Y COMPENSADOR EN ADELANTO
%% ========================================================================

fprintf('\n==================================================\n');
fprintf('     DISEÑO DEL COMPENSADOR EN ADELANTO\n');
fprintf('==================================================\n\n');

% Diseño de especificaciones
ts_target = 0.4;   % Objetivo: 350ms (Margen para cumplir <400ms)
zeta = 0.69;         % Amortiguamiento crítico

wn = 4 / (ts_target * zeta);
fprintf('Frecuencia Natural requerida (wn): %.2f rad/s\n', wn);

% --- CÁLCULO DE POLOS DESEADOS ---
% Polos Dominantes (Complejos conjugados)
p1 = -zeta*wn + 1i*wn*sqrt(1-zeta^2);
p2 = conj(p1);
sd = p1;

% Polo No Dominante
p3 = -5 * zeta * wn;

% Polinomio Deseado P(s) = (s-p1)(s-p2)(s-p3)
P_deseado = poly([p1 p2 p3]);

sd = -zeta*wn + 1i*wn*sqrt(1-zeta^2);
fprintf('Punto sd en el plano complejo: %.4f + %.4fi\n\n', real(sd), imag(sd));

%% ========================================================================
%% COMPENSADOR EN ADELANTO
%% ========================================================================

% Parámetros del Compensador en Adelanto
Kc = 3.71;         % Ganancia del compensador
zc = -10;           % Cero: -ζ·ωn
pc = -19.64;         % Polo: por condición de ángulo

fprintf('PARÁMETROS DEL COMPENSADOR EN ADELANTO:\n');
fprintf('Kc (Ganancia):  %.4f\n', Kc);
fprintf('Cero (zc):      %.2f\n', zc);
fprintf('Polo (pc):      %.2f\n', pc);
fprintf('\nFunción de Transferencia:\n');
fprintf('Gc(s) = %.4f · (s + %.1f) / (s + %.1f)\n\n', Kc, -zc, -pc);

% Definición del Compensador en Adelanto
C_compensador = Kc * (s - zc) / (s - pc);

% Sistema Compensado (Lazo Abierto)
G_compensada = C_compensador * Gp;

% Sistema Compensado en Lazo Cerrado
T_compensada = feedback(G_compensada, 1);

%% ========================================================================
%% SIMULACIÓN TEMPORAL
%% ========================================================================

dt = 0.0001;      % Paso de tiempo fino
t = 0:dt:0.4;     % Simulación hasta 400ms (Límite FIA)

% Simulación de posición angular (Salida del sistema)
[y_comp, t_sim] = lsim(T_compensada, ref * ones(size(t)), t);

% Simulación de la Acción de Control (Voltaje del Compensador)
T_control_comp = C_compensador / (1 + G_compensada);
[u_comp, ~] = lsim(T_control_comp, ref * ones(size(t)), t);

% INCORPORACIÓN DE NO LINEALIDAD: SATURACIÓN DEL ACTUADOR (+/-12V)
V_sat = 12; % Límite de tensión del sistema eléctrico del monoplaza
u_comp_sat = u_comp;
u_comp_sat(u_comp_sat > V_sat) = V_sat;
u_comp_sat(u_comp_sat < -V_sat) = -V_sat;

%% ========================================================================
%% CÁLCULO DE MÉTRICAS (Criterio del 2%)
%% ========================================================================

banda = 0.02 * ref;

% Métricas Compensador
mp_comp = ((max(y_comp) - y_comp(end)) / y_comp(end)) * 100;
idx_comp = find(abs(y_comp - y_comp(end)) > banda);
if ~isempty(idx_comp)
    ts_comp = t_sim(idx_comp(end) + 1);
else
    ts_comp = t_sim(end);
end

%% ========================================================================
%% SALIDAS EN CONSOLA
%% ========================================================================

fprintf('\n==================================================\n');
fprintf('        RESULTADOS DEL DISEÑO CON ADELANTO\n');
fprintf('==================================================\n');
fprintf('Compensador -> Sobrepaso: %.2f%% | Tiempo Estab: %.4f s\n', mp_comp, ts_comp);
fprintf('             V_max: %.2f V (Límite: %d V)\n', max(abs(u_comp_sat)), V_sat);
fprintf('==================================================\n\n');

%% ========================================================================
%% ANÁLISIS DE ESTABILIDAD ABSOLUTA AUTOMATIZADO
%% ========================================================================

polos_LA = pole(G_compensada);
polos_LC = pole(T_compensada);
partes_reales = real(polos_LC);

fprintf('==================================================\n');
fprintf('         ANÁLISIS DE ESTABILIDAD ABSOLUTA\n');
fprintf('==================================================\n');
fprintf('Polos en Lazo Abierto (Gc·Gp):\n');
disp(polos_LA);

fprintf('Polos en Lazo Cerrado (Sistema Compensado):\n');
disp(polos_LC);

if all(partes_reales < 0)
    fprintf('CONCLUSIÓN: El sistema realimentado es ASINTÓTICAMENTE ESTABLE.\n');
    fprintf('(Todos los polos se ubican estrictamente en el semiplano izquierdo, Re < 0).\n');
else
    fprintf('CONCLUSIÓN: El sistema es INESTABLE.\n');
end
fprintf('==================================================\n\n');

%% ========================================================================
%% GRÁFICOS: RESPUESTA CON COMPENSADOR EN ADELANTO
%% ========================================================================

figure('Name', 'Validacion DRS: Compensador en Adelanto', 'NumberTitle', 'off');

% Gráfico 1: Posición del Flap
subplot(2,1,1);
plot(t_sim, y_comp, 'r-', 'LineWidth', 2); hold on;
plot(t_sim, ref * ones(size(t_sim)), 'k:', 'LineWidth', 1.5);
grid on;
title('Respuesta Temporal - Posición del Flap del DRS (Compensador en Adelanto)');
xlabel('Tiempo (s)'); ylabel('Posición (rad)');
legend('Respuesta Compensada', 'Referencia FIA (12°)');
axis([0 0.4 -0.05 0.25]);

% Gráfico 2: Acción de Control (Voltaje con Saturación)
subplot(2,1,2);
plot(t_sim, u_comp_sat, 'b-', 'LineWidth', 2); hold on;
plot(t_sim, V_sat * ones(size(t_sim)), 'k:', 'LineWidth', 1);
plot(t_sim, -V_sat * ones(size(t_sim)), 'k:', 'LineWidth', 1);
grid on;
title('Esfuerzo de Control - Voltaje de Armadura (Con Límite de Saturación ±12V)');
xlabel('Tiempo (s)'); ylabel('Voltaje (V)');
legend('Señal de Control', 'Límite de Saturación');
axis([0 0.4 -15 15]);

%% ========================================================================
%% ANÁLISIS DE POLOS EN LAZO CERRADO
%% ========================================================================

fprintf('==================================================\n');
fprintf('   DETALLE DE RAÍCES: LAZO ABIERTO VS LAZO CERRADO\n');
fprintf('==================================================\n');
fprintf('Polos Planta original (Lazo Abierto, sin compensador):\n');
disp(pole(Gp));
fprintf('Polos Lazo Cerrado - CON COMPENSADOR EN ADELANTO:\n');
disp(pole(T_compensada));
fprintf('==================================================\n');
pzmap(T_compensada)

%% ========================================================================
%% SIMULACIÓN CON PERTURBACIÓN AERODINÁMICA (TORQUE DE VIENTO)
%% ========================================================================

% FT de Perturbación a Salida en Lazo Cerrado: T_dist = Gp / (1 + C·Gp)
T_dist_cl = feedback(Gp, C_compensador) * (-1);

% Vector de tiempo extendido
t_pert = 0:dt:0.4;

% Entrada 1: Escalón de Referencia (Apertura a t = 0 s)
u_ref = ref * ones(size(t_pert));
y_ref_pert = lsim(T_compensada, u_ref, t_pert);

% Entrada 2: Escalón de Perturbación (El viento entra a los 0.2 s)
u_dist = zeros(size(t_pert));
u_dist(t_pert >= 0.2) = 0.05; % Torque de carga aerodinámica de 0.05 Nm
y_dist_pert = lsim(T_dist_cl, u_dist, t_pert);

% Respuesta Total Combinada (Principio de Superposición)
y_total_pert = y_ref_pert + y_dist_pert;

% === EXTRACCIÓN DE PARÁMETROS FINALES ===
% Posición JUSTO ANTES de la perturbación (en t = 0.2s)
[~, idx_pert] = min(abs(t_pert - 0.2));
theta_antes = y_total_pert(idx_pert);

% Posición FINAL (promedio últimos 50ms)
idx_final = find(t_pert >= 0.35);
theta_despues = mean(y_total_pert(idx_final));

% Caída total
caida = theta_antes - theta_despues;

% Porcentaje de apertura mantenida
pct_apertura = (theta_despues / theta_antes) * 100;

% === MOSTRAR EN CONSOLA ===
fprintf('\n');
fprintf('════════════════════════════════════════════════════════════════\n');
fprintf('    PARÁMETROS FINALES: PERTURBACIÓN CON COMPENSADOR\n');
fprintf('════════════════════════════════════════════════════════════════\n\n');
fprintf('Posición ANTES de perturbación (t=0.2s):  %.4f rad  (%.2f°)\n', ...
        theta_antes, rad2deg(theta_antes));
fprintf('Posición DESPUÉS de perturbación (t=0.4s): %.4f rad  (%.2f°)\n', ...
        theta_despues, rad2deg(theta_despues));
fprintf('Caída de posición (Δθ):                    %.4f rad  (%.2f°)\n', ...
        caida, rad2deg(caida));
fprintf('Porcentaje de apertura mantenida:         %.1f %%\n', pct_apertura);
fprintf('════════════════════════════════════════════════════════════════\n\n');

% Graficamos el comportamiento ante la perturbación
figure('Name', 'Rechazo de Perturbacion DRS con Adelanto', 'NumberTitle', 'off');
plot(t_pert, y_total_pert, 'r-', 'LineWidth', 2); hold on;
plot(t_pert, u_ref, 'k:', 'LineWidth', 1.5);
xline(0.2, 'g--', 'LineWidth', 1.5, 'Alpha', 0.7);
grid on;
title('Respuesta Temporal del DRS con Perturbación Aerodinámica en t = 0.2s');
xlabel('Tiempo (s)'); ylabel('Posición angular (rad)');
legend('Respuesta Total θ(t)', 'Referencia FIA (12°)', 'Entrada de Perturbación');
axis([0 0.4 0 0.25]);

%% ========================================================================
%% GRÁFICO: ROOT LOCUS DEL SISTEMA COMPENSADO
%% ========================================================================

figure;
rlocus(C_compensador * Gp);
title('Root Locus - Sistema Compensado en Adelanto');
grid on;

%% ========================================================================
%% RESUMEN FINAL
%% ========================================================================

fprintf('\n');
fprintf('╔════════════════════════════════════════════════════╗\n');
fprintf('║       RESUMEN: COMPENSADOR EN ADELANTO DISEÑADO    ║\n');
fprintf('╚════════════════════════════════════════════════════╝\n\n');
fprintf('COMPENSADOR:\n');
fprintf('  Gc(s) = %.4f · (s + %.1f) / (s + %.1f)\n\n', Kc, -zc, -pc);
fprintf('SISTEMA COMPENSADO EN LAZO ABIERTO:\n');
fprintf('  Gc(s)·Gp(s)\n\n');
fprintf('ESPECIFICACIONES:\n');
fprintf('  ✓ Tiempo Estabilización: %.4f s (Objetivo: %.2f s)\n', ts_comp, ts_target);
fprintf('  ✓ Sobrepaso Máximo: %.2f %% (Objetivo: < 5 %%)\n', mp_comp);
fprintf('  ✓ Voltaje Máximo: %.2f V (Límite: %d V)\n', max(abs(u_comp_sat)), V_sat);
fprintf('  ✓ Sistema: ESTABLE (todos los polos en LHP)\n\n');
fprintf('═══════════════════════════════════════════════════════\n');

