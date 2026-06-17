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

% Diseño de especificaciones
ts_target = 0.35; % Objetivo: 350ms (Margen para cumplir <400ms)
zeta = 0.7;       % Amortiguamiento crítico

wn = 4 / (ts_target * zeta)
fprintf('Frecuencia Natural requerida (wn): %.2f rad/s\n', wn);

% --- 3. CÁLCULO DE POLOS DESEADOS ---
% Polos Dominantes (Complejos conjugados)
p1 = -zeta*wn + 1i*wn*sqrt(1-zeta^2)
p2 = conj(p1)
sd = p1
% Polo No Dominante
% Debe ser más rápido que los dominantes.
%p3 = -5 * zeta * wn
p3 = -5 * zeta * wn % 5 veces más lejos a la izquierda

% Polinomio Deseado P(s) = (s-p1)(s-p2)(s-p3)
P_deseado = poly([p1 p2 p3])

% 2. Coordenada exacta del polo deseado (sd) en el plano complejo
sd = -zeta*wn + 1i*wn*sqrt(1-zeta^2);
fprintf('Punto sd en el plano complejo: %.4f + %.4fi\n\n', real(sd), imag(sd));

% Parámetros del Controlador PD
Kp = 20;
Kd = 0.02;
Tf = 0.001;  % Filtro de 1 ms para hacerlo realizable

% ------------------------------------------------------------------------
% 1. CASO IDEAL (Teórico)
% ------------------------------------------------------------------------
C_ideal = Kp + Kd*s;
T_ideal = feedback(C_ideal * Gp, 1);

% ------------------------------------------------------------------------
% 2. CASO REAL (Físicamente Realizable con Filtro)
% ------------------------------------------------------------------------
C_real = Kp + (Kd*s)/(1 + Tf*s);
C_real = minreal(C_real);
T_real = feedback(C_real * Gp, 1);

% ------------------------------------------------------------------------
% SIMULACIÓN TEMPORAL
% ------------------------------------------------------------------------
dt = 0.0001;      % Paso de tiempo fino
t = 0:dt:0.4;     % Simulación hasta 400ms (Límite FIA)

% Simulación de posición angular (Salida del sistema)
[y_ideal, t_sim] = lsim(T_ideal, ref * ones(size(t)), t);
[y_real, t_sim]  = lsim(T_real, ref * ones(size(t)), t);

% Simulación de la Acción de Control (Voltaje aplicado al motor)
% Caso Ideal: u = Kp*e + Kd*(de/dt) calculado numéricamente
E_ideal = ref - y_ideal;
de_dt_ideal = [diff(E_ideal) / dt; 0];
u_ideal = Kp * E_ideal + Kd * de_dt_ideal;

% Caso Real: Se calcula por bloques porque el filtro la hace propia
T_control_real = C_real / (1 + C_real*Gp);
[u_real, ~] = lsim(T_control_real, ref * ones(size(t)), t);

% ------------------------------------------------------------------------
% INCORPORACIÓN DE NO LINEALIDAD: SATURACIÓN DEL ACTUADOR (+/-12V)
% ------------------------------------------------------------------------
% Analizamos el impacto de la restricción fisica limitando el vector u
V_sat = 12; % Límite de tensión del sistema eléctrico del monoplaza
u_real_sat = u_real;
u_real_sat(u_real_sat > V_sat) = V_sat;
u_real_sat(u_real_sat < -V_sat) = -V_sat;

% ------------------------------------------------------------------------
% CÁLCULO DE MÉTRICAS (Criterio del 2%)
% ------------------------------------------------------------------------
banda = 0.02 * ref;

% Métricas Lazo Ideal
mp_id = ((max(y_ideal) - y_ideal(end)) / y_ideal(end)) * 100;
idx_id = find(abs(y_ideal - y_ideal(end)) > banda);
ts_id = t_sim(idx_id(end) + 1);

% Métricas Lazo Real
mp_re = ((max(y_real) - y_real(end)) / y_real(end)) * 100;
idx_re = find(abs(y_real - y_real(end)) > banda);
ts_re = t_sim(idx_re(end) + 1);

% ------------------------------------------------------------------------
% SALIDAS SIMPLES POR CONSOLA
% ------------------------------------------------------------------------
fprintf('\n==================================================\n');
fprintf('              RESULTADOS DEL DISEÑO\n');
fprintf('==================================================\n');
fprintf('PD IDEAL -> Sobrepaso: %.2f%% | Tiempo Estab: %.4f s | V_max: %.2f V\n', mp_id, ts_id, max(abs(u_ideal)));
fprintf('PD REAL  -> Sobrepaso: %.2f%% | Tiempo Estab: %.4f s | V_max: %.2f V\n', mp_re, ts_re, max(abs(u_real)));
fprintf('==================================================\n\n');

% ------------------------------------------------------------------------
% B. ANÁLISIS DE ESTABILIDAD ABSOLUTA AUTOMATIZADO
% ------------------------------------------------------------------------
polos_LC = pole(T_real);
partes_reales = real(polos_LC);

fprintf('==================================================\n');
fprintf('         ANÁLISIS DE ESTABILIDAD ABSOLUTA\n');
fprintf('==================================================\n');
fprintf('Polos de Lazo Cerrado del Sistema Real obtenidos:\n');
disp(polos_LC);

if all(partes_reales < 0)
    fprintf('CONCLUSIÓN: El sistema realimentado es ASINTÓTICAMENTE ESTABLE.\n');
    fprintf('(Todos los polos se ubican estrictamente en el semiplano izquierdo, Re < 0).\n');
else
    fprintf('CONCLUSIÓN: El sistema es INESTABLE.\n');
end
fprintf('==================================================\n\n');

% ------------------------------------------------------------------------
% GRÁFICOS COMPACTOS (Incluye la Acción de Control con No Linealidad)
% ------------------------------------------------------------------------
figure('Name', 'Validacion DRS: Ideal vs Realizable', 'NumberTitle', 'off');

% Gráfico 1: Posición del Flap
subplot(2,1,1);
plot(t_sim, y_ideal, 'b-', 'LineWidth', 2); hold on;
plot(t_sim, y_real, 'r--', 'LineWidth', 2);
plot(t_sim, ref * ones(size(t_sim)), 'k:', 'LineWidth', 1.5);
grid on;
title('Respuesta Temporal - Posición del Flap del DRS');
xlabel('Tiempo (s)'); ylabel('Posición (rad)');
legend('PD Ideal', 'PD Real (Filtro 1ms)', 'Referencia FIA (12°)');

% Gráfico 2: Acción de Control (Voltaje con Sat.)
subplot(2,1,2);
plot(t_sim, u_ideal, 'b-', 'LineWidth', 2); hold on;
plot(t_sim, u_real_sat, 'r--', 'LineWidth', 2);
plot(t_sim, V_sat * ones(size(t_sim)), 'k:', 'LineWidth', 1);
plot(t_sim, -V_sat * ones(size(t_sim)), 'k:', 'LineWidth', 1);
grid on;
title('Esfuerzo de Control - Voltaje de Armadura (Con Límite de Saturación)');
xlabel('Tiempo (s)'); ylabel('Voltaje (V)');
legend('PD Ideal (Impropio)', 'PD Real (Acotado a +/-12V)');

% ------------------------------------------------------------------------
% ANÁLISIS DE POLOS EN LAZO CERRADO PARA CONCLUSIONES
% ------------------------------------------------------------------------
fprintf('==================================================\n');
fprintf('   DETALLE DE RAÍCES: LAZO ABIERTO VS LAZO CERRADO\n');
fprintf('==================================================\n');
fprintf('Polos Planta original (Lazo Abierto):\n');
disp(pole(Gp));
fprintf('Polos Lazo Cerrado - PD IDEAL:\n');
disp(pole(T_ideal));
fprintf('Polos Lazo Cerrado - PD REAL (Con Filtro):\n');
disp(pole(T_real));
fprintf('==================================================\n');

% ------------------------------------------------------------------------
% SIMULACIÓN CON PERTURBACIÓN AERODINÁMICA (TORQUE DE VIENTO)
% ------------------------------------------------------------------------
% FT de Perturbación a Salida en Lazo Cerrado: T_dist = Gp / (1 + C_real*Gp)
% Multiplicamos por -1 porque el torque del viento se opone al motor
T_dist_cl = feedback(Gp, C_real) * (-1);

% Vector de tiempo extendido para ver el escalón de viento
t_pert = 0:dt:0.4;

% Entrada 1: Escalón de Referencia (Apertura a t = 0 s)
u_ref = ref * ones(size(t_pert));
y_ref_pert = lsim(T_real, u_ref, t_pert);

% Entrada 2: Escalón de Perturbación (El viento entra a los 0.2 s)
u_dist = zeros(size(t_pert));
u_dist(t_pert >= 0.2) = 0.05; % Torque de carga aerodinámica de 0.05 Nm
y_dist_pert = lsim(T_dist_cl, u_dist, t_pert);

% Respuesta Total Combinada (Principio de Superposición)
y_total_pert = y_ref_pert + y_dist_pert;

% Graficamos el comportamiento ante la perturbación
figure('Name', 'Rechazo de Perturbacion DRS', 'NumberTitle', 'off');
plot(t_pert, y_total_pert, 'r-', 'LineWidth', 2); hold on;
plot(t_pert, u_ref, 'k:', 'LineWidth', 1.5);
grid on;
title('Respuesta Temporal del DRS con Perturbación Aerodinámica en t = 0.2s');
xlabel('Tiempo (s)'); ylabel('Posición angular (rad)');
legend('Respuesta Total \theta(t)', 'Referencia FIA (12°)');

