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

pole(Gp)
step(Gp)

figure(1);
rlocus(Gp);
title('Lugar de las Raíces - Sistema Original (Sin Compensar)');
grid on;

% --- ANÁLISIS DE RESPUESTA TEMPORAL (SIN CONTROLADOR) ---
T = feedback(Gp, 1);

% Simulamos la respuesta al escalón.
[y_sin, t_sin] = step(T);

% Como el sistema es lineal, multiplicamos la salida unitaria por 12
% para simular la referencia real de apertura del DRS.
Referencia_Grados = 12;
y_deg = y_sin * Referencia_Grados;
y_final_deg = y_deg(end);

% 1. Sobrepaso (Overshoot)
[y_max_deg, ~] = max(y_deg);
if y_max_deg > y_final_deg
    Overshoot_sin = ((y_max_deg - y_final_deg) / y_final_deg) * 100;
else
    Overshoot_sin = 0;
end

% 2. Tiempo de Establecimiento (Criterio del 2%)
margen_sup = 1.02 * y_final_deg;
margen_inf = 0.98 * y_final_deg;

% Buscamos índices donde la señal se sale de la banda del 2%
idx_fuera = find(y_deg > margen_sup | y_deg < margen_inf);

if isempty(idx_fuera)
    SettlingTime_sin = 0;
else
    % El tiempo de establecimiento es el instante del último punto fuera de banda
    SettlingTime_sin = t_sin(idx_fuera(end));
end

% --- GRÁFICA REESCALADA ---
figure;
plot(t_sin, y_deg, 'r-', 'LineWidth', 2); hold on;
% Dibujamos la línea de referencia en 12 grados
plot([t_sin(1) t_sin(end)], [12 12], 'k--', 'LineWidth', 1.5);
hold off;

grid on;
title('Respuesta Temporal SIN Controlador (Referencia: 12°)');
xlabel('Tiempo (s)');
ylabel('Posición Angular (Grados)');
legend('Respuesta del Sistema', 'Objetivo (12°)', 'Location', 'SouthEast');

% --- MOSTRAR RESULTADOS ---
fprintf('\n--- ANÁLISIS DEL SISTEMA ORIGINAL (SIN CONTROLADOR) ---\n');
fprintf('Valor Final: %.4f grados (Debería llegar a 12)\n', y_final_deg);
fprintf('Sobrepaso: %.2f %%\n', Overshoot_sin);
fprintf('Tiempo de Establecimiento (2%%): %.4f s\n', SettlingTime_sin);
fprintf('------------------------------------------------------\n\n');


% --- 2. REQUISITOS DE DISEÑO ---
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
p3 = -5 * real(p1);% 5 veces más lejos a la izquierda

% Polinomio Deseado P(s) = (s-p1)(s-p2)(s-p3)
P_deseado = poly([p1 p2 p3])

% --- 4. CÁLCULO DEL CONTROLADOR (MÉTODO DE FASE EXACTO) ---
% Estructura: C(s) = K * (s + z_c) / (s + p_c)

% A. Ubicación del Cero (z_c)
% Cancelamos el polo mecánico dominante (-200)
polos_planta = roots(den_plant);
z_c = abs(polos_planta(2)); % z_c = 200

% B. Ubicación del Polo (p_c) usando Condición de Fase
% Ángulo(C*Gp) = -180° en s=sd
% Ángulo(s+z) - Ángulo(s+p) + Ángulo(Gp) = -180°
% Ángulo(s+p) = 180° + Ángulo(Gp) + Ángulo(s+z)

% Evaluamos Gp en sd
Gp_en_sd = polyval(num_plant, sd) / polyval(den_plant, sd);
angulo_Gp = angle(Gp_en_sd);

% Evaluamos el término del cero (sd + zc)
angulo_z = angle(sd + z_c);

% Calculamos el ángulo necesario del polo
phi_polo = pi + angulo_Gp + angulo_z;

% Geometría básica para hallar la posición del polo p_c en el eje real
% tan(phi) = Imag(sd) / (Real(sd) - (-p_c))
% p_c = Real(sd) - Imag(sd) / tan(phi_polo)
p_c_val = real(sd) - imag(sd) / tan(phi_polo);

% Como p_c_val es la coordenada negativa (ej -1200),
% el denominador será (s - p_c_val) -> (s + abs(p_c_val))
p_c = abs(p_c_val);

fprintf('Controlador Recalculado (Método Fase) -> Cero: %.2f, Polo: %.2f\n', z_c, p_c);

% Definimos el Controlador Base
% Nota: (s + p_c) asegura que el polo sea estable
C_base = (s + z_c) / (s + p_c)

% --- 5. CÁLCULO DE K ---

% 1. Definimos el Polo Objetivo Exacto (s_d)
sd = -zeta*wn + 1i*wn*sqrt(1-zeta^2)

% 2. Extraemos los vectores de numerador y denominador
% tfdata(Sistema, 'v') devuelve los vectores numéricos
[num_Gp, den_Gp] = tfdata(Gp, 'v');
[num_C, den_C] = tfdata(C_base, 'v');

% 3. Evaluamos los polinomios en el punto complejo 'sd' usando polyval
val_Gp_complex = polyval(num_Gp, sd) / polyval(den_Gp, sd);
val_C_base_complex = polyval(num_C, sd) / polyval(den_C, sd);

% 4. Aplicamos la Condición de Magnitud: K = 1 / |L(s)|
K_calculada = 1 / abs(val_C_base_complex * val_Gp_complex)

fprintf('Ganancia K calculada analíticamente: %.4f\n', K_calculada);

figure(3);
rlocus(Controlador_Final * Gp);
hold on;
plot(real(sd), imag(sd), 'r*', 'MarkerSize', 10, 'LineWidth', 2); % Marca el polo deseado
hold off;
title(['Lugar de las Raíces - Sistema COMPENSADO (K = ' num2str(K_calculada) ')']);
grid on;
legend('LGR Compensado', 'Polo Deseado (Target)');
% Aquí verás que las ramas ahora pasan exactamente por el punto rojo (tu objetivo)

% --- SIMULACIÓN FINAL ---
% 1. Definición del Controlador Final
Controlador_Final = K_calculada * C_base;

% 2. Lazo Cerrado
T_final = feedback(Controlador_Final * Gp, 1);

% ----------------------------------------------------------
% A. SIMULACIÓN DE RESPUESTA AL ESCALÓN
% ----------------------------------------------------------
t_sim = 0:0.001:1.5;
[y_norm, t] = step(T_final, t_sim);

% --- ESCALADO A GRADOS ---
Ref_Grados = 12;
y_deg = y_norm * Ref_Grados;
y_final_deg = y_deg(end);

% --- CÁLCULO DE MÉTRICAS ---
% 1. Pico Máximo (Overshoot)
[val_max, idx_max] = max(y_deg);
t_peak = t(idx_max);
Overshoot_pct = ((val_max - y_final_deg) / y_final_deg) * 100;

% 2. Tiempo de Establecimiento (2%)
margen_sup = 1.02 * y_final_deg;
margen_inf = 0.98 * y_final_deg;
idx_fuera = find(y_deg > margen_sup | y_deg < margen_inf);

if isempty(idx_fuera)
    ts = 0;
else
    ts = t(idx_fuera(end));
end

% --- GRÁFICA: RESPUESTA TEMPORAL DETALLADA ---
figure;
plot(t, y_deg, 'b-', 'LineWidth', 2); hold on;

% Referencia
plot([t(1) t(end)], [Ref_Grados Ref_Grados], 'k--', 'LineWidth', 1.5);

% Banda del 2% (Visualización de tolerancia)
plot([t(1) t(end)], [margen_sup margen_sup], 'g:', 'LineWidth', 1);
plot([t(1) t(end)], [margen_inf margen_inf], 'g:', 'LineWidth', 1);

% MARCADOR 1: PICO MÁXIMO
plot(t_peak, val_max, 'r*', 'MarkerSize', 10, 'LineWidth', 2);
text(t_peak, val_max*1.02, sprintf('Pico: %.2f°', val_max), 'Color', 'r');

% MARCADOR 2: TIEMPO DE ESTABLECIMIENTO
% Línea vertical donde se estabiliza
plot([ts ts], [0 y_final_deg], 'm-.', 'LineWidth', 1.5);
text(ts, y_final_deg*0.5, sprintf(' Ts = %.3fs', ts), 'Color', 'm', 'FontWeight', 'bold');

hold off;
grid on;
title(['Evaluación Temporal: Respuesta de Lazo Cerrado (Ref: 12°)']);
xlabel('Tiempo (s)'); ylabel('Apertura (Grados)');
legend('Salida', 'Referencia', 'Banda 2%', 'Banda 2%', 'Pico Máx', 'Tiempo Estab.');

% --- MOSTRAR RESULTADOS TRACKING ---
fprintf('\n--------------------------------------------------\n');
fprintf('RESULTADOS DEL DISEÑO (Simulación Final):\n');
fprintf('--------------------------------------------------\n');
fprintf('Valor Final: %.4f grados\n', y_final_deg);
fprintf('Sobrepaso (Overshoot): %.2f %%\n', Overshoot_pct);
fprintf('Tiempo de Establecimiento (2%%): %.4f s\n', ts);

% ----------------------------------------------------------
% B. ANÁLISIS DE ESTABILIDAD ABSOLUTA
% ----------------------------------------------------------
polos_LC = pole(T_final);
partes_reales = real(polos_LC);

fprintf('\n--------------------------------------------------\n');
fprintf('ANÁLISIS DE ESTABILIDAD ABSOLUTA\n');
fprintf('--------------------------------------------------\n');
fprintf('Polos de Lazo Cerrado obtenidos:\n');
disp(polos_LC);

if all(partes_reales < 0)
    fprintf('CONCLUSIÓN: El sistema es ESTABLE.\n');
    fprintf('(Todos los polos tienen parte real negativa).\n');
else
    fprintf('CONCLUSIÓN: El sistema es INESTABLE.\n');
end
fprintf('--------------------------------------------------\n');


% 1. Definimos la Función de Lazo Abierto L(s) = C(s) * Gp(s)
L_open_loop = Controlador_Final * Gp;
% Analizamos el error en estado estable luego del controlador
minreal(L_open_loop)

% ==========================================================
% --- PERTURBACION: COMPARATIVA SIN vs CON CONTROL ---
% ==========================================================

% Parámetro de Carga Realista (Torque reflejado en el motor)
Torque_Viento = 0.05; % Nm (50 mNm)

% 1. Definición de la Perturbación (Física del sistema)
num_dist = [-L, -R];
G_dist_OL = tf(num_dist, den_plant);

% 2. Cálculo del Lazo Cerrado (Con tu controlador)
L_loop = Controlador_Final * Gp;
S = feedback(1, L_loop);
T_perturbacion_CL = G_dist_OL * S;

% --- FIGURA 1: ESCENARIO SIN CONTROLADOR (LAZO ABIERTO) ---
figure(1);
t_short = 0:0.001:0.5; % Solo 0.5s porque se va al infinito
[y_ol, t_ol] = step(G_dist_OL, t_short);

% Escalamos por el torque
y_ol = y_ol * Torque_Viento;

plot(t_ol, y_ol*180/pi, 'LineWidth', 2);
grid on;
title(['Efecto del Viento SIN Control (Torque = ' num2str(Torque_Viento) ' Nm)']);
xlabel('Tiempo (s)'); ylabel('Desviación (Grados)');
legend('Caída libre del flap (Inestable)');

% --- FIGURA 2: ESCENARIO CON CONTROLADOR (LAZO CERRADO) ---
figure(2);
t_long = 0:0.01:2; % 2 segundos para ver la estabilidad
[y_cl, t_cl] = step(T_perturbacion_CL, t_long);

% Escalamos por el torque y pasamos a grados
y_cl = y_cl * Torque_Viento * (180/pi);

plot(t_cl, y_cl, 'b-', 'LineWidth', 2);
grid on;
hold on;

% Dibujamos una línea roja punteada en y = -12 desde t=0 hasta el final
plot([t_long(1) t_long(end)], [-12 -12], 'r--', 'LineWidth', 2);

hold off;

title(['Efecto del Viento CON Control (Torque = ' num2str(Torque_Viento) ' Nm)']);
xlabel('Tiempo (s)'); ylabel('Desviación (Grados)');
legend('Compensador sosteniendo la carga', 'Límite Mecánico (-12°)');

% --- ANÁLISIS AUTOMÁTICO ---
val_final_grados = y_cl(end);

fprintf('\n--- ANÁLISIS COMPARATIVO DE ROBUSTEZ ---\n');
fprintf('Carga de Viento aplicada: %.3f Nm\n', Torque_Viento);
fprintf('1. SIN CONTROL: El sistema es inestable (la gráfica 1 se va al infinito).\n');
fprintf('2. CON CONTROL: Desviación estable de %.2f grados.\n', val_final_grados);

if abs(val_final_grados) < 12
    fprintf('CONCLUSIÓN: El controlador evita el colapso del mecanismo.\n');
else
    fprintf('CONCLUSIÓN: Aún con control, la carga es muy alta.\n');
end
