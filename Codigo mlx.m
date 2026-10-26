Codigos de Proyecto 
Sección: 4
N° de grupo: 8
Titulo del Proyecto: Simulación numérica y control de una suspensión activa de cuarto de vehículo
Integrante 1 (líder): Diego Alonso Carpio Panduro (100% de participación)          
Integrante 2: Sayuri Kiara Lima Pucllas (100% de participación)
Integrante 3: Victor Alonso Tipe Morillo (100% de participación)
Integrante 4: Pedro Franklin Salazar Isique (100% de participación)
Integrante 5: Ricardo André Nobuyoshi Obregón Espinoza (100% de participación)
Observación: 
No olvidar añadir el porcentaje de participación.
Incluir comentarios relevantes para entender el codigo, tener en cuenta que comentarios de IA anularan la entrega
El codigo debe de tener una identacion adecuada para facilitar su lectura
El codigo debe der replicable en cualquier computadora (Uso adecuado de semillas para funciones)
Las funciones principales deben de ir en la parte inferior de este archivo
Titulos de los codigos deben reflejar su funcionalidad
Codigo 1: Parámetros, modelo y simulación con Euler mejorado, RK4 y Dormand-Prince (ode45)
% =========================================================================
% CÓDIGO 1: Parámetros del Sistema, Modelo y Ejecución de Métodos
% =========================================================================

clear; clc; close all;

% 1. Parámetros del Sistema
ms = 290;       % Masa suspendida (kg)
mu = 59;        % Masa no suspendida (kg)
ks = 16812;     % Rigidez de suspensión (N/m)
cs = 1000;      % Amortiguamiento de suspensión (N.s/m)
kt = 190000;    % Rigidez del neumático (N/m)
ct = 0;         % Amortiguamiento del neumático (N.s/m)

% 2. Parámetros del Resalto y Velocidad
a = 0.05;                % Altura del resalto (m)
L = 3.50;                % Longitud del resalto (m)
v = 25 * (1000 / 3600);  % Velocidad (m/s)

t1 = 0.50;          % Inicio del contacto (s)
t2 = t1 + (L / v);  % Salida del resalto (s)

% 3. Constantes del Controlador PID
Kp = 1000;      % Constante proporcional
Ki = 100;       % Constante integral
Kd = 50;        % Constante derivativa

% 4. Perfil de la Superficie zr(t) y su Derivada dzr(t)
get_zr  = @(t) (t >= t1 && t <= t2) * (a / 2) * (1 - cos((2 * pi * v * (t - t1)) / L)); 
get_dzr = @(t) (t >= t1 && t <= t2) * (a / 2) * ((2 * pi * v) / L) * sin((2 * pi * v * (t - t1)) / L);

% 5. Ecuaciones de Estado del Sistema
dot_x = @(t, state, f_act) [
    state(2);
    (-ks/ms)*state(1) - (cs/ms)*state(2) + (ks/ms)*state(3) + (cs/ms)*state(4) + (1/ms)*f_act;
    state(4);
    (ks/mu)*state(1) + (cs/mu)*state(2) - ((ks+kt)/mu)*state(3) - ((cs+ct)/mu)*state(4) - (1/mu)*f_act + (kt/mu)*get_zr(t) + (ct/mu)*get_dzr(t)
    ];

% 6. Configuración de Tiempo
t_ini = 0; t_fin = 4.0; dt = 0.001;
t_vec = t_ini:dt:t_fin;

% 7. Ejecución de los Métodos Numéricos mediante Funciones Principales
[t_e, x_e, f_e] = metodo_euler_mejorado(dot_x, t_vec, Kp, Ki, Kd);
[t_r, x_r, f_r] = metodo_runge_kutta_4(dot_x, t_vec, Kp, Ki, Kd);
[t_p, x_p, f_p] = metodo_dormand_prince(ms, mu, ks, cs, kt, ct, Kp, Ki, Kd, get_zr, get_dzr, [t_ini t_fin]);
Codigo 2: Gráficas de la respuesta con cada método (desplazamiento, deflexión y fuerza)
% ======================================
% CÓDIGO 2: Visualización de Gráficas
% ======================================

% ------------------------------------------
% 1. Gráfica - Euler Mejorado
% ------------------------------------------
zs_e = x_e(1, :);
zu_e = x_e(3, :);
suspension_deflection_e = zs_e - zu_e;

figure('Name', 'Resultados - Euler Mejorado', 'Position', [100, 100, 600, 750]);

subplot(3,1,1);
plot(t_e, zs_e, 'b-', 'LineWidth', 1.5);
ylabel('z_s (m)');
title('Euler Mejorado');
grid on;

subplot(3,1,2);
plot(t_e, suspension_deflection_e, 'r-', 'LineWidth', 1.5);
ylabel('Deflexión Suspensión (m)');
grid on;

subplot(3,1,3);
plot(t_e, f_e, 'g-', 'LineWidth', 1.5);
xlabel('Tiempo (s)');
ylabel('Fuerza Activa f(t) (N)');
grid on;

% ----------------------------------------
% 2. Gráfica - Runge Kutta de 4to Orden
% ----------------------------------------
zs_r = x_r(1, :);
zu_r = x_r(3, :);
suspension_deflection_r = zs_r - zu_r;

figure('Name', 'Resultados - Runge Kutta 4to Orden', 'Position', [100, 100, 600, 750]);

subplot(3,1,1);
plot(t_r, zs_r, 'b-', 'LineWidth', 1.5);
ylabel('z_s (m)');
title('Runge Kutta 4to Orden');
grid on;

subplot(3,1,2);
plot(t_r, suspension_deflection_r, 'r-', 'LineWidth', 1.5);
ylabel('Deflexión Suspensión (m)');
grid on;

subplot(3,1,3);
plot(t_r, f_r, 'g-', 'LineWidth', 1.5);
xlabel('Tiempo (s)');
ylabel('Fuerza Activa f(t) (N)');
grid on;

% ----------------------------------------
% 3. Gráfica - Runge Kutta Dormand Prince
% ----------------------------------------
zs_p = x_p(1, :);
zu_p = x_p(3, :);
suspension_deflection_p = zs_p - zu_p;

figure('Name', 'Resultados - Runge Kutta Dormand Prince', 'Position', [100, 100, 600, 750]);

subplot(3,1,1);
plot(t_p, zs_p, 'b-', 'LineWidth', 1.5);
ylabel('z_s (m)');
title('Runge Kutta Dormand Prince');
grid on;

subplot(3,1,2);
plot(t_p, suspension_deflection_p, 'r-', 'LineWidth', 1.5);
ylabel('Deflexión Suspensión (m)');
grid on;

subplot(3,1,3);
plot(t_p, f_p, 'g-', 'LineWidth', 1.5);
xlabel('Tiempo (s)');
ylabel('Fuerza Activa f(t) (N)');
grid on;


Codigo 3: Análisis paramétrico de k_s, c_s, a y v (suspensión pasiva vs activa)
%Código 3: Análisis paramétrico
factores = [0.8 0.9 1.0 1.1 1.2];
nombres = {'ks', 'cs', 'a', 'v'};
base = [ks, cs, a, v];
filas = [];
for p = 1:4
    for j = 1:numel(factores)
        %aca se cambia solo el parametro p, los demas siguen con su valor base
        par = base;
        par(p) = base(p) * factores(j);
        %la primera llamada es sin control y la segunda con PID
        [acc_pas, defs_pas, deft_pas] = simular_caso(par, 0, 0, 0, ms, mu, kt, ct, L, t1, t_vec);
        [acc_act, defs_act, deft_act] = simular_caso(par, Kp, Ki, Kd, ms, mu, kt, ct, L, t1, t_vec);
        filas = [filas; p, factores(j), acc_pas, defs_pas, deft_pas, acc_act, defs_act, deft_act];
    end
end
Resultados_param = table(nombres(filas(:,1))', filas(:,2), filas(:,3), filas(:,4), filas(:,5), filas(:,6), filas(:,7), filas(:,8), ...
    'VariableNames', {'Parametro', 'Factor', 'AccMax_pas', 'DefSuspMax_pas', 'DefNeumMax_pas', 'AccMax_act', 'DefSuspMax_act', 'DefNeumMax_act'});
disp(Resultados_param);
%esto para que la deflexión de la suspensión no pase de 8cm
limite_ok = all(filas(:, [4 7]) <= 0.08, 'all');
disp(limite_ok); fprintf('Número de corridas: %d casos (%d simulaciones: pasiva y activa)\n', size(filas,1), 2*size(filas,1)); % corridas totales
Codigo 4: Comparación de error y tiempo de los métodos; figuras y tablas APA
% ================================================================
% CÓDIGO 4: Figuras y Tablas de Resultados (Parte 4, formato APA)
% ================================================================
% Usa los parámetros, el perfil del resalto y dot_x definidos en el Código 1.

% 2. Configuración de la presentación (APA 7)
dir_fig = 'Figuras_APA';           % Carpeta de salida (JPEG y tablas)
if ~exist(dir_fig, 'dir'), mkdir(dir_fig); end
fuente = 'Times New Roman';
dpi    = 300;
gris   = [0.45 0.45 0.45];

% 3. Simulaciones
t_ini = 0; t_fin = 4.0;
dt     = 0.001;                    % Paso de los métodos de paso fijo (s)
dt_ref = 0.0001;                   % Paso de la solución de referencia (s)

% Solución de referencia: RK4 con paso muy pequeño
t_ref = t_ini:dt_ref:t_fin;
x_ref = simular('rk4', dot_x, t_ref, dt_ref, Kp, Ki, Kd);

% Métodos de paso fijo
t_c = t_ini:dt:t_fin;

tic;
[x_em, f_em] = simular('em', dot_x, t_c, dt, Kp, Ki, Kd);
tiempo_em = toc;

tic;
[x_rk, f_rk] = simular('rk4', dot_x, t_c, dt, Kp, Ki, Kd);
tiempo_rk = toc;

% Dormand-Prince (ode45), con el integrador del error como 5to estado
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
f_pid = @(s) Kp*(0 - s(1)) + Ki*s(5) + Kd*(-s(2));
sistema = @(time, s) [dot_x(time, s(1:4), f_pid(s)); 0 - s(1)];

tic;
[t_dp, x_dp] = ode45(sistema, [t_ini t_fin], zeros(5,1), options);
tiempo_dp = toc;

f_dp = Kp*(-x_dp(:,1)) + Ki*x_dp(:,5) + Kd*(-x_dp(:,2));

% 4. Errores respecto a la referencia
e_em = abs(x_em(1,:) - interp1(t_ref, x_ref(1,:), t_c));
e_rk = abs(x_rk(1,:) - interp1(t_ref, x_ref(1,:), t_c));
x_dp_i = interp1(t_dp, x_dp(:,1), t_ref, 'spline');
e_dp = abs(x_dp_i - x_ref(1,:));

err_max = [max(e_em); max(e_rk); max(e_dp)];
tiempos = [tiempo_em; tiempo_rk; tiempo_dp];

% 5. Tablas
Metodo = {'Euler mejorado'; 'Runge-Kutta de cuarto orden'; 'Dormand-Prince (ode45)'};
Tabla1 = table(Metodo, err_max, tiempos, ...
    'VariableNames', {'Metodo', 'ErrorMaximo_m', 'TiempoEjecucion_s'});

m_em = metricas(t_c, x_em(1,:), x_em(3,:), f_em);
m_rk = metricas(t_c, x_rk(1,:), x_rk(3,:), f_rk);
m_dp = metricas(t_dp', x_dp(:,1)', x_dp(:,3)', f_dp');
M = [m_em; m_rk; m_dp];
Tabla2 = table(Metodo, M(:,1), M(:,2), M(:,3), M(:,4), M(:,5), M(:,6), ...
    'VariableNames', {'Metodo', 'PicoZs_m', 'TiempoPico_s', 'ValleZs_m', 'DeflexionMax_m', 'FuerzaMax_N', 'Zs_t4s_m'});

disp(Tabla1); disp(Tabla2);
writetable(Tabla1, fullfile(dir_fig, 'Tabla_1_Comparacion.csv'), 'Encoding', 'UTF-8');
writetable(Tabla2, fullfile(dir_fig, 'Tabla_2_Respuesta.csv'), 'Encoding', 'UTF-8');

% 6. Figuras (APA: sin título interno, sin cuadrícula, escala de grises)

% Figura 1. Perfil del resalto
zr = arrayfun(get_zr, t_c);
[f, ax] = nueva_figura(fuente, 3.0);
plot(ax, t_c, zr, 'k-', 'LineWidth', 1.6);
cerrar_figura(f, ax, 'Tiempo (s)', 'Elevación del resalto, \itz\rm_{r} (m)', t1, t2, '', fullfile(dir_fig, 'Figura_1_Perfil_resalto.jpg'), dpi);

% Figura 2. Desplazamiento de la carrocería
[f, ax] = nueva_figura(fuente, 3.6);
plot(ax, t_c, x_em(1,:), 'k-', 'LineWidth', 1.6); hold(ax, 'on');
plot(ax, t_c, x_rk(1,:), '--', 'Color', gris, 'LineWidth', 1.6);
plot(ax, t_dp, x_dp(:,1), 'k:', 'LineWidth', 2);
cerrar_figura(f, ax, 'Tiempo (s)', 'Desplazamiento de la carrocería, \itz\rm_{s} (m)', t1, t2, 'northeast', fullfile(dir_fig, 'Figura_2_Desplazamiento_carroceria.jpg'), dpi);

% Figura 3. Deflexión de la suspensión
[f, ax] = nueva_figura(fuente, 3.6);
plot(ax, t_c, x_em(1,:) - x_em(3,:), 'k-', 'LineWidth', 1.6); hold(ax, 'on');
plot(ax, t_c, x_rk(1,:) - x_rk(3,:), '--', 'Color', gris, 'LineWidth', 1.6);
plot(ax, t_dp, x_dp(:,1) - x_dp(:,3), 'k:', 'LineWidth', 2);
cerrar_figura(f, ax, 'Tiempo (s)', 'Deflexión de la suspensión, \itz\rm_{s} - \itz\rm_{u} (m)', t1, t2, 'northeast', fullfile(dir_fig, 'Figura_3_Deflexion_suspension.jpg'), dpi);

% Figura 4. Fuerza del actuador
[f, ax] = nueva_figura(fuente, 3.6);
plot(ax, t_c, f_em, 'k-', 'LineWidth', 1.6); hold(ax, 'on');
plot(ax, t_c, f_rk, '--', 'Color', gris, 'LineWidth', 1.6);
plot(ax, t_dp, f_dp, 'k:', 'LineWidth', 2);
cerrar_figura(f, ax, 'Tiempo (s)', 'Fuerza del actuador, \itf\rm(\itt\rm) (N)', t1, t2, 'northeast', fullfile(dir_fig, 'Figura_4_Fuerza_activa.jpg'), dpi);

% Figura 5. Error absoluto respecto a la referencia (escala logarítmica)
[f, ax] = nueva_figura(fuente, 3.6);
semilogy(ax, t_c, max(e_em, 1e-12), 'k-', 'LineWidth', 1.6); hold(ax, 'on');
semilogy(ax, t_c, max(e_rk, 1e-12), '--', 'Color', gris, 'LineWidth', 1.6);
semilogy(ax, t_ref, max(e_dp, 1e-12), 'k:', 'LineWidth', 2);
ylim(ax, [1e-9 1e-3]);
cerrar_figura(f, ax, 'Tiempo (s)', 'Error absoluto en \itz\rm_{s} (m)', t1, t2, 'northeast', fullfile(dir_fig, 'Figura_5_Error_absoluto.jpg'), dpi);

% Figura 6. Error máximo y tiempo de ejecución
f = figure('Units', 'inches', 'Position', [1 1 6.5 3.2], 'Color', 'w');
etiquetas = {'Euler mejorado', 'RK4', 'Dormand-Prince'};
colores = [0.25 0.25 0.25; 0.6 0.6 0.6; 0.9 0.9 0.9];

ax1 = subplot(1, 2, 1);
b = bar(ax1, 1:3, err_max, 'FaceColor', 'flat', 'EdgeColor', 'k');
b.CData = colores;
set(ax1, 'YScale', 'log', 'XTick', 1:3, 'XTickLabel', etiquetas, 'FontName', fuente, 'FontSize', 9, 'Box', 'off', 'TickDir', 'out');
ylim(ax1, [1e-7 1e-4]); ylabel(ax1, 'Error máximo (m)');
for k = 1:3, text(ax1, k, err_max(k)*1.2, sprintf('%.2e', err_max(k)), 'HorizontalAlignment', 'center', 'FontName', fuente, 'FontSize', 8); end

ax2 = subplot(1, 2, 2);
b = bar(ax2, 1:3, tiempos, 'FaceColor', 'flat', 'EdgeColor', 'k');
b.CData = colores;
set(ax2, 'XTick', 1:3, 'XTickLabel', etiquetas, 'FontName', fuente, 'FontSize', 9, 'Box', 'off', 'TickDir', 'out');
ylim(ax2, [0 1.2*max(tiempos)]); ylabel(ax2, 'Tiempo de ejecución (s)');
for k = 1:3, text(ax2, k, tiempos(k) + 0.02*max(tiempos), sprintf('%.4f', tiempos(k)), 'HorizontalAlignment', 'center', 'FontName', fuente, 'FontSize', 8); end
exportgraphics(f, fullfile(dir_fig, 'Figura_6_Error_y_tiempo.jpg'), 'Resolution', dpi);

fprintf('Figuras y tablas guardadas en la carpeta "%s".\n', dir_fig);
Funciones Principales
% ======================
% FUNCIONES PRINCIPALES 
% ======================

% 1. Función para Euler Mejorado
function [t_vec, x_out, f_act_vec] = metodo_euler_mejorado(dot_x, t_vec, Kp, Ki, Kd)
N = length(t_vec);
dt = t_vec(2) - t_vec(1);
x_out = zeros(4, N);
f_act_vec = zeros(1, N);
int_e = 0;

for i = 1:N-1
    time = t_vec(i);
    state = x_out(:, i);
    error = 0 - state(1);
    if i > 1, int_e = int_e + error * dt; end
    f_act = Kp * error + Ki * int_e + Kd * (-state(2));
    f_act_vec(i) = f_act;

    k1 = dot_x(time, state, f_act);
    k2 = dot_x(t_vec(i+1), state + dt*k1, f_act);
    x_out(:, i+1) = state + (dt / 2) * (k1 + k2);
end
f_act_vec(N) = f_act_vec(N-1);
end

% 2. Función para Runge Kutta 4to Orden
function [t_vec, x_out, f_act_vec] = metodo_runge_kutta_4(dot_x, t_vec, Kp, Ki, Kd)
N = length(t_vec);
dt = t_vec(2) - t_vec(1);
x_out = zeros(4, N);
f_act_vec = zeros(1, N);
int_r = 0;

for i = 1:N-1
    time = t_vec(i);
    state = x_out(:, i);
    error = 0 - state(1);
    if i > 1, int_r = int_r + error * dt; end
    f_act = Kp * error + Ki * int_r + Kd * (-state(2));
    f_act_vec(i) = f_act;

    k1 = dot_x(time, state, f_act);
    k2 = dot_x(time + dt/2, state + (dt/2)*k1, f_act);
    k3 = dot_x(time + dt/2, state + (dt/2)*k2, f_act);
    k4 = dot_x(time + dt, state + dt*k3, f_act);
    x_out(:, i+1) = state + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4);
end
f_act_vec(N) = f_act_vec(N-1);
end

% 3. Función para Runge Kutta Dormand Prince 
function [t_adap, x_adap, f_act_vec] = metodo_dormand_prince(ms, mu, ks, cs, kt, ct, Kp, Ki, Kd, get_zr, get_dzr, tspan)
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
[t_adap, x_raw] = ode45(@(time, s) system_ode_adap_aux(time, s, ms, mu, ks, cs, kt, ct, Kp, Ki, Kd, get_zr, get_dzr), tspan, [0;0;0;0;0], options);

x_adap = x_raw(:, 1:4)';
int_err = x_raw(:, 5)';
f_act_vec = Kp * (0 - x_adap(1, :)) + Ki * int_err + Kd * (-x_adap(2, :));
end

% 4. Función Auxiliar para Dormand Prince
function ds = system_ode_adap_aux(time, s, ms, mu, ks, cs, kt, ct, Kp, Ki, Kd, get_zr, get_dzr)
z_s = s(1); z_dot_s = s(2); z_u = s(3); z_dot_u = s(4); int_error = s(5);
error = 0 - z_s; der_error = -z_dot_s;
f_act = Kp * error + Ki * int_error + Kd * der_error;

dzs = z_dot_s;
dzs_dot = (-ks/ms)*z_s - (cs/ms)*z_dot_s + (ks/ms)*z_u + (cs/ms)*z_dot_u + (1/ms)*f_act;
dzu = z_dot_u;
dzu_dot = (ks/mu)*z_s + (cs/mu)*z_dot_s - ((ks+kt)/mu)*z_u - ((cs+ct)/mu)*z_dot_u - (1/mu)*f_act + (kt/mu)*get_zr(time) + (ct/mu)*get_dzr(time);
d_interror = error;
ds = [dzs; dzs_dot; dzu; dzu_dot; d_interror];
end

% 5. Función para un caso del análisis paramétrico
%esto devuelve los máximos de aceleración, deflexión de la suspensión y del neumático
function [acc_max, defs_max, deft_max] = simular_caso(par, Kp, Ki, Kd, ms, mu, kt, ct, L, t1, t_vec)
ks_ = par(1); cs_ = par(2); a_ = par(3); v_ = par(4);
t2_ = t1 + L / v_;
zr_f  = @(t) (t >= t1 && t <= t2_) * (a_ / 2) * (1 - cos((2 * pi * v_ * (t - t1)) / L));
dzr_f = @(t) (t >= t1 && t <= t2_) * (a_ / 2) * ((2 * pi * v_) / L) * sin((2 * pi * v_ * (t - t1)) / L);
%estas son las mismas ecuaciones del código1 pero con los parámetros de este caso
dx = @(t, s, f_act) [
    s(2);
    (-ks_/ms)*s(1) - (cs_/ms)*s(2) + (ks_/ms)*s(3) + (cs_/ms)*s(4) + (1/ms)*f_act;
    s(4);
    (ks_/mu)*s(1) + (cs_/mu)*s(2) - ((ks_+kt)/mu)*s(3) - ((cs_+ct)/mu)*s(4) - (1/mu)*f_act + (kt/mu)*zr_f(t) + (ct/mu)*dzr_f(t)
    ];
[~, x, f] = metodo_runge_kutta_4(dx, t_vec, Kp, Ki, Kd);
N = numel(t_vec);
acc = zeros(1, N);
zr_vec = zeros(1, N);
for i = 1:N
    d = dx(t_vec(i), x(:, i), f(i));
    acc(i) = d(2);
    zr_vec(i) = zr_f(t_vec(i));
end
acc_max  = max(abs(acc));
defs_max = max(abs(x(1, :) - x(3, :)));
deft_max = max(abs(x(3, :) - zr_vec));
end

% 6. Funciones de la Parte 4 (Resultados)
function [x, f_control] = simular(metodo, dot_x, t, dt, Kp, Ki, Kd)
% Integra el sistema con PID. metodo: 'em' (Euler mejorado) o 'rk4' (Runge-Kutta 4to orden)
    N = length(t);
    x = zeros(4, N);
    f_control = zeros(1, N);
    integral_error = 0;
    for i = 1:N-1
        time = t(i);
        state = x(:, i);
        error = 0 - state(1);
        if i > 1
            integral_error = integral_error + error * dt;
        end
        f_act = Kp * error + Ki * integral_error + Kd * (-state(2));
        f_control(i) = f_act;
        switch metodo
            case 'em'
                k1 = dot_x(time, state, f_act);
                k2 = dot_x(t(i+1), state + dt * k1, f_act);
                x(:, i+1) = state + (dt / 2) * (k1 + k2);
            case 'rk4'
                k1 = dot_x(time, state, f_act);
                k2 = dot_x(time + dt/2, state + (dt/2) * k1, f_act);
                k3 = dot_x(time + dt/2, state + (dt/2) * k2, f_act);
                k4 = dot_x(time + dt, state + dt * k3, f_act);
                x(:, i+1) = state + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4);
        end
    end
    f_control(N) = f_control(N-1);
end

function m = metricas(t, zs, zu, f)
% [pico zs, tiempo del pico, valle zs, deflexión máx., fuerza máx., zs final]
    [pico, i] = max(zs);
    m = [pico, t(i), min(zs), max(abs(zs - zu)), max(abs(f)), zs(end)];
end

function [f, ax] = nueva_figura(fuente, alto)
    f = figure('Units', 'inches', 'Position', [1 1 6.5 alto], 'Color', 'w');
    ax = axes(f);
    set(ax, 'FontName', fuente, 'FontSize', 11, 'Box', 'off', 'LineWidth', 0.8, 'TickDir', 'out');
end

function cerrar_figura(f, ax, xl, yl, t1, t2, leyenda, archivo, dpi)
    hold(ax, 'on');
    xline(ax, t1, '-.', 'Color', [0.7 0.7 0.7], 'HandleVisibility', 'off');
    xline(ax, t2, '-.', 'Color', [0.7 0.7 0.7], 'HandleVisibility', 'off');
    xlim(ax, [0 4]);
    xlabel(ax, xl); ylabel(ax, yl);
    if ~isempty(leyenda)
        legend(ax, {'Euler mejorado', 'Runge-Kutta de cuarto orden', 'Dormand-Prince (ode45)'}, 'Location', leyenda, 'Box', 'off', 'FontSize', 9);
    end
    exportgraphics(f, archivo, 'Resolution', dpi);
end