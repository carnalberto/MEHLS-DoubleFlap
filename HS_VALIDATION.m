%% HS VALIDATION
clear
close all
clc

%% HSPM2DMP 
% Validation of HSPM2DMP with x-foil for a single airfoil:
xy_0012 = fopen('naca0012.dat');
data_0012 = textscan(xy_0012, '%f %f %f', 'HeaderLines', 3, ...
    'CollectOutput', 1, ...
    'Delimiter','');
fclose(xy_0012);

xpos_MAIN= data_0012{1,1}(:,1);
ypos_MAIN = data_0012{1,1}(:,2);

% modifica xfoil:
A_INPUT = [xpos_MAIN ypos_MAIN];
A_INPUT = flipud(A_INPUT);

V_INF = 67;
ALPHA=deg2rad(2);
GAP=[0.0;0.0];
OUT_SINGLE=HSPM2DMP(V_INF,14,GAP,ALPHA,A_INPUT);

Cp_XFOIL = 'prova0012.dat';
data_Cp = fopen(Cp_XFOIL);
data_tot = textscan(data_Cp, '%f %f %f', 'HeaderLines', 3, ...
    'CollectOutput', 1, ...
    'Delimiter','');
fclose(data_Cp);


% Separate CP data:
x0 = data_tot{1,1}(:,1);
y0 = data_tot{1,1}(:,2);
Cp = data_tot{1,1}(:,3);

figure(1); grid on; hold on;
plot(x0, -Cp, OUT_SINGLE.Cp{1}(:,1), -OUT_SINGLE.Cp{1}(:,2), '--', Linewidth=2); 
xlabel('x [m]', Interpreter='LaTeX', FontSize=14);
ylabel('$-C_p$ [-]', Interpreter='LaTeX', FontSize=14);
legend('Xfoil', '\verb|HSPM2DMP.m|', FontSize=18, Location='North', Interpreter='latex');
set(gca, TickLabelInterpreter='LaTeX', FontSize=18);


%% PM3
% Validation of PM3 trough HSPM2DMP for 3 airfoils
xy_632A015 = fopen('naca63(2)A-015.dat');
data_632A015 = textscan(xy_632A015, '%f %f %f', 'HeaderLines', 3, ...
    'CollectOutput', 1, ...
    'Delimiter','');
fclose(xy_632A015);

xpos_MAIN = data_632A015{1,1}(:,1);
ypos_MAIN = data_632A015{1,1}(:,2);

A_INPUT = [xpos_MAIN ypos_MAIN];
A_INPUT = flipud(A_INPUT);

V_INF = 67;
ALPHA=deg2rad(0);

xy_0012 = fopen('naca0012.dat');
data_0012 = textscan(xy_0012, '%f %f %f', 'HeaderLines', 3, ...
    'CollectOutput', 1, ...
    'Delimiter','');
fclose(xy_0012);

xpos_FIRST = data_0012{1,1}(:,1);
ypos_FIRST = data_0012{1,1}(:,2);

B_INPUT = [xpos_FIRST ypos_FIRST];
B_INPUT = flipud(B_INPUT);

C_INPUT = B_INPUT*0.3;

% ANGLE SET 1
ALPHA=deg2rad(0);
DELTA=deg2rad(-32);
THETA=deg2rad(-60);

GAP=[0.0;-0.05;0.0;-0.025];
OUT_PM3=PM3(ALPHA,DELTA,THETA,A_INPUT,B_INPUT,C_INPUT,GAP,V_INF);

figure(2);
grid on; hold on;
color = get(gca,'ColorOrder');

plot(OUT_PM3.Cp_A(:,1),-OUT_PM3.Cp_A(:,2), Linewidth=2, Color=color(1,:));
plot(OUT_PM3.Cp_B(:,1),-OUT_PM3.Cp_B(:,2), Linewidth=2, Color=color(1,:));
plot(OUT_PM3.Cp_C(:,1),-OUT_PM3.Cp_C(:,2), Linewidth=2, Color=color(1,:));

GAP=[0.0,  0.0 ,   0.0;
     0.0,-0.05 ,-0.025];

ALPHA=-ALPHA;
BETA=DELTA+ALPHA;
GAMMA=THETA+ALPHA;
OUT_HSPM2DMP=HSPM2DMP(V_INF,14,GAP,ALPHA,A_INPUT,BETA,B_INPUT,GAMMA,C_INPUT);

Cl_HSPM2DMP=OUT_HSPM2DMP.Cl{1} + OUT_HSPM2DMP.Cl{2}*1/3 + OUT_HSPM2DMP.Cl{3}*0.3/3;


figure(2);
plot(OUT_HSPM2DMP.Cp{1}(:,1), -OUT_HSPM2DMP.Cp{1}(:,2), '--', Linewidth=2, Color=color(2,:));
plot(OUT_HSPM2DMP.Cp{2}(:,1), -OUT_HSPM2DMP.Cp{2}(:,2), '--', Linewidth=2, Color=color(2,:)); 
plot(OUT_HSPM2DMP.Cp{3}(:,1), -OUT_HSPM2DMP.Cp{3}(:,2), '--', Linewidth=2, Color=color(2,:));

xlabel('x [m]', Interpreter='LaTeX', FontSize=14);
ylabel('$-C_p$ [-]', Interpreter='LaTeX', FontSize=14);
legend('\verb|PM3.m|', '', '', '\verb|HSPM2DMP.m|', '', '', FontSize=18, Location='North', Interpreter='latex');
set(gca, TickLabelInterpreter='LaTeX', FontSize=18);


%% WILLIAMS TEST CASE VALIDATION:

% Main_coord = readmatrix('MainFoil_N=300.csv');
Main_coord = readmatrix('MainFoil_N=200.csv');
% Main_coord = readmatrix('MainFoil_N=100.csv');
% Main_coord = readmatrix('MainFoil_N=50.csv');
x_Main = Main_coord(:,1);
y_Main = Main_coord(:,2);
x_Main = x_Main(2:end-1);
y_Main = y_Main(2:end-1);

% Flap_coord = readmatrix('FlapFoil_N=300.csv');
Flap_coord = readmatrix('FlapFoil_N=200.csv');
% Flap_coord = readmatrix('FlapFoil_N=100.csv');
% Flap_coord = readmatrix('FlapFoil_N=50.csv');
x_Flap = Flap_coord(:,1);
y_Flap = Flap_coord(:,2);
x_Flap = x_Flap(2:end-1);
y_Flap = y_Flap(2:end-1);

Main = readmatrix('Cp_Main_theoretical.csv');
x_MainCp = Main(:,1);
Cp_Main = Main(:,2);

Flap = readmatrix('Cp_Flap_theoretical.csv');
x_FlapCp = Flap(:,1);
Cp_Flap = Flap(:,2);


%% LEs coordinates:

% Flap's L.E. and gaps:
[LE_Flap_x, LE_Flap_I] = min(x_Flap);
LE_Flap_y = y_Flap(LE_Flap_I);
Flap_chord = sqrt((LE_Flap_y - y_Flap(1))^2 + (LE_Flap_x - x_Flap(1))^2);

% Flap's coords. translated in zero:
x_Flap0 = x_Flap - abs(LE_Flap_x);
y_Flap0 = y_Flap + abs(LE_Flap_y);

x_Gap = abs(LE_Flap_x - x_Main(1));
y_Gap = abs(LE_Flap_y - y_Main(1));


% HSPM2DMP
V_INF = 67;

A_INPUT = [x_Main y_Main];
A_INPUT = flipud(A_INPUT);

B_INPUT = [x_Flap0 y_Flap0];
B_INPUT = flipud(B_INPUT);

% ANGLE SET 1
ALPHA = deg2rad(0);
DELTA = deg2rad(0);

GAP = [0.0, -x_Gap; 
       0.0, -y_Gap];

ALPHA = -ALPHA;
BETA = DELTA + ALPHA;

OUT_HSPM2DMP = HSPM2DMP(V_INF, 14, GAP, ALPHA, A_INPUT, BETA, B_INPUT);
Cl_HSPM2DMP = OUT_HSPM2DMP.Cl{1} + OUT_HSPM2DMP.Cl{2} * Flap_chord / 1;


% GRAPHS & OUTPUTS:
fprintf('\n \n Total C_L (Williams) is: [-]: %.5f \n \n', Cl_HSPM2DMP);

% Coordinates
h3 = figure(3);
subplot(121);
grid on; hold on; axis equal;
plot(x_Main, y_Main, Linewidth=2);
plot(x_Flap, y_Flap, Linewidth=2);
xlabel('x [m]', Interpreter='LaTeX', FontSize=14);
ylabel('y [m]', Interpreter='LaTeX', FontSize=14);
legend('Main', 'Flap', FontSize=18, Location='NorthEast', Interpreter='latex');
% title('\textbf{WINK''s experiments (low Re)}', FontSize=14, Interpreter='latex');
set(gca, TickLabelInterpreter='LaTeX', FontSize=18);

% Graphs:
subplot(122); grid on; hold on;
color = get(gca,'ColorOrder');

plot(x_MainCp, -Cp_Main, Linewidth=2, Color=color(1,:));
plot(x_FlapCp, -Cp_Flap, Linewidth=2, Color=color(1,:));

plot(OUT_HSPM2DMP.Cp{1}(:,1), -OUT_HSPM2DMP.Cp{1}(:,2), '--', Linewidth=2, Color=color(2,:));
plot(OUT_HSPM2DMP.Cp{2}(:,1), -OUT_HSPM2DMP.Cp{2}(:,2), '--', Linewidth=2, Color=color(2,:)); 

xlabel('x [m]', Interpreter='LaTeX', FontSize=14);
ylabel('$-C_p$ [-]', Interpreter='LaTeX', FontSize=14);

legend('Williams', '', '\verb|HSPM2DMP.m|', '', FontSize=18, Location='North', Interpreter='latex');
set(gca, TickLabelInterpreter='LaTeX', FontSize=18);








