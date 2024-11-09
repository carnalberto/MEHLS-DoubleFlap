%% OPTIMIZATION
clear
close all
clc

tic;

% Ranges of variation of the gaps
GAPS=struct('B_H', [ -0.1  -0.05  0.0  0.05  0.1], 'B_V',  -[0.05  0.1  0.15  0.2  0.25], ...
            'C_H', [-0.05  -0.025  0.0  0.025  0.05], 'C_V', -[0.025  0.05  0.075  0.1  0.125]);
       
RANGE=length(GAPS.B_H);

% Allocating memory 
CELL=cell(RANGE^4,1);
ANALYSIS=struct('Cp_A',CELL,'Cp_B', CELL, 'Cp_C', CELL, 'Cl', CELL,...
                'GAP', zeros(4,1), 'STALL', zeros(3,1));
%% AIRFOILS
xy_0012 = fopen('naca0012.dat');
data_0012 = textscan(xy_0012, '%f %f %f', 'HeaderLines', 3, ...
    'CollectOutput', 1, ...
    'Delimiter','');
fclose(xy_0012);

xy_632A015 = fopen('naca63(2)A-015.dat');
data_632A015 = textscan(xy_632A015, '%f %f %f', 'HeaderLines', 3, ...
    'CollectOutput', 1, ...
    'Delimiter','');
fclose(xy_632A015);

xpos_MAIN = data_632A015{1,1}(:,1);
ypos_MAIN = data_632A015{1,1}(:,2);

xpos_FIRST = data_0012{1,1}(:,1);
ypos_FIRST = data_0012{1,1}(:,2);

% modifica xfoil:
A_INPUT = [xpos_MAIN ypos_MAIN];
A_INPUT = flipud(A_INPUT);

B_INPUT = [xpos_FIRST ypos_FIRST];
B_INPUT = flipud(B_INPUT);

C_INPUT = 0.3*B_INPUT;
  
% Aerodynamic parameters
% Reynolds Number
V_INF = 67;

%% SIMULATIONS
% % ANGLE SET 1
% ALPHA=deg2rad(0.0);
% DELTA=deg2rad(-32);
% THETA=deg2rad(-60);

% % ANGLE SET 2
% ALPHA=deg2rad(0.0);
% DELTA=deg2rad(-20);
% THETA=deg2rad(-40);

% ANGLE SET 3
ALPHA=deg2rad(0.0);
DELTA=deg2rad(-15);
THETA=deg2rad(-30);

STALL_MATRIX=zeros(RANGE^4,3);
MK=0;
I_P=0;

POINTER=[];

att = waitbar(0,'Wait, numeric integration in progress...');

for a=1:RANGE
    for b=1:RANGE
        for c=1:RANGE
            for d=1:RANGE
                MK=MK+1;
                GAP=[GAPS.B_H(a); GAPS.B_V(b); GAPS.C_H(c); GAPS.C_V(d)];
                OUTPUT=PM3(ALPHA,DELTA,THETA,A_INPUT,B_INPUT,C_INPUT,GAP,V_INF);
                ANALYSIS(MK).Cp_A=OUTPUT.Cp_A;
                ANALYSIS(MK).Cp_B=OUTPUT.Cp_B;
                ANALYSIS(MK).Cp_C=OUTPUT.Cp_C;
                ANALYSIS(MK).Cl=OUTPUT.Cl;
                ANALYSIS(MK).GAP=GAP;
                ANALYSIS(MK).STALL=OUTPUT.STALL;
                STALL_MATRIX(MK,:)=OUTPUT.STALL;
                STALL=ANALYSIS(MK).STALL;

                if STALL==[0,0,0]
                    I_P=I_P+1;
                    POINTER(I_P)=MK;              
                end

                waitbar(d/RANGE);
            end
        end
    end
end

close (att);


%% EXTRACTION OF THE NON STALLED CONFIGURATIONS
I_P_MAX=length(POINTER);
Cl=zeros(I_P_MAX,1);
CELL=cell(I_P_MAX,1);
FINAL_DRAFT=struct('Cp_A',CELL,'Cp_B', CELL, 'Cp_C', CELL, 'Cl', CELL,...
                'GAP', zeros(4,1), 'STALL', zeros(3,1));
            for I_P=1:length(POINTER)
                FINAL_DRAFT(I_P).Cp_A=ANALYSIS(POINTER(I_P)).Cp_A;
                FINAL_DRAFT(I_P).Cp_B=ANALYSIS(POINTER(I_P)).Cp_B;
                FINAL_DRAFT(I_P).Cp_C=ANALYSIS(POINTER(I_P)).Cp_C;
                FINAL_DRAFT(I_P).Cl=ANALYSIS(POINTER(I_P)).Cl;
                FINAL_DRAFT(I_P).STALL=ANALYSIS(POINTER(I_P)).STALL;
                FINAL_DRAFT(I_P).GAP=ANALYSIS(POINTER(I_P)).GAP;
                Cl(I_P)=FINAL_DRAFT(I_P).Cl;
            end
            
[Cl_MAX,I_Cl_MAX]=max(Cl);
GAP_Cl_MAX=FINAL_DRAFT(I_Cl_MAX).GAP;


% Plot optimum Cp
Cp_A=FINAL_DRAFT(I_Cl_MAX).Cp_A;
Cp_B=FINAL_DRAFT(I_Cl_MAX).Cp_B;
Cp_C=FINAL_DRAFT(I_Cl_MAX).Cp_C;

figure(2);
%subplot(224);
plot(Cp_A(:,1),-Cp_A(:,2), Cp_B(:,1),-Cp_B(:,2), Cp_C(:,1),-Cp_C(:,2), 'Linewidth',2); grid on;
xlabel('x [m]', 'FontSize',18, 'Interpreter','latex'); 
ylabel('$-C_p$', 'FontSize',18, 'Interpreter','latex');
title('\textbf{$C_p$}', 'FontSize',14, 'Interpreter','latex');
legend('Main', 'First', 'Second', 'FontSize',18, 'Location','North', 'Interpreter','latex');

fprintf('\n Time spent for this simulation [sec]: %.5f \n \n', toc);

fprintf('\n \n Max Cl for this configuration [-]: %.5f', Cl_MAX);
fprintf('\n \n Check Cl [-]: %.5f \n \n', ANALYSIS(POINTER(I_Cl_MAX)).Cl);
disp(' At gap')
disp(GAP_Cl_MAX);


fprintf('\n \n Number of configurations: %5d', RANGE^4);
fprintf('\n \n Number of cases without stall: %5d \n \n', I_P_MAX);










