function[OUTPUT]=PM3(ALPHA,DELTA,THETA,A_INPUT,B_INPUT,C_INPUT,GAP,V_INF)

% DESCRIPTION:  provides the potential, incompressible, flow solution (Cp,
%               Vt) on the surface of each profile. 
%               The profiles' surfaces are discretized with panels.
%               The flow solution is given at each panel's midpoint.
%               It takes as input the arrays of the coordinates of the
%               nodes of each profile with respect to a
%               reference frame centered in the leading edge of each
%               profile and aligned with the chord of each profile.
%               The sources and vortexes strengths per unit length 
%               distributions are uniform on each panel.             
%               The sources strength per unit length varies from panel to
%               panel and and for different profiles.
%               The vortexes strength per unit length is the same for 
%               all panels of a specific profile but different from profile
%               to profile. 
%               This is the 'Hess-Smith' choice for panel method.
%               It solves a linear system to find all the sources and
%               vortexes strengths per unit length.
%               It implemets the Valarezo-Chin criterion (DELTA_Cp < 14
%               between LE and TE) to see if the flow separates on one or
%               more profiles.

%% INPUT 

% % Aerodynamic parameters
% % Reynolds Number
% V_INF = 67;
% 
% % Specify the numbers of profiles n
% n=3;

%% PARTE DA MODIFICARE DURANTE I CICLI DI OTTIMIZZAZIONE 
% (non cambiare i nomi delle variabili)

% Gaps between the profiles
GAP_B_H=GAP(1);
GAP_B_V=GAP(2);
GAP_C_H=GAP(3);
GAP_C_V=GAP(4);

% Data: arrays each containing the coordinates of the nodes of a profile
% referred to the leading edge and at zero AoA

% xy_0012 = fopen('naca0012.dat');
% data_0012 = textscan(xy_0012, '%f %f %f', 'HeaderLines', 3, ...
%     'CollectOutput', 1, ...
%     'Delimiter','');
% fclose(xy_0012);
% 
% xy_632A015 = fopen('naca63(2)A-015.dat');
% data_632A015 = textscan(xy_632A015, '%f %f %f', 'HeaderLines', 3, ...
%     'CollectOutput', 1, ...
%     'Delimiter','');
% fclose(xy_632A015);
% 
% xpos_MAIN = data_632A015{1,1}(:,1);
% ypos_MAIN = data_632A015{1,1}(:,2);
% 
% xpos_FIRST = data_0012{1,1}(:,1);
% ypos_FIRST = data_0012{1,1}(:,2);
% 
% % modifica xfoil:
% A_INPUT = [xpos_MAIN ypos_MAIN];
% A_INPUT = flipud(A_INPUT);
% 
% B_INPUT = [xpos_FIRST ypos_FIRST];
% B_INPUT = flipud(B_INPUT);
% 
% C_INPUT = B_INPUT/3;
% FINE PARTE DA MODIFICARE

LENGTH_A=length(A_INPUT);
LENGTH_B=length(B_INPUT);
LENGTH_C=length(C_INPUT);

chord_A=A_INPUT(1,1);
chord_B=B_INPUT(1,1);
chord_C=C_INPUT(1,1);

TRIG_ALPHA=[cos(ALPHA),sin(ALPHA)];

% Deflections rotation matrixes
R_DELTA=[cos(DELTA),-sin(DELTA);
         sin(DELTA), cos(DELTA)];

R_THETA=[cos(THETA),-sin(THETA);
         sin(THETA), cos(THETA)];

% Angles with respect to flow velocity (ccw+)
% BETA=ALPHA-DELTA;
% GAMMA=ALPHA-THETA;

%% COMPUTATION OF GEOMETRY PARAMETERS
% Leading edges positions with respect to the absolute reference frame
% (pole= main's LE, axis: x-> main's chord to the right, y-> vertical)
B_LE=[chord_A+GAP_B_H , GAP_B_V];

C_LE = B_LE+[GAP_C_H,GAP_C_V]+(R_DELTA*[chord_B; 0.0])';

% Nodes coordinate rotated 
B_ROT=zeros(LENGTH_B,2);

C_ROT=zeros(LENGTH_C,2);

 for i=1:LENGTH_B
     LOCAL=B_INPUT(i,:);
     B_ROT(i,:)=(R_DELTA*LOCAL')';
 end
 
 for i=1:LENGTH_C
     LOCAL=C_INPUT(i,:);
     C_ROT(i,:)=(R_THETA*LOCAL')';
 end
 
% Nodes coordinates in absolute refererence frame (translation)
A=A_INPUT; % Does not change

B=B_ROT+B_LE;

C=C_ROT+C_LE;

% Number of panels of each airfol
N_P_A=LENGTH_A-1;
N_P_B=LENGTH_B-1;
N_P_C=LENGTH_C-1;

% Panels' lengths
L_P_A=zeros(N_P_A,1);
for i=1:N_P_A
    L_P_A(i)=sqrt((A(i+1,1)-A(i,1))^2+(A(i+1,2)-A(i,2))^2); 
end

L_P_B=zeros(N_P_B,2);
for i=1:N_P_B
    L_P_B(i)=sqrt((B(i+1,1)-B(i,1))^2+(B(i+1,2)-B(i,2))^2); 
end
    
L_P_C=zeros(N_P_C,2);
for i=1:N_P_C
    L_P_C(i)=sqrt((C(i+1,1)-C(i,1))^2+(C(i+1,2)-C(i,2))^2); 
end

% Trigonometric functions of panels' angles 
TRIG_A=zeros(N_P_A,2);
for i=1:N_P_A
    TRIG_A(i,1)=(A(i+1,1)-A(i,1))/L_P_A(i); % cos
    TRIG_A(i,2)=(A(i+1,2)-A(i,2))/L_P_A(i); % sin
end

TRIG_B=zeros(N_P_B,2);
for i=1:N_P_B
    TRIG_B(i,1)=(B(i+1,1)-B(i,1))/L_P_B(i); % cos
    TRIG_B(i,2)=(B(i+1,2)-B(i,2))/L_P_B(i); % sin
    
end

TRIG_C=zeros(N_P_C,2);
for i=1:N_P_C
    TRIG_C(i,1)=(C(i+1,1)-C(i,1))/L_P_C(i); % cos
    TRIG_C(i,2)=(C(i+1,2)-C(i,2))/L_P_C(i); % sin
end

% Panels' midpoints
MID_A=zeros(N_P_A,2);
for i=1:N_P_A
    MID_A(i,1)=(A(i+1,1)+A(i,1))/2; 
    MID_A(i,2)=(A(i+1,2)+A(i,2))/2; 
end

MID_B=zeros(N_P_B,2);
for i=1:N_P_B
    MID_B(i,1)=(B(i+1,1)+B(i,1))/2; 
    MID_B(i,2)=(B(i+1,2)+B(i,2))/2; 
end

MID_C=zeros(N_P_C,2);
for i=1:N_P_C
    MID_C(i,1)=(C(i+1,1)+C(i,1))/2; 
    MID_C(i,2)=(C(i+1,2)+C(i,2))/2; 
end

%% Nodes-midpoints angles and distances
% K and L are the profiles indexes A, B, C.
% i and j are the panels' indexes.
% i-->K
% j-->L

% B_KL(i,j) angle formed by the segments connecting node j on profile L to
% midpoint i on profile K and node j+1 on profile L

% R_KL(i,j) distance of midpoint i on profile K from node j on profile L. 

% B_AA & R_AA 
N_P_K=N_P_A; % Number of panels on profile K --> rows
N_P_L=N_P_A; % Number of panels on profile L --> columns
MID_K=MID_A;
L=A;
TRIG_L=TRIG_A;
% Not to be modified vvv
B_KL=zeros(N_P_K,N_P_L);
R_KL=zeros(N_P_K,N_P_L+1);
for i=1:N_P_K
    for j=1:N_P_L
        % Relative positions midpoint i wrt node j and wrt node j+1  
        % in local components (panel's j reference frame)
        P1_ABSOLUTE=[MID_K(i,1)-L(j,1); MID_K(i,2)-L(j,2)];
        P2_ABSOLUTE=[MID_K(i,1)-L(j+1,1); MID_K(i,2)-L(j+1,2)];
        R_INV=[TRIG_L(j,1) , TRIG_L(j,2); -TRIG_L(j,2) , TRIG_L(j,1)];
        P1_LOCAL=R_INV*P1_ABSOLUTE;
        P2_LOCAL=R_INV*P2_ABSOLUTE;
        THETA1_LOCAL=atan2(P1_LOCAL(2),P1_LOCAL(1));
        THETA2_LOCAL=atan2(P2_LOCAL(2),P2_LOCAL(1));
        % Fix in case of autoinduction
        if (abs(THETA1_LOCAL)<10^(-12) && abs(THETA2_LOCAL)>3)
            THETA1_LOCAL=0; 
            THETA2_LOCAL=pi;
        end
        if (abs(THETA2_LOCAL)<10^(-12) && abs(THETA1_LOCAL)>3)
            THETA2_LOCAL=0; 
            THETA1_LOCAL=-pi; 
        end
        B_KL(i,j)= THETA2_LOCAL-THETA1_LOCAL;  
        R_KL(i,j)=sqrt( ( MID_K(i,1)-L(j,1) )^2 + ( MID_K(i,2)-L(j,2) )^2 );
    end
    R_KL(i,N_P_L+1)=sqrt( ( MID_K(i,1)-L(N_P_L+1,1) )^2 + ( MID_K(i,2)-L(N_P_L+1,2) )^2 );
end
% Not to be modified ^^^
B_AA=B_KL;
R_AA=R_KL;

% B_AB & R_AB
N_P_K=N_P_A;    % Number of panels on profile K --> rows
N_P_L=N_P_B;    % Number of panels on profile L --> columns
MID_K=MID_A;    
L=B;            
TRIG_L=TRIG_B;  
% Not to be modified vvv
B_KL=zeros(N_P_K,N_P_L);
R_KL=zeros(N_P_K,N_P_L+1);
for i=1:N_P_K
    for j=1:N_P_L
        % Relative positions midpoint i wrt node j and wrt node j+1  
        % in local components (panel's j reference frame)
        P1_ABSOLUTE=[MID_K(i,1)-L(j,1); MID_K(i,2)-L(j,2)];
        P2_ABSOLUTE=[MID_K(i,1)-L(j+1,1); MID_K(i,2)-L(j+1,2)];
        R_INV=[TRIG_L(j,1) , TRIG_L(j,2); -TRIG_L(j,2) , TRIG_L(j,1)];
        P1_LOCAL=R_INV*P1_ABSOLUTE;
        P2_LOCAL=R_INV*P2_ABSOLUTE;
        THETA1_LOCAL=atan2(P1_LOCAL(2),P1_LOCAL(1));
        THETA2_LOCAL=atan2(P2_LOCAL(2),P2_LOCAL(1));
        % Fix in case of autoinduction
        if (abs(THETA1_LOCAL)<10^(-12) && abs(THETA2_LOCAL)>3)
            THETA1_LOCAL=0; 
            THETA2_LOCAL=pi;
        end
        if (abs(THETA2_LOCAL)<10^(-12) && abs(THETA1_LOCAL)>3)
            THETA2_LOCAL=0; 
            THETA1_LOCAL=-pi; 
        end
        B_KL(i,j)= THETA2_LOCAL-THETA1_LOCAL;  
        R_KL(i,j)=sqrt( ( MID_K(i,1)-L(j,1) )^2 + ( MID_K(i,2)-L(j,2) )^2 );
    end
    R_KL(i,N_P_L+1)=sqrt( ( MID_K(i,1)-L(N_P_L+1,1) )^2 + ( MID_K(i,2)-L(N_P_L+1,2) )^2 );
end
% Not to be modified ^^^
B_AB=B_KL;
R_AB=R_KL;

% B_AC & R_AC
N_P_K=N_P_A; % Number of panels on profile K --> rows
N_P_L=N_P_C; % Number of panels on profile L --> columns
MID_K=MID_A;
L=C;
TRIG_L=TRIG_C;
% Not to be modified vvv
B_KL=zeros(N_P_K,N_P_L);
R_KL=zeros(N_P_K,N_P_L+1);
for i=1:N_P_K
    for j=1:N_P_L
        % Relative positions midpoint i wrt node j and wrt node j+1  
        % in local components (panel's j reference frame)
        P1_ABSOLUTE=[MID_K(i,1)-L(j,1); MID_K(i,2)-L(j,2)];
        P2_ABSOLUTE=[MID_K(i,1)-L(j+1,1); MID_K(i,2)-L(j+1,2)];
        R_INV=[TRIG_L(j,1) , TRIG_L(j,2); -TRIG_L(j,2) , TRIG_L(j,1)];
        P1_LOCAL=R_INV*P1_ABSOLUTE;
        P2_LOCAL=R_INV*P2_ABSOLUTE;
        THETA1_LOCAL=atan2(P1_LOCAL(2),P1_LOCAL(1));
        THETA2_LOCAL=atan2(P2_LOCAL(2),P2_LOCAL(1));
        % Fix in case of autoinduction
        if (abs(THETA1_LOCAL)<10^(-12) && abs(THETA2_LOCAL)>3)
            THETA1_LOCAL=0; 
            THETA2_LOCAL=pi;
        end
        if (abs(THETA2_LOCAL)<10^(-12) && abs(THETA1_LOCAL)>3)
            THETA2_LOCAL=0; 
            THETA1_LOCAL=-pi; 
        end
        B_KL(i,j)= THETA2_LOCAL-THETA1_LOCAL;  
        R_KL(i,j)=sqrt( ( MID_K(i,1)-L(j,1) )^2 + ( MID_K(i,2)-L(j,2) )^2 );
    end
    R_KL(i,N_P_L+1)=sqrt( ( MID_K(i,1)-L(N_P_L+1,1) )^2 + ( MID_K(i,2)-L(N_P_L+1,2) )^2 );
end
% Not to be modified ^^^
B_AC=B_KL;
R_AC=R_KL;

% B_BA & R_BA
N_P_K=N_P_B; % Number of panels on profile K --> rows
N_P_L=N_P_A; % Number of panels on profile L --> columns
MID_K=MID_B;
L=A;
TRIG_L=TRIG_A;
% Not to be modified vvv
B_KL=zeros(N_P_K,N_P_L);
R_KL=zeros(N_P_K,N_P_L+1);
for i=1:N_P_K
    for j=1:N_P_L
        % Relative positions midpoint i wrt node j and wrt node j+1  
        % in local components (panel's j reference frame)
        P1_ABSOLUTE=[MID_K(i,1)-L(j,1); MID_K(i,2)-L(j,2)];
        P2_ABSOLUTE=[MID_K(i,1)-L(j+1,1); MID_K(i,2)-L(j+1,2)];
        R_INV=[TRIG_L(j,1) , TRIG_L(j,2); -TRIG_L(j,2) , TRIG_L(j,1)];
        P1_LOCAL=R_INV*P1_ABSOLUTE;
        P2_LOCAL=R_INV*P2_ABSOLUTE;
        THETA1_LOCAL=atan2(P1_LOCAL(2),P1_LOCAL(1));
        THETA2_LOCAL=atan2(P2_LOCAL(2),P2_LOCAL(1));
        % Fix in case of autoinduction
        if (abs(THETA1_LOCAL)<10^(-12) && abs(THETA2_LOCAL)>3)
            THETA1_LOCAL=0; 
            THETA2_LOCAL=pi;
        end
        if (abs(THETA2_LOCAL)<10^(-12) && abs(THETA1_LOCAL)>3)
            THETA2_LOCAL=0; 
            THETA1_LOCAL=-pi; 
        end
        B_KL(i,j)= THETA2_LOCAL-THETA1_LOCAL;  
        R_KL(i,j)=sqrt( ( MID_K(i,1)-L(j,1) )^2 + ( MID_K(i,2)-L(j,2) )^2 );
    end
    R_KL(i,N_P_L+1)=sqrt( ( MID_K(i,1)-L(N_P_L+1,1) )^2 + ( MID_K(i,2)-L(N_P_L+1,2) )^2 );
end
% Not to be modified ^^^
B_BA=B_KL;
R_BA=R_KL;

% B_BB & R_BB
N_P_K=N_P_B; % Number of panels on profile K --> rows
N_P_L=N_P_B; % Number of panels on profile L --> columns
MID_K=MID_B;
L=B;
TRIG_L=TRIG_B;
% Not to be modified vvv
B_KL=zeros(N_P_K,N_P_L);
R_KL=zeros(N_P_K,N_P_L+1);
for i=1:N_P_K
    for j=1:N_P_L
        % Relative positions midpoint i wrt node j and wrt node j+1  
        % in local components (panel's j reference frame)
        P1_ABSOLUTE=[MID_K(i,1)-L(j,1); MID_K(i,2)-L(j,2)];
        P2_ABSOLUTE=[MID_K(i,1)-L(j+1,1); MID_K(i,2)-L(j+1,2)];
        R_INV=[TRIG_L(j,1) , TRIG_L(j,2); -TRIG_L(j,2) , TRIG_L(j,1)];
        P1_LOCAL=R_INV*P1_ABSOLUTE;
        P2_LOCAL=R_INV*P2_ABSOLUTE;
        THETA1_LOCAL=atan2(P1_LOCAL(2),P1_LOCAL(1));
        THETA2_LOCAL=atan2(P2_LOCAL(2),P2_LOCAL(1));
        % Fix in case of autoinduction
        if (abs(THETA1_LOCAL)<10^(-12) && abs(THETA2_LOCAL)>3)
            THETA1_LOCAL=0; 
            THETA2_LOCAL=pi;
        end
        if (abs(THETA2_LOCAL)<10^(-12) && abs(THETA1_LOCAL)>3)
            THETA2_LOCAL=0; 
            THETA1_LOCAL=-pi; 
        end
        B_KL(i,j)= THETA2_LOCAL-THETA1_LOCAL;  
        R_KL(i,j)=sqrt( ( MID_K(i,1)-L(j,1) )^2 + ( MID_K(i,2)-L(j,2) )^2 );
    end
    R_KL(i,N_P_L+1)=sqrt( ( MID_K(i,1)-L(N_P_L+1,1) )^2 + ( MID_K(i,2)-L(N_P_L+1,2) )^2 );
end
% Not to be modified ^^^
B_BB=B_KL;
R_BB=R_KL;

% B_BC & R_BC
N_P_K=N_P_B; % Number of panels on profile K --> rows
N_P_L=N_P_C; % Number of panels on profile L --> columns
MID_K=MID_B;
L=C;
TRIG_L=TRIG_C;
% Not to be modified vvv
B_KL=zeros(N_P_K,N_P_L);
R_KL=zeros(N_P_K,N_P_L+1);
for i=1:N_P_K
    for j=1:N_P_L
        % Relative positions midpoint i wrt node j and wrt node j+1  
        % in local components (panel's j reference frame)
        P1_ABSOLUTE=[MID_K(i,1)-L(j,1); MID_K(i,2)-L(j,2)];
        P2_ABSOLUTE=[MID_K(i,1)-L(j+1,1); MID_K(i,2)-L(j+1,2)];
        R_INV=[TRIG_L(j,1) , TRIG_L(j,2); -TRIG_L(j,2) , TRIG_L(j,1)];
        P1_LOCAL=R_INV*P1_ABSOLUTE;
        P2_LOCAL=R_INV*P2_ABSOLUTE;
        THETA1_LOCAL=atan2(P1_LOCAL(2),P1_LOCAL(1));
        THETA2_LOCAL=atan2(P2_LOCAL(2),P2_LOCAL(1));
        % Fix in case of autoinduction
        if (abs(THETA1_LOCAL)<10^(-12) && abs(THETA2_LOCAL)>3)
            THETA1_LOCAL=0; 
            THETA2_LOCAL=pi;
        end
        if (abs(THETA2_LOCAL)<10^(-12) && abs(THETA1_LOCAL)>3)
            THETA2_LOCAL=0; 
            THETA1_LOCAL=-pi; 
        end
        B_KL(i,j)= THETA2_LOCAL-THETA1_LOCAL;  
        R_KL(i,j)=sqrt( ( MID_K(i,1)-L(j,1) )^2 + ( MID_K(i,2)-L(j,2) )^2 );
    end
    R_KL(i,N_P_L+1)=sqrt( ( MID_K(i,1)-L(N_P_L+1,1) )^2 + ( MID_K(i,2)-L(N_P_L+1,2) )^2 );
end
% Not to be modified ^^^
B_BC=B_KL;
R_BC=R_KL;

% B_CA & R_CA
N_P_K=N_P_C; % Number of panels on profile K --> rows
N_P_L=N_P_A; % Number of panels on profile L --> columns
MID_K=MID_C;
L=A;
TRIG_L=TRIG_A;
% Not to be modified vvv
B_KL=zeros(N_P_K,N_P_L);
R_KL=zeros(N_P_K,N_P_L+1);
for i=1:N_P_K
    for j=1:N_P_L
        % Relative positions midpoint i wrt node j and wrt node j+1  
        % in local components (panel's j reference frame)
        P1_ABSOLUTE=[MID_K(i,1)-L(j,1); MID_K(i,2)-L(j,2)];
        P2_ABSOLUTE=[MID_K(i,1)-L(j+1,1); MID_K(i,2)-L(j+1,2)];
        R_INV=[TRIG_L(j,1) , TRIG_L(j,2); -TRIG_L(j,2) , TRIG_L(j,1)];
        P1_LOCAL=R_INV*P1_ABSOLUTE;
        P2_LOCAL=R_INV*P2_ABSOLUTE;
        THETA1_LOCAL=atan2(P1_LOCAL(2),P1_LOCAL(1));
        THETA2_LOCAL=atan2(P2_LOCAL(2),P2_LOCAL(1));
        % Fix in case of autoinduction
        if (abs(THETA1_LOCAL)<10^(-12) && abs(THETA2_LOCAL)>3)
            THETA1_LOCAL=0; 
            THETA2_LOCAL=pi;
        end
        if (abs(THETA2_LOCAL)<10^(-12) && abs(THETA1_LOCAL)>3)
            THETA2_LOCAL=0; 
            THETA1_LOCAL=-pi; 
        end
        B_KL(i,j)= THETA2_LOCAL-THETA1_LOCAL;  
        R_KL(i,j)=sqrt( ( MID_K(i,1)-L(j,1) )^2 + ( MID_K(i,2)-L(j,2) )^2 );
    end
    R_KL(i,N_P_L+1)=sqrt( ( MID_K(i,1)-L(N_P_L+1,1) )^2 + ( MID_K(i,2)-L(N_P_L+1,2) )^2 );
end
% Not to be modified ^^^
B_CA=B_KL;
R_CA=R_KL;

% B_CB & R_CB
N_P_K=N_P_C; % Number of panels on profile K --> rows
N_P_L=N_P_B; % Number of panels on profile L --> columns
MID_K=MID_C;
L=B;
TRIG_L=TRIG_B;
% Not to be modified vvv
B_KL=zeros(N_P_K,N_P_L);
R_KL=zeros(N_P_K,N_P_L+1);
for i=1:N_P_K
    for j=1:N_P_L
        % Relative positions midpoint i wrt node j and wrt node j+1  
        % in local components (panel's j reference frame)
        P1_ABSOLUTE=[MID_K(i,1)-L(j,1); MID_K(i,2)-L(j,2)];
        P2_ABSOLUTE=[MID_K(i,1)-L(j+1,1); MID_K(i,2)-L(j+1,2)];
        R_INV=[TRIG_L(j,1) , TRIG_L(j,2); -TRIG_L(j,2) , TRIG_L(j,1)];
        P1_LOCAL=R_INV*P1_ABSOLUTE;
        P2_LOCAL=R_INV*P2_ABSOLUTE;
        THETA1_LOCAL=atan2(P1_LOCAL(2),P1_LOCAL(1));
        THETA2_LOCAL=atan2(P2_LOCAL(2),P2_LOCAL(1));
        % Fix in case of autoinduction
        if (abs(THETA1_LOCAL)<10^(-12) && abs(THETA2_LOCAL)>3)
            THETA1_LOCAL=0; 
            THETA2_LOCAL=pi;
        end
        if (abs(THETA2_LOCAL)<10^(-12) && abs(THETA1_LOCAL)>3)
            THETA2_LOCAL=0; 
            THETA1_LOCAL=-pi; 
        end
        B_KL(i,j)= THETA2_LOCAL-THETA1_LOCAL;  
        R_KL(i,j)=sqrt( ( MID_K(i,1)-L(j,1) )^2 + ( MID_K(i,2)-L(j,2) )^2 );
    end
    R_KL(i,N_P_L+1)=sqrt( ( MID_K(i,1)-L(N_P_L+1,1) )^2 + ( MID_K(i,2)-L(N_P_L+1,2) )^2 );
end
% Not to be modified ^^^
B_CB=B_KL;
R_CB=R_KL;

% B_CC & R_CC
N_P_K=N_P_C; % Number of panels on profile K --> rows
N_P_L=N_P_C; % Number of panels on profile L --> columns
MID_K=MID_C;
L=C;
TRIG_L=TRIG_C;
% Not to be modified vvv
B_KL=zeros(N_P_K,N_P_L);
R_KL=zeros(N_P_K,N_P_L+1);
for i=1:N_P_K
    for j=1:N_P_L
        % Relative positions midpoint i wrt node j and wrt node j+1  
        % in local components (panel's j reference frame)
        P1_ABSOLUTE=[MID_K(i,1)-L(j,1); MID_K(i,2)-L(j,2)];
        P2_ABSOLUTE=[MID_K(i,1)-L(j+1,1); MID_K(i,2)-L(j+1,2)];
        R_INV=[TRIG_L(j,1) , TRIG_L(j,2); -TRIG_L(j,2) , TRIG_L(j,1)];
        P1_LOCAL=R_INV*P1_ABSOLUTE;
        P2_LOCAL=R_INV*P2_ABSOLUTE;
        THETA1_LOCAL=atan2(P1_LOCAL(2),P1_LOCAL(1));
        THETA2_LOCAL=atan2(P2_LOCAL(2),P2_LOCAL(1));
        % Fix in case of autoinduction
        if (abs(THETA1_LOCAL)<10^(-12) && abs(THETA2_LOCAL)>3)
            THETA1_LOCAL=0; 
            THETA2_LOCAL=pi;
        end
        if (abs(THETA2_LOCAL)<10^(-12) && abs(THETA1_LOCAL)>3)
            THETA2_LOCAL=0; 
            THETA1_LOCAL=-pi; 
        end
        B_KL(i,j)= THETA2_LOCAL-THETA1_LOCAL;  
        R_KL(i,j)=sqrt( ( MID_K(i,1)-L(j,1) )^2 + ( MID_K(i,2)-L(j,2) )^2 );
    end
    R_KL(i,N_P_L+1)=sqrt( ( MID_K(i,1)-L(N_P_L+1,1) )^2 + ( MID_K(i,2)-L(N_P_L+1,2) )^2 );
end
% Not to be modified ^^^
B_CC=B_KL;
R_CC=R_KL;

% % Plot geometry
% figure
% plot(A(:,1),A(:,2),'k','linewidth',2)
% axis equal
% hold
% plot(B(:,1),B(:,2),'k','linewidth',2)
% plot(C(:,1),C(:,2),'k','linewidth',2)
% grid on
% title('GEOMETRY CHECK','Color','k');

%% COMPUTATION OF SYSTEM'S COEFFICIENTS
% The sistem of linear equations is obtained by imposing the flow tangency
% condition at each panel's midpoint for each profile and the Kutta
% condition for each profile's trailing edge (first and last panels).
% 
% STRUCTURE OF THE SYSTEM: 
%            [ S_TAN , V_TAN ; S_KUT , V_KUT ]*[Q ; NI]==[B_TAN ; B_TAN];
%
% DIMENSION: (N_P_1 +...+ N_P_K +...+ N_P_n + n)^2                            

% S_TAN_AA & V_TAN_AA
N_P_K=N_P_A;
N_P_L=N_P_A;
TRIG_K=TRIG_A;
TRIG_L=TRIG_A;
R_KL=R_AA;
B_KL=B_AA;
% Not to be modified vvv
S_TAN_KL=zeros(N_P_K,N_P_L);
V_TAN_KL=zeros(N_P_K,1);
for i=1:N_P_K
    for j=1:N_P_L
        S_TAN_KL(i,j)=(TRIG_L(j,2)*TRIG_K(i,1)-TRIG_L(j,1)*TRIG_K(i,2))*...
                      (-1)/(2*pi)*log( R_KL(i,j+1)/R_KL(i,j) ) +...
                      (TRIG_L(j,2)*TRIG_K(i,2)+TRIG_L(j,1)*TRIG_K(i,1))*...
                      B_KL(i,j)/(2*pi);
        
        V_TAN_KL(i)=V_TAN_KL(i)+... 
                     (TRIG_L(j,2)*TRIG_K(i,1)-TRIG_L(j,1)*TRIG_K(i,2))*...
                     B_KL(i,j)/(2*pi) + ...
                     (TRIG_L(j,1)*TRIG_K(i,1)+TRIG_L(j,2)*TRIG_K(i,2))*...
                     1/(2*pi)*log(R_KL(i,j+1)/R_KL(i,j));
                     
    end
end
% Not to be modified ^^^
S_TAN_AA=S_TAN_KL;
V_TAN_AA=V_TAN_KL;

% S_TAN_AB & V_TAN_AB
N_P_K=N_P_A;
N_P_L=N_P_B;
TRIG_K=TRIG_A;
TRIG_L=TRIG_B;
R_KL=R_AB;
B_KL=B_AB;
% Not to be modified vvv
S_TAN_KL=zeros(N_P_K,N_P_L);
V_TAN_KL=zeros(N_P_K,1);
for i=1:N_P_K
    for j=1:N_P_L
        S_TAN_KL(i,j)=(TRIG_L(j,2)*TRIG_K(i,1)-TRIG_L(j,1)*TRIG_K(i,2))*...
                      (-1)/(2*pi)*log( R_KL(i,j+1)/R_KL(i,j) ) +...
                      (TRIG_L(j,2)*TRIG_K(i,2)+TRIG_L(j,1)*TRIG_K(i,1))*...
                      B_KL(i,j)/(2*pi);
        
        V_TAN_KL(i)=V_TAN_KL(i)+... 
                     (TRIG_L(j,2)*TRIG_K(i,1)-TRIG_L(j,1)*TRIG_K(i,2))*...
                     B_KL(i,j)/(2*pi) + ...
                     (TRIG_L(j,1)*TRIG_K(i,1)+TRIG_L(j,2)*TRIG_K(i,2))*...
                     1/(2*pi)*log(R_KL(i,j+1)/R_KL(i,j));
                     
    end
end
% Not to be modified ^^^
S_TAN_AB=S_TAN_KL;
V_TAN_AB=V_TAN_KL;

% S_TAN_AC & V_TAN_AC
N_P_K=N_P_A;
N_P_L=N_P_C;
TRIG_K=TRIG_A;
TRIG_L=TRIG_C;
R_KL=R_AC;
B_KL=B_AC;
% Not to be modified vvv
S_TAN_KL=zeros(N_P_K,N_P_L);
V_TAN_KL=zeros(N_P_K,1);
for i=1:N_P_K
    for j=1:N_P_L
        S_TAN_KL(i,j)=(TRIG_L(j,2)*TRIG_K(i,1)-TRIG_L(j,1)*TRIG_K(i,2))*...
                      (-1)/(2*pi)*log( R_KL(i,j+1)/R_KL(i,j) ) +...
                      (TRIG_L(j,2)*TRIG_K(i,2)+TRIG_L(j,1)*TRIG_K(i,1))*...
                      B_KL(i,j)/(2*pi);
        
        V_TAN_KL(i)=V_TAN_KL(i)+... 
                     (TRIG_L(j,2)*TRIG_K(i,1)-TRIG_L(j,1)*TRIG_K(i,2))*...
                     B_KL(i,j)/(2*pi) + ...
                     (TRIG_L(j,1)*TRIG_K(i,1)+TRIG_L(j,2)*TRIG_K(i,2))*...
                     1/(2*pi)*log(R_KL(i,j+1)/R_KL(i,j));
                     
    end
end
% Not to be modified ^^^
S_TAN_AC=S_TAN_KL;
V_TAN_AC=V_TAN_KL;

% S_TAN_BA & V_TAN_BA
N_P_K=N_P_B;
N_P_L=N_P_A;
TRIG_K=TRIG_B;
TRIG_L=TRIG_A;
R_KL=R_BA;
B_KL=B_BA;
% Not to be modified vvv
S_TAN_KL=zeros(N_P_K,N_P_L);
V_TAN_KL=zeros(N_P_K,1);
for i=1:N_P_K
    for j=1:N_P_L
        S_TAN_KL(i,j)=(TRIG_L(j,2)*TRIG_K(i,1)-TRIG_L(j,1)*TRIG_K(i,2))*...
                      (-1)/(2*pi)*log( R_KL(i,j+1)/R_KL(i,j) ) +...
                      (TRIG_L(j,2)*TRIG_K(i,2)+TRIG_L(j,1)*TRIG_K(i,1))*...
                      B_KL(i,j)/(2*pi);
        
        V_TAN_KL(i)=V_TAN_KL(i)+... 
                     (TRIG_L(j,2)*TRIG_K(i,1)-TRIG_L(j,1)*TRIG_K(i,2))*...
                     B_KL(i,j)/(2*pi) + ...
                     (TRIG_L(j,1)*TRIG_K(i,1)+TRIG_L(j,2)*TRIG_K(i,2))*...
                     1/(2*pi)*log(R_KL(i,j+1)/R_KL(i,j));
                     
    end
end
% Not to be modified ^^^
S_TAN_BA=S_TAN_KL;
V_TAN_BA=V_TAN_KL;

% S_TAN_BB & V_TAN_BB
N_P_K=N_P_B;
N_P_L=N_P_B;
TRIG_K=TRIG_B;
TRIG_L=TRIG_B;
R_KL=R_BB;
B_KL=B_BB;
% Not to be modified vvv
S_TAN_KL=zeros(N_P_K,N_P_L);
V_TAN_KL=zeros(N_P_K,1);
for i=1:N_P_K
    for j=1:N_P_L
        S_TAN_KL(i,j)=(TRIG_L(j,2)*TRIG_K(i,1)-TRIG_L(j,1)*TRIG_K(i,2))*...
                      (-1)/(2*pi)*log( R_KL(i,j+1)/R_KL(i,j) ) +...
                      (TRIG_L(j,2)*TRIG_K(i,2)+TRIG_L(j,1)*TRIG_K(i,1))*...
                      B_KL(i,j)/(2*pi);
        
        V_TAN_KL(i)=V_TAN_KL(i)+... 
                     (TRIG_L(j,2)*TRIG_K(i,1)-TRIG_L(j,1)*TRIG_K(i,2))*...
                     B_KL(i,j)/(2*pi) + ...
                     (TRIG_L(j,1)*TRIG_K(i,1)+TRIG_L(j,2)*TRIG_K(i,2))*...
                     1/(2*pi)*log(R_KL(i,j+1)/R_KL(i,j));
                     
    end
end
% Not to be modified ^^^
S_TAN_BB=S_TAN_KL;
V_TAN_BB=V_TAN_KL;

% S_TAN_BC & V_TAN_BC
N_P_K=N_P_B;
N_P_L=N_P_C;
TRIG_K=TRIG_B;
TRIG_L=TRIG_C;
R_KL=R_BC;
B_KL=B_BC;
% Not to be modified vvv
S_TAN_KL=zeros(N_P_K,N_P_L);
V_TAN_KL=zeros(N_P_K,1);
for i=1:N_P_K
    for j=1:N_P_L
        S_TAN_KL(i,j)=(TRIG_L(j,2)*TRIG_K(i,1)-TRIG_L(j,1)*TRIG_K(i,2))*...
                      (-1)/(2*pi)*log( R_KL(i,j+1)/R_KL(i,j) ) +...
                      (TRIG_L(j,2)*TRIG_K(i,2)+TRIG_L(j,1)*TRIG_K(i,1))*...
                      B_KL(i,j)/(2*pi);
        
        V_TAN_KL(i)=V_TAN_KL(i)+... 
                     (TRIG_L(j,2)*TRIG_K(i,1)-TRIG_L(j,1)*TRIG_K(i,2))*...
                     B_KL(i,j)/(2*pi) + ...
                     (TRIG_L(j,1)*TRIG_K(i,1)+TRIG_L(j,2)*TRIG_K(i,2))*...
                     1/(2*pi)*log(R_KL(i,j+1)/R_KL(i,j));
                     
    end
end
% Not to be modified ^^^
S_TAN_BC=S_TAN_KL;
V_TAN_BC=V_TAN_KL;

% S_TAN_CA & V_TAN_CA
N_P_K=N_P_C;
N_P_L=N_P_A;
TRIG_K=TRIG_C;
TRIG_L=TRIG_A;
R_KL=R_CA;
B_KL=B_CA;
% Not to be modified vvv
S_TAN_KL=zeros(N_P_K,N_P_L);
V_TAN_KL=zeros(N_P_K,1);
for i=1:N_P_K
    for j=1:N_P_L
        S_TAN_KL(i,j)=(TRIG_L(j,2)*TRIG_K(i,1)-TRIG_L(j,1)*TRIG_K(i,2))*...
                      (-1)/(2*pi)*log( R_KL(i,j+1)/R_KL(i,j) ) +...
                      (TRIG_L(j,2)*TRIG_K(i,2)+TRIG_L(j,1)*TRIG_K(i,1))*...
                      B_KL(i,j)/(2*pi);
        
        V_TAN_KL(i)=V_TAN_KL(i)+... 
                     (TRIG_L(j,2)*TRIG_K(i,1)-TRIG_L(j,1)*TRIG_K(i,2))*...
                     B_KL(i,j)/(2*pi) + ...
                     (TRIG_L(j,1)*TRIG_K(i,1)+TRIG_L(j,2)*TRIG_K(i,2))*...
                     1/(2*pi)*log(R_KL(i,j+1)/R_KL(i,j));
                     
    end
end
% Not to be modified ^^^
S_TAN_CA=S_TAN_KL;
V_TAN_CA=V_TAN_KL;

% S_TAN_CB & V_TAN_CB
N_P_K=N_P_C;
N_P_L=N_P_B;
TRIG_K=TRIG_C;
TRIG_L=TRIG_B;
R_KL=R_CB;
B_KL=B_CB;
% Not to be modified vvv
S_TAN_KL=zeros(N_P_K,N_P_L);
V_TAN_KL=zeros(N_P_K,1);
for i=1:N_P_K
    for j=1:N_P_L
        S_TAN_KL(i,j)=(TRIG_L(j,2)*TRIG_K(i,1)-TRIG_L(j,1)*TRIG_K(i,2))*...
                      (-1)/(2*pi)*log( R_KL(i,j+1)/R_KL(i,j) ) +...
                      (TRIG_L(j,2)*TRIG_K(i,2)+TRIG_L(j,1)*TRIG_K(i,1))*...
                      B_KL(i,j)/(2*pi);
        
        V_TAN_KL(i)=V_TAN_KL(i)+... 
                     (TRIG_L(j,2)*TRIG_K(i,1)-TRIG_L(j,1)*TRIG_K(i,2))*...
                     B_KL(i,j)/(2*pi) + ...
                     (TRIG_L(j,1)*TRIG_K(i,1)+TRIG_L(j,2)*TRIG_K(i,2))*...
                     1/(2*pi)*log(R_KL(i,j+1)/R_KL(i,j));
                     
    end
end
% Not to be modified ^^^
S_TAN_CB=S_TAN_KL;
V_TAN_CB=V_TAN_KL;

% S_TAN_CC & V_TAN_CC
N_P_K=N_P_C;
N_P_L=N_P_C;
TRIG_K=TRIG_C;
TRIG_L=TRIG_C;
R_KL=R_CC;
B_KL=B_CC;
% Not to be modified vvv
S_TAN_KL=zeros(N_P_K,N_P_L);
V_TAN_KL=zeros(N_P_K,1);
for i=1:N_P_K
    for j=1:N_P_L
        S_TAN_KL(i,j)=(TRIG_L(j,2)*TRIG_K(i,1)-TRIG_L(j,1)*TRIG_K(i,2))*...
                      (-1)/(2*pi)*log( R_KL(i,j+1)/R_KL(i,j) ) +...
                      (TRIG_L(j,2)*TRIG_K(i,2)+TRIG_L(j,1)*TRIG_K(i,1))*...
                      B_KL(i,j)/(2*pi);
        
        V_TAN_KL(i)=V_TAN_KL(i)+... 
                     (TRIG_L(j,2)*TRIG_K(i,1)-TRIG_L(j,1)*TRIG_K(i,2))*...
                     B_KL(i,j)/(2*pi) + ...
                     (TRIG_L(j,1)*TRIG_K(i,1)+TRIG_L(j,2)*TRIG_K(i,2))*...
                     1/(2*pi)*log(R_KL(i,j+1)/R_KL(i,j));
                     
    end
end
% Not to be modified ^^^
S_TAN_CC=S_TAN_KL;
V_TAN_CC=V_TAN_KL;

% S_KUT_AA & V_KUT_AA
N_P_K=N_P_A;
TRIG_K=TRIG_A;
TRIG_L=TRIG_A;
R_KL=R_AA;
B_KL=B_AA;
% Not to be modified vvv
S_KUT_KL=zeros(1,N_P_L);
V_KUT_KL=0.0;
for j=1:N_P_L
    S_KUT_KL(j)=(TRIG_L(j,1)*TRIG_K(1,1)+TRIG_L(j,2)*TRIG_K(1,2))*...
                (-1)/(2*pi)*log( R_KL(1,j+1)/R_KL(1,j) ) +...
                (-TRIG_L(j,2)*TRIG_K(1,1)+TRIG_L(j,1)*TRIG_K(1,2))*...
                B_KL(1,j)/(2*pi)+...
                (TRIG_L(j,1)*TRIG_K(N_P_K,1)+TRIG_L(j,2)*TRIG_K(N_P_K,2))*...
                (-1)/(2*pi)*log( R_KL(N_P_K,j+1)/R_KL(N_P_K,j) ) +...
                (-TRIG_L(j,2)*TRIG_K(N_P_K,1)+TRIG_L(j,1)*TRIG_K(N_P_K,2))*...
                B_KL(N_P_K,j)/(2*pi);
                
    V_KUT_KL=V_KUT_KL+(TRIG_L(j,1)*TRIG_K(1,1)+TRIG_L(j,2)*TRIG_K(1,2))*...
                B_KL(1,j)/(2*pi) +...
                (-TRIG_L(j,2)*TRIG_K(1,1)+TRIG_L(j,1)*TRIG_K(1,2))*...
                (1)/(2*pi)*log( R_KL(1,j+1)/R_KL(1,j) )+...
                (TRIG_L(j,1)*TRIG_K(N_P_K,1)+TRIG_L(j,2)*TRIG_K(N_P_K,2))*...
                B_KL(N_P_K,j)/(2*pi) +...
                (-TRIG_L(j,2)*TRIG_K(N_P_K,1)+TRIG_L(j,1)*TRIG_K(N_P_K,2))*...
                (1)/(2*pi)*log( R_KL(N_P_K,j+1)/R_KL(N_P_K,j) );
end 
% Not to be modified ^^^
S_KUT_AA=S_KUT_KL;
V_KUT_AA=V_KUT_KL;

% S_KUT_AB & V_KUT_AB
N_P_K=N_P_A;
TRIG_K=TRIG_A;
TRIG_L=TRIG_B;
R_KL=R_AB;
B_KL=B_AB;
% Not to be modified vvv
S_KUT_KL=zeros(1,N_P_L);
V_KUT_KL=0.0;
for j=1:N_P_L
    S_KUT_KL(j)=(TRIG_L(j,1)*TRIG_K(1,1)+TRIG_L(j,2)*TRIG_K(1,2))*...
                (-1)/(2*pi)*log( R_KL(1,j+1)/R_KL(1,j) ) +...
                (-TRIG_L(j,2)*TRIG_K(1,1)+TRIG_L(j,1)*TRIG_K(1,2))*...
                B_KL(1,j)/(2*pi)+...
                (TRIG_L(j,1)*TRIG_K(N_P_K,1)+TRIG_L(j,2)*TRIG_K(N_P_K,2))*...
                (-1)/(2*pi)*log( R_KL(N_P_K,j+1)/R_KL(N_P_K,j) ) +...
                (-TRIG_L(j,2)*TRIG_K(N_P_K,1)+TRIG_L(j,1)*TRIG_K(N_P_K,2))*...
                B_KL(N_P_K,j)/(2*pi);
                
    V_KUT_KL=V_KUT_KL+(TRIG_L(j,1)*TRIG_K(1,1)+TRIG_L(j,2)*TRIG_K(1,2))*...
                B_KL(1,j)/(2*pi) +...
                (-TRIG_L(j,2)*TRIG_K(1,1)+TRIG_L(j,1)*TRIG_K(1,2))*...
                (1)/(2*pi)*log( R_KL(1,j+1)/R_KL(1,j) )+...
                (TRIG_L(j,1)*TRIG_K(N_P_K,1)+TRIG_L(j,2)*TRIG_K(N_P_K,2))*...
                B_KL(N_P_K,j)/(2*pi) +...
                (-TRIG_L(j,2)*TRIG_K(N_P_K,1)+TRIG_L(j,1)*TRIG_K(N_P_K,2))*...
                (1)/(2*pi)*log( R_KL(N_P_K,j+1)/R_KL(N_P_K,j) );
end 
% Not to be modified ^^^
S_KUT_AB=S_KUT_KL;
V_KUT_AB=V_KUT_KL;

% S_KUT_AC & V_KUT_AC
N_P_K=N_P_A;
TRIG_K=TRIG_A;
TRIG_L=TRIG_C;
R_KL=R_AC;
B_KL=B_AC;
% Not to be modified vvv
S_KUT_KL=zeros(1,N_P_L);
V_KUT_KL=0.0;
for j=1:N_P_L
    S_KUT_KL(j)=(TRIG_L(j,1)*TRIG_K(1,1)+TRIG_L(j,2)*TRIG_K(1,2))*...
                (-1)/(2*pi)*log( R_KL(1,j+1)/R_KL(1,j) ) +...
                (-TRIG_L(j,2)*TRIG_K(1,1)+TRIG_L(j,1)*TRIG_K(1,2))*...
                B_KL(1,j)/(2*pi)+...
                (TRIG_L(j,1)*TRIG_K(N_P_K,1)+TRIG_L(j,2)*TRIG_K(N_P_K,2))*...
                (-1)/(2*pi)*log( R_KL(N_P_K,j+1)/R_KL(N_P_K,j) ) +...
                (-TRIG_L(j,2)*TRIG_K(N_P_K,1)+TRIG_L(j,1)*TRIG_K(N_P_K,2))*...
                B_KL(N_P_K,j)/(2*pi);
                
    V_KUT_KL=V_KUT_KL+(TRIG_L(j,1)*TRIG_K(1,1)+TRIG_L(j,2)*TRIG_K(1,2))*...
                B_KL(1,j)/(2*pi) +...
                (-TRIG_L(j,2)*TRIG_K(1,1)+TRIG_L(j,1)*TRIG_K(1,2))*...
                (1)/(2*pi)*log( R_KL(1,j+1)/R_KL(1,j) )+...
                (TRIG_L(j,1)*TRIG_K(N_P_K,1)+TRIG_L(j,2)*TRIG_K(N_P_K,2))*...
                B_KL(N_P_K,j)/(2*pi) +...
                (-TRIG_L(j,2)*TRIG_K(N_P_K,1)+TRIG_L(j,1)*TRIG_K(N_P_K,2))*...
                (1)/(2*pi)*log( R_KL(N_P_K,j+1)/R_KL(N_P_K,j) );
end 
% Not to be modified ^^^
S_KUT_AC=S_KUT_KL;
V_KUT_AC=V_KUT_KL;

% S_KUT_BA & V_KUT_BA
N_P_K=N_P_B;
TRIG_K=TRIG_B;
TRIG_L=TRIG_A;
R_KL=R_BA;
B_KL=B_BA;
% Not to be modified vvv
S_KUT_KL=zeros(1,N_P_L);
V_KUT_KL=0.0;
for j=1:N_P_L
    S_KUT_KL(j)=(TRIG_L(j,1)*TRIG_K(1,1)+TRIG_L(j,2)*TRIG_K(1,2))*...
                (-1)/(2*pi)*log( R_KL(1,j+1)/R_KL(1,j) ) +...
                (-TRIG_L(j,2)*TRIG_K(1,1)+TRIG_L(j,1)*TRIG_K(1,2))*...
                B_KL(1,j)/(2*pi)+...
                (TRIG_L(j,1)*TRIG_K(N_P_K,1)+TRIG_L(j,2)*TRIG_K(N_P_K,2))*...
                (-1)/(2*pi)*log( R_KL(N_P_K,j+1)/R_KL(N_P_K,j) ) +...
                (-TRIG_L(j,2)*TRIG_K(N_P_K,1)+TRIG_L(j,1)*TRIG_K(N_P_K,2))*...
                B_KL(N_P_K,j)/(2*pi);
                
    V_KUT_KL=V_KUT_KL+(TRIG_L(j,1)*TRIG_K(1,1)+TRIG_L(j,2)*TRIG_K(1,2))*...
                B_KL(1,j)/(2*pi) +...
                (-TRIG_L(j,2)*TRIG_K(1,1)+TRIG_L(j,1)*TRIG_K(1,2))*...
                (1)/(2*pi)*log( R_KL(1,j+1)/R_KL(1,j) )+...
                (TRIG_L(j,1)*TRIG_K(N_P_K,1)+TRIG_L(j,2)*TRIG_K(N_P_K,2))*...
                B_KL(N_P_K,j)/(2*pi) +...
                (-TRIG_L(j,2)*TRIG_K(N_P_K,1)+TRIG_L(j,1)*TRIG_K(N_P_K,2))*...
                (1)/(2*pi)*log( R_KL(N_P_K,j+1)/R_KL(N_P_K,j) );
end 
% Not to be modified ^^^
S_KUT_BA=S_KUT_KL;
V_KUT_BA=V_KUT_KL;

% S_KUT_BB & V_KUT_BB
N_P_K=N_P_B;
TRIG_K=TRIG_B;
TRIG_L=TRIG_B;
R_KL=R_BB;
B_KL=B_BB;
% Not to be modified vvv
S_KUT_KL=zeros(1,N_P_L);
V_KUT_KL=0.0;
for j=1:N_P_L
    S_KUT_KL(j)=(TRIG_L(j,1)*TRIG_K(1,1)+TRIG_L(j,2)*TRIG_K(1,2))*...
                (-1)/(2*pi)*log( R_KL(1,j+1)/R_KL(1,j) ) +...
                (-TRIG_L(j,2)*TRIG_K(1,1)+TRIG_L(j,1)*TRIG_K(1,2))*...
                B_KL(1,j)/(2*pi)+...
                (TRIG_L(j,1)*TRIG_K(N_P_K,1)+TRIG_L(j,2)*TRIG_K(N_P_K,2))*...
                (-1)/(2*pi)*log( R_KL(N_P_K,j+1)/R_KL(N_P_K,j) ) +...
                (-TRIG_L(j,2)*TRIG_K(N_P_K,1)+TRIG_L(j,1)*TRIG_K(N_P_K,2))*...
                B_KL(N_P_K,j)/(2*pi);
                
    V_KUT_KL=V_KUT_KL+(TRIG_L(j,1)*TRIG_K(1,1)+TRIG_L(j,2)*TRIG_K(1,2))*...
                B_KL(1,j)/(2*pi) +...
                (-TRIG_L(j,2)*TRIG_K(1,1)+TRIG_L(j,1)*TRIG_K(1,2))*...
                (1)/(2*pi)*log( R_KL(1,j+1)/R_KL(1,j) )+...
                (TRIG_L(j,1)*TRIG_K(N_P_K,1)+TRIG_L(j,2)*TRIG_K(N_P_K,2))*...
                B_KL(N_P_K,j)/(2*pi) +...
                (-TRIG_L(j,2)*TRIG_K(N_P_K,1)+TRIG_L(j,1)*TRIG_K(N_P_K,2))*...
                (1)/(2*pi)*log( R_KL(N_P_K,j+1)/R_KL(N_P_K,j) );
end 
% Not to be modified ^^^
S_KUT_BB=S_KUT_KL;
V_KUT_BB=V_KUT_KL;

% S_KUT_BC & V_KUT_BC
N_P_K=N_P_B;
TRIG_K=TRIG_B;
TRIG_L=TRIG_C;
R_KL=R_BC;
B_KL=B_BC;
% Not to be modified vvv
S_KUT_KL=zeros(1,N_P_L);
V_KUT_KL=0.0;
for j=1:N_P_L
    S_KUT_KL(j)=(TRIG_L(j,1)*TRIG_K(1,1)+TRIG_L(j,2)*TRIG_K(1,2))*...
                (-1)/(2*pi)*log( R_KL(1,j+1)/R_KL(1,j) ) +...
                (-TRIG_L(j,2)*TRIG_K(1,1)+TRIG_L(j,1)*TRIG_K(1,2))*...
                B_KL(1,j)/(2*pi)+...
                (TRIG_L(j,1)*TRIG_K(N_P_K,1)+TRIG_L(j,2)*TRIG_K(N_P_K,2))*...
                (-1)/(2*pi)*log( R_KL(N_P_K,j+1)/R_KL(N_P_K,j) ) +...
                (-TRIG_L(j,2)*TRIG_K(N_P_K,1)+TRIG_L(j,1)*TRIG_K(N_P_K,2))*...
                B_KL(N_P_K,j)/(2*pi);
                
    V_KUT_KL=V_KUT_KL+(TRIG_L(j,1)*TRIG_K(1,1)+TRIG_L(j,2)*TRIG_K(1,2))*...
                B_KL(1,j)/(2*pi) +...
                (-TRIG_L(j,2)*TRIG_K(1,1)+TRIG_L(j,1)*TRIG_K(1,2))*...
                (1)/(2*pi)*log( R_KL(1,j+1)/R_KL(1,j) )+...
                (TRIG_L(j,1)*TRIG_K(N_P_K,1)+TRIG_L(j,2)*TRIG_K(N_P_K,2))*...
                B_KL(N_P_K,j)/(2*pi) +...
                (-TRIG_L(j,2)*TRIG_K(N_P_K,1)+TRIG_L(j,1)*TRIG_K(N_P_K,2))*...
                (1)/(2*pi)*log( R_KL(N_P_K,j+1)/R_KL(N_P_K,j) );
end 
% Not to be modified ^^^
S_KUT_BC=S_KUT_KL;
V_KUT_BC=V_KUT_KL;

% S_KUT_CA & V_KUT_CA
N_P_K=N_P_C;
TRIG_K=TRIG_C;
TRIG_L=TRIG_A;
R_KL=R_CA;
B_KL=B_CA;
% Not to be modified vvv
S_KUT_KL=zeros(1,N_P_L);
V_KUT_KL=0.0;
for j=1:N_P_L
    S_KUT_KL(j)=(TRIG_L(j,1)*TRIG_K(1,1)+TRIG_L(j,2)*TRIG_K(1,2))*...
                (-1)/(2*pi)*log( R_KL(1,j+1)/R_KL(1,j) ) +...
                (-TRIG_L(j,2)*TRIG_K(1,1)+TRIG_L(j,1)*TRIG_K(1,2))*...
                B_KL(1,j)/(2*pi)+...
                (TRIG_L(j,1)*TRIG_K(N_P_K,1)+TRIG_L(j,2)*TRIG_K(N_P_K,2))*...
                (-1)/(2*pi)*log( R_KL(N_P_K,j+1)/R_KL(N_P_K,j) ) +...
                (-TRIG_L(j,2)*TRIG_K(N_P_K,1)+TRIG_L(j,1)*TRIG_K(N_P_K,2))*...
                B_KL(N_P_K,j)/(2*pi);
                
    V_KUT_KL=V_KUT_KL+(TRIG_L(j,1)*TRIG_K(1,1)+TRIG_L(j,2)*TRIG_K(1,2))*...
                B_KL(1,j)/(2*pi) +...
                (-TRIG_L(j,2)*TRIG_K(1,1)+TRIG_L(j,1)*TRIG_K(1,2))*...
                (1)/(2*pi)*log( R_KL(1,j+1)/R_KL(1,j) )+...
                (TRIG_L(j,1)*TRIG_K(N_P_K,1)+TRIG_L(j,2)*TRIG_K(N_P_K,2))*...
                B_KL(N_P_K,j)/(2*pi) +...
                (-TRIG_L(j,2)*TRIG_K(N_P_K,1)+TRIG_L(j,1)*TRIG_K(N_P_K,2))*...
                (1)/(2*pi)*log( R_KL(N_P_K,j+1)/R_KL(N_P_K,j) );
end 
% Not to be modified ^^^
S_KUT_CA=S_KUT_KL;
V_KUT_CA=V_KUT_KL;

% S_KUT_CB & V_KUT_CB
N_P_K=N_P_C;
TRIG_K=TRIG_C;
TRIG_L=TRIG_B;
R_KL=R_CB;
B_KL=B_CB;
% Not to be modified vvv
S_KUT_KL=zeros(1,N_P_L);
V_KUT_KL=0.0;
for j=1:N_P_L
    S_KUT_KL(j)=(TRIG_L(j,1)*TRIG_K(1,1)+TRIG_L(j,2)*TRIG_K(1,2))*...
                (-1)/(2*pi)*log( R_KL(1,j+1)/R_KL(1,j) ) +...
                (-TRIG_L(j,2)*TRIG_K(1,1)+TRIG_L(j,1)*TRIG_K(1,2))*...
                B_KL(1,j)/(2*pi)+...
                (TRIG_L(j,1)*TRIG_K(N_P_K,1)+TRIG_L(j,2)*TRIG_K(N_P_K,2))*...
                (-1)/(2*pi)*log( R_KL(N_P_K,j+1)/R_KL(N_P_K,j) ) +...
                (-TRIG_L(j,2)*TRIG_K(N_P_K,1)+TRIG_L(j,1)*TRIG_K(N_P_K,2))*...
                B_KL(N_P_K,j)/(2*pi);
                
    V_KUT_KL=V_KUT_KL+(TRIG_L(j,1)*TRIG_K(1,1)+TRIG_L(j,2)*TRIG_K(1,2))*...
                B_KL(1,j)/(2*pi) +...
                (-TRIG_L(j,2)*TRIG_K(1,1)+TRIG_L(j,1)*TRIG_K(1,2))*...
                (1)/(2*pi)*log( R_KL(1,j+1)/R_KL(1,j) )+...
                (TRIG_L(j,1)*TRIG_K(N_P_K,1)+TRIG_L(j,2)*TRIG_K(N_P_K,2))*...
                B_KL(N_P_K,j)/(2*pi) +...
                (-TRIG_L(j,2)*TRIG_K(N_P_K,1)+TRIG_L(j,1)*TRIG_K(N_P_K,2))*...
                (1)/(2*pi)*log( R_KL(N_P_K,j+1)/R_KL(N_P_K,j) );
end 
% Not to be modified ^^^
S_KUT_CB=S_KUT_KL;
V_KUT_CB=V_KUT_KL;

% S_KUT_CC & V_KUT_CC
N_P_K=N_P_C;
TRIG_K=TRIG_C;
TRIG_L=TRIG_C;
R_KL=R_CC;
B_KL=B_CC;
% Not to be modified vvv
S_KUT_KL=zeros(1,N_P_L);
V_KUT_KL=0.0;
for j=1:N_P_L
    S_KUT_KL(j)=(TRIG_L(j,1)*TRIG_K(1,1)+TRIG_L(j,2)*TRIG_K(1,2))*...
                (-1)/(2*pi)*log( R_KL(1,j+1)/R_KL(1,j) ) +...
                (-TRIG_L(j,2)*TRIG_K(1,1)+TRIG_L(j,1)*TRIG_K(1,2))*...
                B_KL(1,j)/(2*pi)+...
                (TRIG_L(j,1)*TRIG_K(N_P_K,1)+TRIG_L(j,2)*TRIG_K(N_P_K,2))*...
                (-1)/(2*pi)*log( R_KL(N_P_K,j+1)/R_KL(N_P_K,j) ) +...
                (-TRIG_L(j,2)*TRIG_K(N_P_K,1)+TRIG_L(j,1)*TRIG_K(N_P_K,2))*...
                B_KL(N_P_K,j)/(2*pi);
                
    V_KUT_KL=V_KUT_KL+(TRIG_L(j,1)*TRIG_K(1,1)+TRIG_L(j,2)*TRIG_K(1,2))*...
                B_KL(1,j)/(2*pi) +...
                (-TRIG_L(j,2)*TRIG_K(1,1)+TRIG_L(j,1)*TRIG_K(1,2))*...
                (1)/(2*pi)*log( R_KL(1,j+1)/R_KL(1,j) )+...
                (TRIG_L(j,1)*TRIG_K(N_P_K,1)+TRIG_L(j,2)*TRIG_K(N_P_K,2))*...
                B_KL(N_P_K,j)/(2*pi) +...
                (-TRIG_L(j,2)*TRIG_K(N_P_K,1)+TRIG_L(j,1)*TRIG_K(N_P_K,2))*...
                (1)/(2*pi)*log( R_KL(N_P_K,j+1)/R_KL(N_P_K,j) );
end 
% Not to be modified ^^^
S_KUT_CC=S_KUT_KL;
V_KUT_CC=V_KUT_KL;

% B_TAN_A
N_P_K=N_P_A;
TRIG_K=TRIG_A;
% Not to be modified vvv
B_TAN_K=zeros(N_P_K,1);
for i=1:N_P_K
    B_TAN_K(i)=V_INF*(TRIG_K(i,2)*TRIG_ALPHA(1)-TRIG_K(i,1)*TRIG_ALPHA(2));
end
% Not to be modified ^^^
B_TAN_A=B_TAN_K;

% B_TAN_B
N_P_K=N_P_B;
TRIG_K=TRIG_B;
% Not to be modified vvv
B_TAN_K=zeros(N_P_K,1);
for i=1:N_P_K
    B_TAN_K(i)=V_INF*(TRIG_K(i,2)*TRIG_ALPHA(1)-TRIG_K(i,1)*TRIG_ALPHA(2));
end
% Not to be modified ^^^
B_TAN_B=B_TAN_K;


% B_TAN_C
N_P_K=N_P_C;
TRIG_K=TRIG_C;
% Not to be modified vvv
B_TAN_K=zeros(N_P_K,1);
for i=1:N_P_K
    B_TAN_K(i)=V_INF*(TRIG_K(i,2)*TRIG_ALPHA(1)-TRIG_K(i,1)*TRIG_ALPHA(2));
end
% Not to be modified ^^^
B_TAN_C=B_TAN_K;

% B_KUT_A
TRIG_K=TRIG_A;
N_P_K=N_P_A;
% NOT to be modified vvv
B_KUT_K=-V_INF*(TRIG_ALPHA(1)*TRIG_K(N_P_K,1)+TRIG_ALPHA(2)*TRIG_K(N_P_K,2)+...
         TRIG_ALPHA(1)*TRIG_K(1,1)+TRIG_ALPHA(2)*TRIG_K(1,2));
% Not to be modified ^^^
B_KUT_A=B_KUT_K;
     
% B_KUT_B
TRIG_K=TRIG_B;
N_P_K=N_P_B;
% NOT to be modified vvv
B_KUT_K=-V_INF*(TRIG_ALPHA(1)*TRIG_K(N_P_K,1)+TRIG_ALPHA(2)*TRIG_K(N_P_K,2)+...
         TRIG_ALPHA(1)*TRIG_K(1,1)+TRIG_ALPHA(2)*TRIG_K(1,2));
% Not to be modified ^^^
B_KUT_B=B_KUT_K;

% B_KUT_C
TRIG_K=TRIG_C;
N_P_K=N_P_C;
% NOT to be modified vvv
B_KUT_K=-V_INF*(TRIG_ALPHA(1)*TRIG_K(N_P_K,1)+TRIG_ALPHA(2)*TRIG_K(N_P_K,2)+...
         TRIG_ALPHA(1)*TRIG_K(1,1)+TRIG_ALPHA(2)*TRIG_K(1,2));
% Not to be modified ^^^
B_KUT_C=B_KUT_K;

%% ASSEMBLY AND SOLUTION OF THE LINEAR SYSTEM
SYSTEM=[S_TAN_AA , S_TAN_AB , S_TAN_AC , V_TAN_AA , V_TAN_AB , V_TAN_AC ;
        S_TAN_BA , S_TAN_BB , S_TAN_BC , V_TAN_BA , V_TAN_BB , V_TAN_BC ;
        S_TAN_CA , S_TAN_CB , S_TAN_CC , V_TAN_CA , V_TAN_CB , V_TAN_CC ;         
        S_KUT_AA , S_KUT_AB , S_KUT_AC , V_KUT_AA , V_KUT_AB , V_KUT_AC ;
        S_KUT_BA , S_KUT_BB , S_KUT_BC , V_KUT_BA , V_KUT_BB , V_KUT_BC ;
        S_KUT_CA , S_KUT_CB , S_KUT_CC , V_KUT_CA , V_KUT_CB , V_KUT_CC];

RHS= [B_TAN_A;B_TAN_B;B_TAN_C;B_KUT_A;B_KUT_B;B_KUT_C];

%% POTENTIAL FLOW FIELD SOLUTION
X=SYSTEM\RHS;
Q_A=X(1:N_P_A);
Q_B=X(N_P_A+1:N_P_A+N_P_B);
Q_C=X(N_P_A+N_P_B+1:N_P_A+N_P_B+N_P_C);

NI_A=X(N_P_A+N_P_B+N_P_C+1);
NI_B=X(N_P_A+N_P_B+N_P_C+2);
NI_C=X(N_P_A+N_P_B+N_P_C+3);

% Tangential velocity at midpoints

% Vt_A
N_P_K=N_P_A;
TRIG_K=TRIG_A;

% AA
N_P_L=N_P_A;
TRIG_L=TRIG_A;
R_KL=R_AA;
B_KL=B_AA;
% Not to be modified vvv
S_Vt_KL=zeros(N_P_K,N_P_L);
V_Vt_KL=zeros(N_P_K,1);
for i=1:N_P_K
    for j=1:N_P_L
        S_Vt_KL(i,j)=( TRIG_K(i,2)*TRIG_L(j,1)-TRIG_K(i,1)*TRIG_L(j,2) )*...
                     B_KL(i,j)/(2*pi)-...
                     ( TRIG_K(i,1)*TRIG_L(j,1)+TRIG_K(i,2)*TRIG_L(j,2) )*...
                     log( R_KL(i,j+1)/R_KL(i,j) )/(2*pi);    
        
        
        V_Vt_KL(i)=V_Vt_KL(i)+( TRIG_K(i,2)*TRIG_L(j,1)-TRIG_K(i,1)*TRIG_L(j,2) )*...
                              log( R_KL(i,j+1)/R_KL(i,j) )/(2*pi)+...
                              ( TRIG_K(i,1)*TRIG_L(j,1)+TRIG_K(i,2)*TRIG_L(j,2) )*...
                              B_KL(i,j)/(2*pi);        
    end
end
% Not to be modified ^^^
S_Vt_AA=S_Vt_KL;
V_Vt_AA=V_Vt_KL;

% AB
N_P_L=N_P_B;
TRIG_L=TRIG_B;
R_KL=R_AB;
B_KL=B_AB;
% Not to be modified vvv
S_Vt_KL=zeros(N_P_K,N_P_L);
V_Vt_KL=zeros(N_P_K,1);
for i=1:N_P_K
    for j=1:N_P_L
        S_Vt_KL(i,j)=( TRIG_K(i,2)*TRIG_L(j,1)-TRIG_K(i,1)*TRIG_L(j,2) )*...
                     B_KL(i,j)/(2*pi)-...
                     ( TRIG_K(i,1)*TRIG_L(j,1)+TRIG_K(i,2)*TRIG_L(j,2) )*...
                     log( R_KL(i,j+1)/R_KL(i,j) )/(2*pi);    
        
        
        V_Vt_KL(i)=V_Vt_KL(i)+( TRIG_K(i,2)*TRIG_L(j,1)-TRIG_K(i,1)*TRIG_L(j,2) )*...
                              log( R_KL(i,j+1)/R_KL(i,j) )/(2*pi)+...
                              ( TRIG_K(i,1)*TRIG_L(j,1)+TRIG_K(i,2)*TRIG_L(j,2) )*...
                              B_KL(i,j)/(2*pi);        
    end
end
% Not to be modified ^^^
S_Vt_AB=S_Vt_KL;
V_Vt_AB=V_Vt_KL;


% AC
N_P_L=N_P_C;
TRIG_L=TRIG_C;
R_KL=R_AC;
B_KL=B_AC;
% Not to be modified vvv
S_Vt_KL=zeros(N_P_K,N_P_L);
V_Vt_KL=zeros(N_P_K,1);
for i=1:N_P_K
    for j=1:N_P_L
        S_Vt_KL(i,j)=( TRIG_K(i,2)*TRIG_L(j,1)-TRIG_K(i,1)*TRIG_L(j,2) )*...
                     B_KL(i,j)/(2*pi)-...
                     ( TRIG_K(i,1)*TRIG_L(j,1)+TRIG_K(i,2)*TRIG_L(j,2) )*...
                     log( R_KL(i,j+1)/R_KL(i,j) )/(2*pi);    
        
        
        V_Vt_KL(i)=V_Vt_KL(i)+( TRIG_K(i,2)*TRIG_L(j,1)-TRIG_K(i,1)*TRIG_L(j,2) )*...
                              log( R_KL(i,j+1)/R_KL(i,j) )/(2*pi)+...
                              ( TRIG_K(i,1)*TRIG_L(j,1)+TRIG_K(i,2)*TRIG_L(j,2) )*...
                              B_KL(i,j)/(2*pi);        
    end
end
% Not to be modified ^^^
S_Vt_AC=S_Vt_KL;
V_Vt_AC=V_Vt_KL;

Vt_A=V_INF*( TRIG_A(:,1).*TRIG_ALPHA(1)+TRIG_A(:,2).*TRIG_ALPHA(2) )+...
    S_Vt_AA*Q_A + S_Vt_AB*Q_B + S_Vt_AC*Q_C +...
    V_Vt_AA*NI_A+V_Vt_AB*NI_B+V_Vt_AC*NI_C;

% Vt_B
N_P_K=N_P_B;
TRIG_K=TRIG_B;

% BA
N_P_L=N_P_A;
TRIG_L=TRIG_A;
R_KL=R_BA;
B_KL=B_BA;
% Not to be modified vvv
S_Vt_KL=zeros(N_P_K,N_P_L);
V_Vt_KL=zeros(N_P_K,1);
for i=1:N_P_K
    for j=1:N_P_L
        S_Vt_KL(i,j)=( TRIG_K(i,2)*TRIG_L(j,1)-TRIG_K(i,1)*TRIG_L(j,2) )*...
                     B_KL(i,j)/(2*pi)-...
                     ( TRIG_K(i,1)*TRIG_L(j,1)+TRIG_K(i,2)*TRIG_L(j,2) )*...
                     log( R_KL(i,j+1)/R_KL(i,j) )/(2*pi);    
        
        
        V_Vt_KL(i)=V_Vt_KL(i)+( TRIG_K(i,2)*TRIG_L(j,1)-TRIG_K(i,1)*TRIG_L(j,2) )*...
                              log( R_KL(i,j+1)/R_KL(i,j) )/(2*pi)+...
                              ( TRIG_K(i,1)*TRIG_L(j,1)+TRIG_K(i,2)*TRIG_L(j,2) )*...
                              B_KL(i,j)/(2*pi);        
    end
end
% Not to be modified ^^^

S_Vt_BA=S_Vt_KL;
V_Vt_BA=V_Vt_KL;

% BB
N_P_L=N_P_B;
TRIG_L=TRIG_B;
R_KL=R_BB;
B_KL=B_BB;
% Not to be modified vvv
S_Vt_KL=zeros(N_P_K,N_P_L);
V_Vt_KL=zeros(N_P_K,1);
for i=1:N_P_K
    for j=1:N_P_L
        S_Vt_KL(i,j)=( TRIG_K(i,2)*TRIG_L(j,1)-TRIG_K(i,1)*TRIG_L(j,2) )*...
                     B_KL(i,j)/(2*pi)-...
                     ( TRIG_K(i,1)*TRIG_L(j,1)+TRIG_K(i,2)*TRIG_L(j,2) )*...
                     log( R_KL(i,j+1)/R_KL(i,j) )/(2*pi);    
        
        
        V_Vt_KL(i)=V_Vt_KL(i)+( TRIG_K(i,2)*TRIG_L(j,1)-TRIG_K(i,1)*TRIG_L(j,2) )*...
                              log( R_KL(i,j+1)/R_KL(i,j) )/(2*pi)+...
                              ( TRIG_K(i,1)*TRIG_L(j,1)+TRIG_K(i,2)*TRIG_L(j,2) )*...
                              B_KL(i,j)/(2*pi);        
    end
end
% Not to be modified ^^^

S_Vt_BB=S_Vt_KL;
V_Vt_BB=V_Vt_KL;


% BC
N_P_L=N_P_C;
TRIG_L=TRIG_C;
R_KL=R_BC;
B_KL=B_BC;
% Not to be modified vvv
S_Vt_KL=zeros(N_P_K,N_P_L);
V_Vt_KL=zeros(N_P_K,1);
for i=1:N_P_K
    for j=1:N_P_L
        S_Vt_KL(i,j)=( TRIG_K(i,2)*TRIG_L(j,1)-TRIG_K(i,1)*TRIG_L(j,2) )*...
                     B_KL(i,j)/(2*pi)-...
                     ( TRIG_K(i,1)*TRIG_L(j,1)+TRIG_K(i,2)*TRIG_L(j,2) )*...
                     log( R_KL(i,j+1)/R_KL(i,j) )/(2*pi);    
        
        
        V_Vt_KL(i)=V_Vt_KL(i)+( TRIG_K(i,2)*TRIG_L(j,1)-TRIG_K(i,1)*TRIG_L(j,2) )*...
                              log( R_KL(i,j+1)/R_KL(i,j) )/(2*pi)+...
                              ( TRIG_K(i,1)*TRIG_L(j,1)+TRIG_K(i,2)*TRIG_L(j,2) )*...
                              B_KL(i,j)/(2*pi);        
    end
end
% Not to be modified ^^^

S_Vt_BC=S_Vt_KL;
V_Vt_BC=V_Vt_KL;

Vt_B=V_INF*(TRIG_B(:,1).*TRIG_ALPHA(1)+TRIG_B(:,2).*TRIG_ALPHA(2))+...
       S_Vt_BA*Q_A + S_Vt_BB*Q_B + S_Vt_BC*Q_C +...
       V_Vt_BA*NI_A + V_Vt_BB*NI_B + V_Vt_BC*NI_C;

% Vt_C
N_P_K=N_P_C;
TRIG_K=TRIG_C;
% CA
N_P_L=N_P_A;
TRIG_L=TRIG_A;
R_KL=R_CA;
B_KL=B_CA;
% Not to be modified vvv
S_Vt_KL=zeros(N_P_K,N_P_L);
V_Vt_KL=zeros(N_P_K,1);
for i=1:N_P_K
    for j=1:N_P_L
        S_Vt_KL(i,j)=( TRIG_K(i,2)*TRIG_L(j,1)-TRIG_K(i,1)*TRIG_L(j,2) )*...
                     B_KL(i,j)/(2*pi)-...
                     ( TRIG_K(i,1)*TRIG_L(j,1)+TRIG_K(i,2)*TRIG_L(j,2) )*...
                     log( R_KL(i,j+1)/R_KL(i,j) )/(2*pi);    
        
        
        V_Vt_KL(i)=V_Vt_KL(i)+( TRIG_K(i,2)*TRIG_L(j,1)-TRIG_K(i,1)*TRIG_L(j,2) )*...
                              log( R_KL(i,j+1)/R_KL(i,j) )/(2*pi)+...
                              ( TRIG_K(i,1)*TRIG_L(j,1)+TRIG_K(i,2)*TRIG_L(j,2) )*...
                              B_KL(i,j)/(2*pi);        
    end
end
% Not to be modified ^^^

S_Vt_CA=S_Vt_KL;
V_Vt_CA=V_Vt_KL;

% CB
N_P_L=N_P_B;
TRIG_L=TRIG_B;
R_KL=R_CB;
B_KL=B_CB;
% Not to be modified vvv
S_Vt_KL=zeros(N_P_K,N_P_L);
V_Vt_KL=zeros(N_P_K,1);
for i=1:N_P_K
    for j=1:N_P_L
        S_Vt_KL(i,j)=( TRIG_K(i,2)*TRIG_L(j,1)-TRIG_K(i,1)*TRIG_L(j,2) )*...
                     B_KL(i,j)/(2*pi)-...
                     ( TRIG_K(i,1)*TRIG_L(j,1)+TRIG_K(i,2)*TRIG_L(j,2) )*...
                     log( R_KL(i,j+1)/R_KL(i,j) )/(2*pi);    
        
        
        V_Vt_KL(i)=V_Vt_KL(i)+( TRIG_K(i,2)*TRIG_L(j,1)-TRIG_K(i,1)*TRIG_L(j,2) )*...
                              log( R_KL(i,j+1)/R_KL(i,j) )/(2*pi)+...
                              ( TRIG_K(i,1)*TRIG_L(j,1)+TRIG_K(i,2)*TRIG_L(j,2) )*...
                              B_KL(i,j)/(2*pi);        
    end
end
% Not to be modified ^^^

S_Vt_CB=S_Vt_KL;
V_Vt_CB=V_Vt_KL;


% CC
N_P_L=N_P_C;
TRIG_L=TRIG_C;
R_KL=R_CC;
B_KL=B_CC;
% Not to be modified vvv
S_Vt_KL=zeros(N_P_K,N_P_L);
V_Vt_KL=zeros(N_P_K,1);
for i=1:N_P_K
    for j=1:N_P_L
        S_Vt_KL(i,j)=( TRIG_K(i,2)*TRIG_L(j,1)-TRIG_K(i,1)*TRIG_L(j,2) )*...
                     B_KL(i,j)/(2*pi)-...
                     ( TRIG_K(i,1)*TRIG_L(j,1)+TRIG_K(i,2)*TRIG_L(j,2) )*...
                     log( R_KL(i,j+1)/R_KL(i,j) )/(2*pi);    
        
        
        V_Vt_KL(i)=V_Vt_KL(i)+( TRIG_K(i,2)*TRIG_L(j,1)-TRIG_K(i,1)*TRIG_L(j,2) )*...
                              log( R_KL(i,j+1)/R_KL(i,j) )/(2*pi)+...
                              ( TRIG_K(i,1)*TRIG_L(j,1)+TRIG_K(i,2)*TRIG_L(j,2) )*...
                              B_KL(i,j)/(2*pi);        
    end
end
% Not to be modified ^^^

S_Vt_CC=S_Vt_KL;
V_Vt_CC=V_Vt_KL;

Vt_C=V_INF*(TRIG_C(:,1).*TRIG_ALPHA(1)+TRIG_C(:,2).*TRIG_ALPHA(2))+...
    S_Vt_CA*Q_A + S_Vt_CB*Q_B + S_Vt_CC*Q_C +...
    V_Vt_CA*NI_A+V_Vt_CB*NI_B+V_Vt_CC*NI_C;

% Pressure coefficient at midpoints (Bernoulli incompressible)
Cp_A=1-(Vt_A.^2)./V_INF^2;

Cp_B=1-(Vt_B.^2)./V_INF^2;

Cp_C=1-(Vt_C.^2)./V_INF^2;

% Lift coefficient
Cl_A=0.0;
for i=1:N_P_A
    Cl_A=Cl_A-Cp_A(i)*L_P_A(i)/chord_A*TRIG_A(i,1);
end

Cl_B=0.0;
for i=1:N_P_B
    Cl_B=Cl_B-Cp_B(i)*L_P_B(i)/chord_B*TRIG_B(i,1);
end

Cl_C=0.0;
for i=1:N_P_C
    Cl_C=Cl_C-Cp_C(i)*L_P_C(i)/chord_C*TRIG_C(i,1);
end

Cl=Cl_A+Cl_B*chord_B/chord_A+Cl_C*chord_C/chord_A;

%Solution plot
% figure
% plot(MID_A(:,1),-Cp_A,'g',MID_B(:,1),-Cp_B,'b',MID_C(:,1),-Cp_C,'r','linewidth',2);
% hold on
% grid on
% title('-CP','Color','k');
%% STALL PREDICTION ACCORDING TO VALAREZO-CHIN
DELTA_Cp_min_A=abs(min(Cp_A)-Cp_A(N_P_A));
DELTA_Cp_min_B=abs(min(Cp_B)-Cp_B(N_P_B));
DELTA_Cp_min_C=abs(min(Cp_C)-Cp_C(N_P_C)); % WHAT HAPPEN IF I CHOOSE FIRST ...
%                                            PANEL INSTEAD OF LAST? NOTHING
%                                            SHOULD CHANGE!

Cp_VALCHIN=14;

STALL=zeros(1,3);
if DELTA_Cp_min_A>Cp_VALCHIN
    %disp(' MAIN STALL!')
    STALL(1,1)=1;
end

if DELTA_Cp_min_B>Cp_VALCHIN
   % disp(' SECOND STALL!')
    STALL(1,2)=1;
end

if DELTA_Cp_min_C>Cp_VALCHIN
    %disp(' THIRD STALL!')
    STALL(1,3)=1;
end

% disp('AERODYNAMIC ANGLES OF A,B,C')
% disp(rad2deg(ALPHA))
% disp(rad2deg(BETA))
% disp(rad2deg(GAMMA))

%% OUTPUT
OUTPUT.Cp_A=[MID_A(:,1),Cp_A];
OUTPUT.Cp_B=[MID_B(:,1),Cp_B];
OUTPUT.Cp_C=[MID_C(:,1),Cp_C];
OUTPUT.Cl=Cl;
OUTPUT.STALL=STALL;














