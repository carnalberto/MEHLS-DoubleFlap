function[OUTPUT]=HSPM2DMP(V_INF,Cp_VALCHIN,GAP,varargin)
%% HESS-SMITH PANEL METHOD FOR 2D MULTIPLE PROFILES
% INPUT DESCRIPTION
% V_INF: freestream velocity [m/s]

% Cp_VALCHIN: pressure coefficient of maximum lift according to
% Valarezo-Chin criterium

% GAP: array (2 X number of airfoils) [m]. Each column represents x and y
% components of offset of the leading edge of one airfoil with respect to
% the trailing edge of the previous airfoil. The first value is the offset
% of the first leading edge from the origin of the absolute reference frame 
% (see later).

% varargin: ALPHA,A_INPUT,BETA,B_INPUT,GAMMA,C_INPUT,DELTA,D_INPUT...
% Couples of:  angles of attack [deg] , (number of nodes X 2)-arrays [m].   
% The arrays are column-arrays with the coordinates x and y (on each row)
%               of the nodes of each airfoil with respect to a
%               reference frame centered in the leading edge of each
%               airfoil and aligned with the chord of each airfoil. 

% DESCRIPTION:  provides the potential, incompressible, flow solution (Cp,
%               Vt) on the surface of each profile. 
%               The profiles' surfaces are discretized with panels.
%               The flow solution is given at each panel's midpoint.
%               The sources and vortexes strengths per unit length 
%               distributions are uniform on each panel.             
%               The sources strength per unit length varies from panel to
%               panel and and for different profiles.
%               The vortexes strength per unit length is the same for 
%               all panels of a specific profile but different from profile
%               to profile. 
%               Boundary conditions are imposed at panel's midpoints:
%               non-penetration condition for all panles,
%               Kutta's condition for first and last panles (trailing edge)
%               of each airfoil.
%               This is the 'Hess-Smith' choice for panel method.
%               It implemets the Valarezo-Chin criterium (DELTA_Cp < 14
%               between LE and TE) to see if the flow separates on one or
%               more profiles. Computs the Cl.


% Profiles' number n
n=length(varargin)/2;

% Absolute reference frame: 
% triad: right hand triad, x aligned with the freestream towards right 
%                          y vertical                                 
% origin:                                     

% k and l are the profiles' indexes.
% i and j are the panels' indexes.
% i-->k 
% j-->l

%% INPUT
% Initialization
R=zeros(2,2*n);
AOAS=zeros(1,n);
NODES=cell(1,n);
CHORDS=zeros(1,n);
N=zeros(1,n);
ALPHA=0.0; % Because the absolute reference frame is aligned with V_INF
SIZE=0;

for k=1:n
    % Profiles angles with respect to freestream
    AOAS(k)=varargin{2*k-1};

    % Rotation matrixes
    R(:,2*k-1:2*k)=[cos(AOAS(k)),-sin(AOAS(k));
                    sin(AOAS(k)), cos(AOAS(k))];
    
    % Nodes
    NODES{k}=varargin{2*k};
    
    % Number of panels for each airfoil
    N(k)=length(varargin{2*k})-1;
    SIZE=SIZE+N(k);
    
    % Chords
    CHORDS(k)=NODES{k}(1,1);
    
end

% Dimension of the system
SIZE=SIZE+n;

%% GEOMETRY SETUP
% Initialization
LE=zeros(2,n);

% Leading edges positions absolute frame
LE(:,1)=[0.0;0.0]+GAP(:,1);
for k=2:n
    LE(:,k)=LE(:,k-1) + R(:,2*(k-1)-1:2*(k-1))*[CHORDS(k-1);0.0]+GAP(:,k);
end

% Nodes coordinates in absolute frame 
for k=1:n
    for i=1:N(k)+1
        % Rotation
        NODES{k}(i,:)=(R(:,2*k-1:2*k)*NODES{k}(i,:)')';
        
        %Translation
        NODES{k}(i,:)=LE(:,k)' + NODES{k}(i,:);
        
    end
end

% Plot geometry
% figure
% for k=1:n
% plot(NODES{k}(:,1),NODES{k}(:,2),'k','linewidth',2)
% hold on
% end
% axis equal
% grid on
% title('GEOMETRY','Color','k');

%% COMPUTATION OF GEOMETRY PARAMETERS
% Initialization
PANELS=cell(1,n);
COS=cell(1,n);
SIN=cell(1,n);
MID=cell(1,n);
B=cell(n,n);
r=cell(n,n);

for k=1:n
    % Preallocation
    PANELS{k}=zeros(N(k),1);
    COS{k}=zeros(N(k),1);
    SIN{k}=zeros(N(k),1);
    MID{k}=zeros(N(k),2);
    
    for i=1:N(k)
    % Panels' lengths
    PANELS{k}(i)=sqrt( (NODES{k}(i+1,1)-NODES{k}(i,1))^2 + (NODES{k}(i+1,2)-NODES{k}(i,2))^2 ); 
    
    % Trigonometric functions of panels' angles 
    COS{k}(i)=(NODES{k}(i+1,1)-NODES{k}(i,1))/PANELS{k}(i);
    SIN{k}(i)=(NODES{k}(i+1,2)-NODES{k}(i,2))/PANELS{k}(i);
    
    % Panels' midpoints coordinates 
    MID{k}(i,:)=[NODES{k}(i+1,1)+NODES{k}(i,1),NODES{k}(i+1,2)+NODES{k}(i,2)]/2;
    
    end       
end

for k=1:n
    for l=1:n
        % Preallocation
        B{k,l}=zeros(N(k),N(l));
        r{k,l}=zeros(N(k),N(l)+1);
        
        for i=1:N(k)
            for j=1:N(l)
            % Relative positions of midpoint i wrt node j and wrt node j+1  
            % in local components (panel's j reference frame)
            P1_ABSOLUTE=[MID{k}(i,1)-NODES{l}(j,1); MID{k}(i,2)-NODES{l}(j,2)];
            P2_ABSOLUTE=[MID{k}(i,1)-NODES{l}(j+1,1); MID{k}(i,2)-NODES{l}(j+1,2)];
            INV=[COS{l}(j) , SIN{l}(j); -SIN{l}(j) , COS{l}(j)];
            P1_LOCAL=INV*P1_ABSOLUTE;
            P2_LOCAL=INV*P2_ABSOLUTE;
            THETA1=atan2(P1_LOCAL(2),P1_LOCAL(1));
            THETA2=atan2(P2_LOCAL(2),P2_LOCAL(1));
            % Fix in case of autoinduction
            if (abs(THETA1)<10^(-12) && abs(THETA2)>3)
                THETA1=0; 
                THETA2=pi;
            end
            if (abs(THETA2)<10^(-12) && abs(THETA1)>3)
                THETA2=0; 
                THETA1=-pi; 
            end
            
            % B{k,l}(i,j) angle formed by the segments connecting node j on profile l to
            % midpoint i on profile k and node j+1 on profile l.
            B{k,l}(i,j)= THETA2-THETA1;  
            
            % r{k,l}(i,j) distance of midpoint i on profile k from node j on profile l. 
            r{k,l}(i,j)=sqrt( ( MID{k}(i,1)-NODES{l}(j,1) )^2 +...
                              ( MID{k}(i,2)-NODES{l}(j,2) )^2 );
            end
        r{k,l}(i,N(l)+1)=sqrt( ( MID{k}(i,1)-NODES{l}(N(l)+1,1) )^2 +...
                               ( MID{k}(i,2)-NODES{l}(N(l)+1,2) )^2 );
        end
    end
end

%% COMPUTATION OF SYSTEM'S SUBMATRICES
% The sistem of linear equations is obtained by imposing the flow tangency
% condition at each panel's midpoint for each profile and the Kutta
% condition for each profile's trailing edge (at midpoints of first and last panels).
% 
% STRUCTURE OF THE SYSTEM: 
%            [ As , Av ; Ks , Kv ]*[q ; gamma]==[bA ; bK];
%
% DIMENSION: (N(1) +...+ N(k) +...+ N(n) + n) X (N(1) +...+ N(k) +...+ N(n) + n)

% Initialization
As=cell(n,n);
Av=cell(n,n);
Ks=cell(n,n);
Kv=zeros(n,n);
bA=cell(1,n);
bK=zeros(n,1);

% Submatrices computation
for k=1:n
    % Preallocation
    bA{k}=zeros(N(k),1);
    
    for l=1:n
        %Preallocation
        As{k,l}=zeros(N(k),N(l));
        Av{k,l}=zeros(N(k),1);
        Ks{k,l}=zeros(1,N(l));
        
        for i=1:N(k)
             for j=1:N(l)
                % A Submatrices
                As{k,l}(i,j)=(SIN{l}(j)*COS{k}(i)-COS{l}(j)*SIN{k}(i))*...
                              (-1)/(2*pi)*log( r{k,l}(i,j+1)/r{k,l}(i,j) ) +...
                              (SIN{l}(j)*SIN{k}(i)+COS{l}(j)*COS{k}(i))*...
                              B{k,l}(i,j)/(2*pi);

                Av{k,l}(i)= Av{k,l}(i)+... 
                             (SIN{l}(j)*COS{k}(i)-COS{l}(j)*SIN{k}(i))*...
                             B{k,l}(i,j)/(2*pi) + ...
                             (COS{l}(j)*COS{k}(i)+SIN{l}(j)*SIN{k}(i))*...
                             1/(2*pi)*log(r{k,l}(i,j+1)/r{k,l}(i,j));

             end
             % b array
             bA{k}(i)=V_INF*(SIN{k}(i)*cos(ALPHA)-COS{k}(i)*sin(ALPHA));
             
        end
        
        for j=1:N(l)
         % K submatrices
                Ks{k,l}(j)=(COS{l}(j)*COS{k}(1)+SIN{l}(j)*SIN{k}(1))*...
                (-1)/(2*pi)*log( r{k,l}(1,j+1)/r{k,l}(1,j) ) +...
                (-SIN{l}(j)*COS{k}(1)+COS{l}(j)*SIN{k}(1))*...
                B{k,l}(1,j)/(2*pi)+...
                (COS{l}(j)*COS{k}(N(k))+SIN{l}(j)*SIN{k}(N(k)))*...
                (-1)/(2*pi)*log( r{k,l}(N(k),j+1)/r{k,l}(N(k),j) ) +...
                (-SIN{l}(j)*COS{k}(N(k))+COS{l}(j)*SIN{k}(N(k)))*...
                B{k,l}(N(k),j)/(2*pi);
                
                Kv(k,l)=Kv(k,l)+(COS{l}(j)*COS{k}(1)+SIN{l}(j)*SIN{k}(1))*...
                B{k,l}(1,j)/(2*pi) +...
                (-SIN{l}(j)*COS{k}(1)+COS{l}(j)*SIN{k}(1))*...
                (1)/(2*pi)*log( r{k,l}(1,j+1)/r{k,l}(1,j) )+...
                (COS{l}(j)*COS{k}(N(k))+SIN{l}(j)*SIN{k}(N(k)))*...
                B{k,l}(N(k),j)/(2*pi) +...
                (-SIN{l}(j)*COS{k}(N(k))+COS{l}(j)*SIN{k}(N(k)))*...
                (1)/(2*pi)*log( r{k,l}(N(k),j+1)/r{k,l}(N(k),j) );
        end
            
    end
    
    bK(k)=-V_INF*(cos(ALPHA)*COS{k}(N(k))+sin(ALPHA)*SIN{k}(N(k))+...
                    cos(ALPHA)*COS{k}(1)+sin(ALPHA)*SIN{k}(1));

end

%% ASSEMBLY AND SOLUTION OF THE LINEAR SYSTEM
% Initialization
A=zeros(SIZE,SIZE);
b=zeros(SIZE,1);

I=0;
for k=1:n
    J=0;
    for l=1:n
         A(I+1:I+N(k),J+1:J+N(l))=As{k,l};
         A(I+1:I+N(k),SIZE-n+l)=Av{k,l};
         A(SIZE-n+k,J+1:J+N(l))=Ks{k,l};
         A(SIZE-n+k,SIZE-n+l)=Kv(k,l);
         J=J+N(l);
         
    end
    b(I+1:I+N(k))=bA{k};
    b(SIZE-n+k)=bK(k);
    I=I+N(k);
   
end

x=A\b;
q=x(1:SIZE-n);
gamma=x(SIZE-n+1:SIZE);

%% FLOW VELOCITY AND PRESSURE COEFFICIENT
% at midpoints of each panel of each airfoil
% Initialization
Vts=cell(n,n);
Vtv=cell(n,n);

for k=1:n
    for l=1:n
        % Preallocation
        Vts{k,l}=zeros(N(k),N(l));
        Vtv{k,l}=zeros(N(k),1);
        
        for i=1:N(k)
            for j=1:N(l)
            Vts{k,l}(i,j)=( SIN{k}(i)*COS{l}(j)-COS{k}(i)*SIN{l}(j) )*...
                         B{k,l}(i,j)/(2*pi)-...
                         ( COS{k}(i)*COS{l}(j)+SIN{k}(i)*SIN{l}(j) )*...
                         log( r{k,l}(i,j+1)/r{k,l}(i,j) )/(2*pi);    

            Vtv{k,l}(i)=Vtv{k,l}(i)+( SIN{k}(i)*COS{l}(j)-COS{k}(i)*SIN{l}(j) )*...
                         log( r{k,l}(i,j+1)/r{k,l}(i,j) )/(2*pi)+...
                         ( COS{k}(i)*COS{l}(j)+SIN{k}(i)*SIN{l}(j) )*...
                         B{k,l}(i,j)/(2*pi);    
            
            end
        end
    end
end

Vt=cell(1,n);
Cp=cell(1,n);

for k=1:n
    % Preallocation
    Vt{k}=zeros(N(k),1);
    Cp{k}=zeros(N(k),1);
    
    % Initialization
    SUM=zeros(N(k),1);
    J=0;
        for l=1:n
            SUM=SUM+Vts{k,l}*q(J+1:J+N(l))+Vtv{k,l}*gamma(l);
            J=J+N(l);
        end
        
    % Flow velocity    
    Vt{k}=V_INF*( COS{k}.*cos(ALPHA)+SIN{k}.*sin(ALPHA) )+SUM;
    
    % Pressure coefficient at panel's midpoints (Bernoulli incompressible)
    Cp{k}=1-(Vt{k}.^2)./V_INF^2;

end

% Solution plot
% figure
% for k=1:n
% plot(MID{k}(:,1),-Cp{k},'b','linewidth',2);
% grid on
% hold on
% title('-Cp')
% end

%% STALL PREDICTION ACCORDING TO VALAREZO-CHIN
%Initialization
DELTA_Cp_min=zeros(1,n);
STALL=zeros(1,n);

for k=1:n
DELTA_Cp_min(k)=abs(min(Cp{k})-Cp{k}(N(k)));

    if DELTA_Cp_min(k)>Cp_VALCHIN
       STALL(k)=1;
    end
end

%% LIFT COEFFICIENT
% Initialization
Cl=zeros(1,n);

for k=1:n
    for i=1:N(k)
        Cl(k)=Cl(k)-Cp{k}(i)*PANELS{k}(i)/CHORDS(k)*COS{k}(i);
    
    end
end

%% OUTPUT
for k=1:n
OUTPUT.Cp{k}=[MID{k}(:,1),Cp{k}];
OUTPUT.Cl{k}=Cl(k);

end

OUTPUT.STALL=STALL;

