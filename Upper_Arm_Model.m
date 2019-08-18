%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%                ROBOTICS - Dynamics and Control                  %%%%%
%%%%%                         Final Project                           %%%%%
%%%%%                                                                 %%%%%
%%%%%                  Humberto J De las Casas                        %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% All bodies in global coordinates:
%%%%% Axis on Schematic (Anterior View): Left "Z", Up "Y", Out "X".
ground=[0 0 0];
clavicle=[0.006325 0.00693 0.025465];
scapula=[-0.01433 0.02007 0.135535]+clavicle;
humerus=[-0.00955 -0.034 0.009]+scapula;
ulna=[0.0061 -0.2904 -0.0123]+humerus;
radius=[0.0004 -0.011503 0.019999]+ulna;
proximal_row=[0.018 -0.242 0.025]+radius;
hand=[0.003992 -0.015054 0.002327]+proximal_row;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Bodies used for the project: "ground" and "humerus". PLOTS.
[x,y,z]=sphere;
figure; alpha 0.15; axis equal; view(90,0); hold on;
xlabel('X'); ylabel('Y'); zlabel('Z'); camroll(90)
surf(x*0.03,y*0.03,z*0.03,'FaceColor','k')
surf(x*0.03+humerus(1),y*0.03+humerus(2),z*0.03+humerus(3),'FaceColor','b')
surf(x*0.03+ulna(1),y*0.03+ulna(2),z*0.03+ulna(3),...
     'FaceColor','r');
surf(x*0.03+hand(1),y*0.03+hand(2),z*0.03+hand(3),...
     'FaceColor','g'); 
plot3([ground(1) humerus(1)],[ground(2) humerus(2)],...
      [ground(3) humerus(3)],'k','LineWidth',10)
plot3([humerus(1) ulna(1)],[humerus(2) ulna(2)],...
      [humerus(3) ulna(3)],'k','LineWidth',10)
plot3([ulna(1) hand(1)],[ulna(2) hand(2)],...
      [ulna(3) hand(3)],'k','LineWidth',10)
%plot3([ground(1) 0.30],[ground(2) 0],[ground(3) 0],'r','LineWidth',2)
%plot3([ground(1) 0],[ground(2) 0.30],[ground(3) 0],'y','LineWidth',2)
%plot3([ground(1) 0],[ground(2) 0],[ground(3) 0.30],'g','LineWidth',2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Muscles Attachments (10) - INITIAL POSITIONS.
DELT1_A=[0.009 -0.119 0.006]+humerus;           % Deltoid - 1.
DELT1_B=[-0.014 0.011 0.08]+clavicle;
DELT2_A=[0.005 -0.136 0.006]+humerus;           % Deltoid - 2.
DELT2_B=[-0.011 0 0.006]+scapula;
DELT3_A=[-0.056 0.001 -0.025]+scapula;          % Deltoid - 3.
DELT3_B=[0.002 -0.076 0.01]+humerus;
SUPSP_A=[0.003 0.011 0.026]+humerus;            % Supraspinatus.
SUPSP_B=[-0.044 -0.015 -0.059]+scapula;
INFSP_A=[-0.009 0.005 0.024]+humerus;           % Infraspinatus.
INFSP_B=[-0.074 -0.055 -0.048]+scapula;
SUBSC_A=[0.014 0.008 -0.013]+humerus;           % Subscapularis.
SUBSC_B=[-0.072 -0.039 -0.065]+scapula;
TMIN_A=[-0.001 -0.013 0.022]+humerus;           % Teres Minor.
TMIN_B=[-0.096 -0.081 -0.053]+scapula;
TMAJ_A=[0.01 -0.054 -0.006]+humerus;            % Teres Major.
TMAJ_B=[-0.105 -0.103 -0.058]+scapula;
CORB_A=[0.012 -0.041 -0.027]+scapula;           % Coracobra Chialis.
CORB_B=[0.007 -0.15 -0.008]+humerus;
TRIlong_A=[-0.046 -0.041 -0.014]+scapula;       % Triceps Brachi Long Head.
TRIlong_B=[-0.022 0.01 -0.001]+ulna;
%%%%% Muscles Attachments (14) - INITIAL POSITIONS.
PEC1_A=[0.012 -0.042 0.008]+humerus;          % Pectoralis Major - 1.
PEC1_B=[0.003 0 0.051]+clavicle;
PEC2_A=[0.013 -0.043 0.008]+humerus;          % Pectoralis Major - 2.
PEC2_B=[0.028 -0.045 0.023];
PEC3_A=[0.013 -0.044 0.008]+humerus;          % Pectoralis Major - 3.
PEC3_B=[0.057 -0.117 0.038];
LAT1_A=[0.01 -0.034 -0.007]+humerus;          % Latissimus Dorsi - 1.
LAT1_B=[-0.096 -0.117 0.009];
LAT2_A=[0.01 -0.041 -0.006]+humerus;          % Latissimus Dorsi - 2.
LAT2_B=[-0.072 -0.188 0.008];
LAT3_A=[0.012 -0.039 -0.004]+humerus;          % Latissimus Dorsi - 3.
LAT3_B=[-0.071 -0.249 0.009];
TRIlat_A=[-0.006 -0.126 0.004]+humerus;        % Triceps Brachi Lateral Head.
TRIlat_B=[-0.022 0.01 -0.001]+ulna;
TRImed_A=[-0.006 -0.126 0.004]+humerus;        % Triceps Brachi Medium Head.
TRImed_B=[-0.022 0.01 -0.001]+ulna;
BIClong_A=[-0.031 -0.024 -0.013]+scapula;       % Biceps - Long.
BIClong_B=[-0.002 -0.038 -0.002]+radius;
BICshort_A=[-0.031 -0.024 -0.013]+scapula;      % Biceps - Short.
BICshort_B=[-0.002 -0.038 -0.002]+radius;
BRA_A=[0.007 -0.174 -0.004]+humerus;           % Brachialis.
BRA_B=[-0.003 -0.024 0.001]+ulna;
BRD_A=[-0.01 -0.2 0.002]+humerus;           % Brachioradialis.
BRD_B=[0.042 -0.221 0.022]+radius;
ANC_A=[-0.007 -0.284 0.01]+humerus;           % Anconeus.
ANC_B=[-0.025 -0.001 0.006]+ulna;
PT_A=[0.004 -0.276 -0.036]+humerus;            % Pronator Teres.
PT_B=[0.025 -0.109 0.02]+radius;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Ploting the Muscles Attachments.
plot3([DELT1_A(1) DELT1_B(1)],[DELT1_A(2) DELT1_B(2)],...
      [DELT1_A(3) DELT1_B(3)],'r','LineWidth',3);
plot3([DELT2_A(1) DELT2_B(1)],[DELT2_A(2) DELT2_B(2)],...
      [DELT2_A(3) DELT2_B(3)],'r','LineWidth',3);
plot3([DELT3_A(1) DELT3_B(1)],[DELT3_A(2) DELT3_B(2)],...
      [DELT3_A(3) DELT3_B(3)],'r','LineWidth',3);
plot3([SUPSP_A(1) SUPSP_B(1)],[SUPSP_A(2) SUPSP_B(2)],...
      [SUPSP_A(3) SUPSP_B(3)],'r','LineWidth',3);
plot3([INFSP_A(1) INFSP_B(1)],[INFSP_A(2) INFSP_B(2)],...
      [INFSP_A(3) INFSP_B(3)],'r','LineWidth',3);
plot3([SUBSC_A(1) SUBSC_B(1)],[SUBSC_A(2) SUBSC_B(2)],...
      [SUBSC_A(3) SUBSC_B(3)],'r','LineWidth',3);
plot3([TMIN_A(1) TMIN_B(1)],[TMIN_A(2) TMIN_B(2)],...
      [TMIN_A(3) TMIN_B(3)],'r','LineWidth',3);
plot3([TMAJ_A(1) TMAJ_B(1)],[TMAJ_A(2) TMAJ_B(2)],...
      [TMAJ_A(3) TMAJ_B(3)],'r','LineWidth',3);
plot3([CORB_A(1) CORB_B(1)],[CORB_A(2) CORB_B(2)],...
      [CORB_A(3) CORB_B(3)],'r','LineWidth',3);
plot3([TRIlong_A(1) TRIlong_B(1)],[TRIlong_A(2) TRIlong_B(2)],...
      [TRIlong_A(3) TRIlong_B(3)],'r','LineWidth',3);
plot3([PEC1_A(1) PEC1_B(1)],[PEC1_A(2) PEC1_B(2)],...
      [PEC1_A(3) PEC1_B(3)],'r','LineWidth',3);
plot3([PEC2_A(1) PEC2_B(1)],[PEC2_A(2) PEC2_B(2)],...
      [PEC2_A(3) PEC2_B(3)],'r','LineWidth',3);
plot3([PEC3_A(1) PEC3_B(1)],[PEC3_A(2) PEC3_B(2)],...
      [PEC3_A(3) PEC3_B(3)],'r','LineWidth',3);
plot3([LAT1_A(1) LAT1_B(1)],[LAT1_A(2) LAT1_B(2)],...
      [LAT1_A(3) LAT1_B(3)],'r','LineWidth',3);
plot3([LAT2_A(1) LAT2_B(1)],[LAT2_A(2) LAT2_B(2)],...
      [LAT2_A(3) LAT2_B(3)],'r','LineWidth',3);
plot3([LAT3_A(1) LAT3_B(1)],[LAT3_A(2) LAT3_B(2)],...
      [LAT3_A(3) LAT3_B(3)],'r','LineWidth',3);
plot3([TRIlat_A(1) TRIlat_B(1)],[TRIlat_A(2) TRIlat_B(2)],...
      [TRIlat_A(3) TRIlat_B(3)],'r','LineWidth',3);
plot3([TRImed_A(1) TRImed_B(1)],[TRImed_A(2) TRImed_B(2)],...
      [TRImed_A(3) TRImed_B(3)],'r','LineWidth',3);
plot3([BIClong_A(1) BIClong_B(1)],[BIClong_A(2) BIClong_B(2)],...
      [BIClong_A(3) BIClong_B(3)],'r','LineWidth',3);
plot3([BICshort_A(1) BICshort_B(1)],[BICshort_A(2) BICshort_B(2)],...
      [BICshort_A(3) BICshort_B(3)],'r','LineWidth',3);
plot3([BRA_A(1) BRA_B(1)],[BRA_A(2) BRA_B(2)],...
      [BRA_A(3) BRA_B(3)],'r','LineWidth',3);
plot3([BRD_A(1) BRD_B(1)],[BRD_A(2) BRD_B(2)],...
      [BRD_A(3) BRD_B(3)],'r','LineWidth',3);
plot3([ANC_A(1) ANC_B(1)],[ANC_A(2) ANC_B(2)],...
      [ANC_A(3) ANC_B(3)],'r','LineWidth',3);
plot3([PT_A(1) PT_B(1)],[PT_A(2) PT_B(2)],...
      [PT_A(3) PT_B(3)],'r','LineWidth',3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Symbolic Variables. 
syms q1 q2 q3 q4 q5 q1dot q2dot q3dot q4dot q5dot 
syms q1ddot q2ddot q3ddot q4ddot q5ddot g
qdot=[q1dot;q2dot;q3dot;q4dot;q5dot];
qddot=[q1ddot;q2ddot;q3ddot;q4ddot;q5ddot];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Frames transformation on function of q=[q1 q2 q3 q4 q5].  
Q_00=eye(4)+[zeros(4,3) [ground(1);ground(2);ground(3);0]];
%%%%% Shoulder:
Q0=(eye(4)+[zeros(4,3) [humerus(1);humerus(2);humerus(3);0]]);
Q1=Q0*[cos(q1) -sin(q1) 0 0;sin(q1) cos(q1) 0 0 ...
        ;0 0 1 0;0 0 0 1]*[1 0 0 0;0 0 1 0;0 -1 0 0;0 0 0 1];
Q2=Q1*...
        [cos(q2) -sin(q2) 0 0;sin(q2) cos(q2) ...
        0 0;0 0 1 0;0 0 0 1]*[0 -1 0 0;1 0 0 0;0 0 1 0; 0 0 0 1]*...
        [1 0 0 0;0 0 -1 0;0 1 0 0; 0 0 0 1];
Q3=Q2*[cos(q3) -sin(q3) ...
        0 0;sin(q3) cos(q3) 0 0;0 0 1 0;0 0 0 1];
%%%%% Elbow:
Q4=Q3*(eye(4)+[zeros(4,3) [-ulna(3)+humerus(3);...
   ulna(2)-humerus(2);ulna(1)-humerus(1);0]])*[0 -1 0 0;1 0 0 0;...
   0 0 1 0;0 0 0 1]*[1 0 0 0;0 0 1 0;0 -1 0 0;0 0 0 1]*[cos(q4) -sin(q4)...
   0 0;sin(q4) cos(q4) 0 0;0 0 1 0;0 0 0 1];
%%%%% Forearm:
Q5=Q4*[0 -1 0 0;1 0 0 0;0 0 1 0;0 0 0 1]*[1 0 0 0;0 0 1 0;0 -1 0 0;...
   0 0 0 1]*[cos(q5) -sin(q5) 0 0;sin(q5) cos(q5) 0 0;0 0 1 0;0 0 0 1];
Q_Radius=Q5*(eye(4)+[zeros(4,3) [-radius(1)+ulna(1);...
        -radius(3)+ulna(3);-radius(2)+ulna(2);0]]);
Q_PR=Q_Radius*(eye(4)+[zeros(4,3) [-proximal_row(1)+radius(1);...
        -proximal_row(3)+radius(3);-proximal_row(2)+radius(2);0]]);
Q_Hand=Q_PR*(eye(4)+[zeros(4,3) [-hand(1)+proximal_row(1);...
        -hand(3)+proximal_row(3);-hand(2)+proximal_row(2);0]]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Jacobian respect to the center of mass.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Physical parameters - Humerus:
Humerus_mass=1.99757;
Humerus_mass_center=Q3*(eye(4)+[0.012746;-0.140141;0.018064;0]);
Humerus_mass_center=Humerus_mass_center(1:3,4);
Humerus_I=[0.01227763 -3.4741E-4 -2.325E-4;-3.4741E-4 0.00255133 0.0012293; ...
          -2.325E-4 0.0012293 0.01257888];                                   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Physical parameters - Ulna, Radius and Hand:
Ulna_mass=1.1053;
Radius_mass=0.23359;
Hand_mass=0.5819;
Ulna_mass_center=Q5*(eye(4)+[-0.00971783;-0.024286;0.0959509;0]);
Radius_mass_center=Q_Radius*(eye(4)+[-0.0336341;-0.0156;0.181559;0]);
Hand_mass_center=Q_Hand*(eye(4)+[0.00301314;0.00112205;0.0424993;0]);
Ulna_mass_center=Ulna_mass_center(1:3,4);
Radius_mass_center=Radius_mass_center(1:3,4);
Hand_mass_center=Hand_mass_center(1:3,4);
Ulna_I=[0.00541309 3.1686E-4 -7.615E-5;3.1686E-4 ...
        0.00115318 0.00109169;-7.615E-5 0.00109169 0.00494361];
Radius_I=[4.3855E-4 3.014E-5 -4.24E-6; 3.014E-5 8.859E-5 ...
        6.418E-5; -4.24E-6 6.418E-5 4.0258E-4];
Hand_I=[1.1E-4 9.0E-7 -2.0E-7;9.0E-7 6.0E-5 1.2E-5; ...
        -2.0E-7 1.2E-5 1.5E-4];
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Due to the model only has 3 DOFs on the humerus socket-ball joint, 
%%%%% the Z vectors direction are the same for all frames.
Z_0=Q0(1:3,3);
Z_1=Q1(1:3,3);
Z_2=Q2(1:3,3);
Z_3=Q3(1:3,3);
Z_4=Q4(1:3,3);
Z_5=Q5(1:3,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Jacobian:
J1=cross(Z_0,Humerus_mass_center-Q0(1:3,4));
J1(1:3,2:5)=zeros(3,4);
J1(4:6,1:5)=[Z_0 zeros(3,4)];
J1=simplify(J1);

J2=cross(Z_0,Humerus_mass_center-Q0(1:3,4));
J2(1:3,2)=cross(Z_1,Humerus_mass_center-Q0(1:3,4));
J2(1:3,3:5)=zeros(3,3);
J2(4:6,1:5)=[Z_0 Z_1 zeros(3,3)];
J2=simplify(J2);

J3=cross(Z_0,Humerus_mass_center-Q0(1:3,4));
J3(1:3,2)=cross(Z_1,Humerus_mass_center-Q0(1:3,4));
J3(1:3,3)=cross(Z_2,Humerus_mass_center-Q0(1:3,4));
J3(1:3,4:5)=zeros(3,2);
J3(4:6,1:5)=[Z_0 Z_1 Z_2 zeros(3,2)];
J3=simplify(J3);

J4_1=cross(Z_0,Ulna_mass_center-Q0(1:3,4));
J4_1(1:3,2)=cross(Z_1,Ulna_mass_center-Q0(1:3,4));
J4_1(1:3,3)=cross(Z_2,Ulna_mass_center-Q0(1:3,4));
J4_1(1:3,4)=cross(Z_3,Ulna_mass_center-Q4(1:3,4));
J4_1(1:3,5)=zeros(3,1);
J4_1(4:6,1:5)=[Z_0 Z_1 Z_2 Z_3 zeros(3,1)];
J4_1=simplify(J4_1);

J4_2=cross(Z_0,Radius_mass_center-Q0(1:3,4));
J4_2(1:3,2)=cross(Z_1,Radius_mass_center-Q0(1:3,4));
J4_2(1:3,3)=cross(Z_2,Radius_mass_center-Q0(1:3,4));
J4_2(1:3,4)=cross(Z_3,Radius_mass_center-Q4(1:3,4));
J4_2(1:3,5)=zeros(3,1);
J4_2(4:6,1:5)=[Z_0 Z_1 Z_2 Z_3 zeros(3,1)];
J4_2=simplify(J4_2);

J4_3=cross(Z_0,Hand_mass_center-Q0(1:3,4));
J4_3(1:3,2)=cross(Z_1,Hand_mass_center-Q0(1:3,4));
J4_3(1:3,3)=cross(Z_2,Hand_mass_center-Q0(1:3,4));
J4_3(1:3,4)=cross(Z_3,Hand_mass_center-Q4(1:3,4));
J4_3(1:3,5)=zeros(3,1);
J4_3(4:6,1:5)=[Z_0 Z_1 Z_2 Z_3 zeros(3,1)];
J4_3=simplify(J4_3);

J5_1=cross(Z_0,Ulna_mass_center-Q0(1:3,4));
J5_1(1:3,2)=cross(Z_1,Ulna_mass_center-Q0(1:3,4));
J5_1(1:3,3)=cross(Z_2,Ulna_mass_center-Q0(1:3,4));
J5_1(1:3,4)=cross(Z_3,Ulna_mass_center-Q4(1:3,4));
J5_1(1:3,5)=cross(Z_4,Ulna_mass_center-Q4(1:3,4));
J5_1(4:6,1:5)=[Z_0 Z_1 Z_2 Z_3 Z_4];
J5_1=simplify(J5_1);

J5_2=cross(Z_0,Radius_mass_center-Q0(1:3,4));
J5_2(1:3,2)=cross(Z_1,Radius_mass_center-Q0(1:3,4));
J5_2(1:3,3)=cross(Z_2,Radius_mass_center-Q0(1:3,4));
J5_2(1:3,4)=cross(Z_3,Radius_mass_center-Q4(1:3,4));
J5_2(1:3,5)=cross(Z_4,Radius_mass_center-Q4(1:3,4));
J5_2(4:6,1:5)=[Z_0 Z_1 Z_2 Z_3 Z_4];
J5_2=simplify(J5_2);

J5_3=cross(Z_0,Hand_mass_center-Q0(1:3,4));
J5_3(1:3,2)=cross(Z_1,Hand_mass_center-Q0(1:3,4));
J5_3(1:3,3)=cross(Z_2,Hand_mass_center-Q0(1:3,4));
J5_3(1:3,4)=cross(Z_3,Hand_mass_center-Q4(1:3,4));
J5_3(1:3,5)=cross(Z_4,Hand_mass_center-Q4(1:3,4));
J5_3(4:6,1:5)=[Z_0 Z_1 Z_2 Z_3 Z_4];
J5_3=simplify(J5_3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Mass Matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Mass Matrix.
R1=Q1(1:3,1:3);
R2=Q2(1:3,1:3);
R3=Q3(1:3,1:3);
R4=Q4(1:3,1:3);
R5=Q5(1:3,1:3);
M_q=simplify(...
            J3(4:6,:).'*R3*Humerus_I*R3.'*J3(4:6,:)+...
            J5_1(4:6,:).'*R5*Ulna_I*R5.'*J5_1(4:6,:)+...
            J5_2(4:6,:).'*R5*Radius_I*R5.'*J5_1(4:6,:)+...
            J5_3(4:6,:).'*R5*Hand_I*R5.'*J5_1(4:6,:)+...
            Humerus_mass*J3(1:3,:).'*J3(1:3,:)+...
            Ulna_mass*J5_1(1:3,:).'*J5_1(1:3,:)+...
            Radius_mass*J5_2(1:3,:).'*J5_2(1:3,:)+...
            Hand_mass*J5_3(1:3,:).'*J5_3(1:3,:));
               
% isequal(M_q,M_q.');                     % It is symmetric!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Coriolis Matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             
q=[q1;q2;q3;q4;q5];
for i=1:5
    for j=1:5
        for k=1:5
c(i,j,k)=((diff(M_q(k,j),q(i)))+(diff(M_q(k,i),q(j)))...
        -(diff(M_q(i,j),q(k))))/2;
        end
    end
end
for k=1:5
    for j=1:5
C_qqdot(k,j)=c(1,j,k)*q1dot+c(2,j,k)*q2dot+c(3,j,k)*q3dot+...
             c(4,j,k)*q4dot+c(5,j,k)*q5dot;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Skew-symmetric Test.
% for i=1:5
%     for j=1:5
% M_DOT(i,j)=diff(M_q(i,j),q1)*q1dot+diff(M_q(i,j),q2)*...
%            q2dot+diff(M_q(i,j),q3)*q3dot+diff(M_q(i,j),q4)*q4dot+...
%            +diff(M_q(i,j),q5)*q5dot;
%     end
% end
% N_qqdot=simplify(M_DOT-2*C_qqdot);
% isequal(N_qqdot,-N_qqdot.');                % It is skew-symmetric!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Gravity Matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
g_v=[0;0;-g];
P_Humerus=-Humerus_mass*g_v.'*Humerus_mass_center;
P_Ulna=-Ulna_mass*g_v.'*Ulna_mass_center;
P_Radius=-Radius_mass*g_v.'*Radius_mass_center;
P_Hand=-Hand_mass*g_v.'*Hand_mass_center;
P_Total=P_Humerus+P_Ulna+P_Radius+P_Hand;

g_q(1)=simplify(diff(P_Total,q1));
g_q(2,1)=simplify(diff(P_Total,q2));
g_q(3,1)=simplify(diff(P_Total,q3)); 
g_q(4,1)=simplify(diff(P_Total,q4));
g_q(5,1)=simplify(diff(P_Total,q5));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% PART II
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Muscles Attachments (10) - IN FCN OF q=[q1 q2 q3].
DELT1_A=Q3*(eye(4)+[zeros(4,3) [-0.006;-0.119;0.009;0]]);
DELT1_B=[-0.014 0.011 0.08]+clavicle;
DELT2_A=Q3*(eye(4)+[zeros(4,3) [-0.006;-0.136;0.005;0]]);
DELT2_B=[-0.011 0 0.006]+scapula;
DELT3_A=[-0.056 0.001 -0.025]+scapula;
DELT3_B=Q3*(eye(4)+[zeros(4,3) [-0.01;-0.076;0.002;0]]);
SUPSP_A=Q3*(eye(4)+[zeros(4,3) [-0.026;0.011;0.003;0]]);
SUPSP_B=[-0.044 -0.015 -0.059]+scapula;
INFSP_A=Q3*(eye(4)+[zeros(4,3) [-0.024;0.005;-0.009;0]]);
INFSP_B=[-0.074 -0.055 -0.048]+scapula;
SUBSC_A=Q3*(eye(4)+[zeros(4,3) [0.013;0.008;0.014;0]]);
SUBSC_B=[-0.072 -0.039 -0.065]+scapula;
TMIN_A=Q3*(eye(4)+[zeros(4,3) [-0.022;-0.013;-0.001;0]]);
TMIN_B=[-0.096 -0.081 -0.053]+scapula;
TMAJ_A=Q3*(eye(4)+[zeros(4,3) [0.006;-0.054;0.01;0]]);
TMAJ_B=[-0.105 -0.103 -0.058]+scapula;
CORB_A=[0.012 -0.041 -0.027]+scapula;
CORB_B=Q3*(eye(4)+[zeros(4,3) [0.008;-0.15;0.007;0]]);
TRIlong_A=[-0.046 -0.041 -0.014]+scapula;
TRIlong_B=Q5*(eye(4)+[zeros(4,3) [0.022;0.001;-0.01;0]]);
PEC1_A=Q3*(eye(4)+[zeros(4,3) [-0.008;-0.042;0.012;0]]);
PEC1_B=[0.003 0 0.051]+clavicle;
PEC2_A=Q3*(eye(4)+[zeros(4,3) [-0.008;-0.043;0.013;0]]);
PEC2_B=[0.028 -0.045 0.023];
PEC3_A=Q3*(eye(4)+[zeros(4,3) [-0.008;-0.044;0.013;0]]);
PEC3_B=[0.057 -0.117 0.038];
LAT1_A=Q3*(eye(4)+[zeros(4,3) [0.007;-0.034;0.01;0]]);
LAT1_B=[-0.096 -0.117 0.009];
LAT2_A=Q3*(eye(4)+[zeros(4,3) [0.006;-0.041;0.01;0]]);
LAT2_B=[-0.072 -0.188 0.008];
LAT3_A=Q3*(eye(4)+[zeros(4,3) [0.004;-0.039;0.012;0]]);
LAT3_B=[-0.071 -0.249 0.009];
TRIlat_A=Q3*(eye(4)+[zeros(4,3) [-0.004;-0.126;-0.006;0]]);
TRIlat_B=Q5*(eye(4)+[zeros(4,3) [0.022;0.001;-0.01;0]]);
TRImed_A=Q3*(eye(4)+[zeros(4,3) [-0.004;-0.126;-0.006;0]]);
TRImed_B=Q5*(eye(4)+[zeros(4,3) [0.022;0.001;-0.01;0]]);
BIClong_A=[-0.031 -0.024 -0.013]+scapula;
BIClong_B=Q5*(eye(4)+[zeros(4,3) [0.002;0.002;0.038;0]]);
BICshort_A=[-0.031 -0.024 -0.013]+scapula; 
BICshort_B=Q5*(eye(4)+[zeros(4,3) [0.002;0.002;0.038;0]]);
BRA_A=Q3*(eye(4)+[zeros(4,3) [0.004;-0.174;0.007;0]]);       
BRA_B=Q5*(eye(4)+[zeros(4,3) [0.003;-0.001;0.024;0]]);
BRD_A=Q3*(eye(4)+[zeros(4,3) [-0.002;-0.2;-0.01;0]]);          
BRD_B=Q5*(eye(4)+[zeros(4,3) [-0.042;-0.022;0.221;0]]);
ANC_A=Q3*(eye(4)+[zeros(4,3) [-0.01;-0.284;-0.007;0]]);        
ANC_B=Q5*(eye(4)+[zeros(4,3) [0.025;-0.006;0.001;0]]);
PT_A=Q3*(eye(4)+[zeros(4,3) [0.036;-0.276;0.004;0]]);        
PT_B=Q5*(eye(4)+[zeros(4,3) [-0.025;-0.02;0.109;0]]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Computing muscle lengths:
L_DELT1=simplify(norm(DELT1_A(1:3,4)-DELT1_B.'));
L_DELT2=simplify(norm(DELT2_A(1:3,4)-DELT2_B.'));
L_DELT3=simplify(norm(DELT3_A.'-DELT3_B(1:3,4)));
L_SUPSP=simplify(norm(SUPSP_A(1:3,4)-SUPSP_B.'));
L_INFSP=simplify(norm(INFSP_A(1:3,4)-INFSP_B.'));
L_SUBSC=simplify(norm(SUBSC_A(1:3,4)-SUBSC_B.'));
L_TMIN=simplify(norm(TMIN_A(1:3,4)-TMIN_B.'));
L_TMAJ=simplify(norm(TMAJ_A(1:3,4)-TMAJ_B.'));
L_CORB=simplify(norm(CORB_A.'-CORB_B(1:3,4)));
L_TRIlong=simplify(norm(TRIlong_A.'-TRIlong_B(1:3,4)));
L_PEC1=simplify(norm(PEC1_A(1:3,4)-PEC1_B.'));
L_PEC2=simplify(norm(PEC2_A(1:3,4)-PEC2_B.'));
L_PEC3=simplify(norm(PEC3_A(1:3,4)-PEC3_B.'));
L_LAT1=simplify(norm(LAT1_A(1:3,4)-LAT1_B.'));
L_LAT2=simplify(norm(LAT2_A(1:3,4)-LAT2_B.'));
L_LAT3=simplify(norm(LAT3_A(1:3,4)-LAT3_B.'));
L_TRIlat=simplify(norm(TRIlat_A(1:3,4)-TRIlat_B(1:3,4)));
L_TRImed=simplify(norm(TRImed_A(1:3,4)-TRImed_B(1:3,4)));
L_BIClong=simplify(norm(BIClong_A.'-BIClong_B(1:3,4)));
L_BICshort=simplify(norm(BICshort_A.'-BICshort_B(1:3,4)));
L_BRA=simplify(norm(BRA_A(1:3,4)-BRA_B(1:3,4)));
L_BRD=simplify(norm(BRD_A(1:3,4)-BRD_B(1:3,4)));
L_ANC=simplify(norm(ANC_A(1:3,4)-ANC_B(1:3,4)));
L_PT=simplify(norm(PT_A(1:3,4)-PT_B(1:3,4)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Computing unit vector of muscles:
Unit_DELT1=simplify((DELT1_A(1:3,4)-DELT1_B.')/L_DELT1);
Unit_DELT2=simplify((DELT2_A(1:3,4)-DELT2_B.')/L_DELT2);
Unit_DELT3=simplify((DELT3_A.'-DELT3_B(1:3,4))/L_DELT3);
Unit_SUPSP=simplify((SUPSP_A(1:3,4)-SUPSP_B.')/L_SUPSP);
Unit_INFSP=simplify((INFSP_A(1:3,4)-INFSP_B.')/L_INFSP);
Unit_SUBSC=simplify((SUBSC_A(1:3,4)-SUBSC_B.')/L_SUBSC);
Unit_TMIN=simplify((TMIN_A(1:3,4)-TMIN_B.')/L_TMIN);
Unit_TMAJ=simplify((TMAJ_A(1:3,4)-TMAJ_B.')/L_TMAJ);
Unit_CORB=simplify((CORB_A.'-CORB_B(1:3,4))/L_CORB);
Unit_TRIlong=simplify((+TRIlong_A.'-TRIlong_B(1:3,4))/L_TRIlong);
Unit_PEC1=simplify((PEC1_A(1:3,4)-PEC1_B.')/L_PEC1);
Unit_PEC2=simplify((PEC2_A(1:3,4)-PEC2_B.')/L_PEC2);
Unit_PEC3=simplify((PEC3_A(1:3,4)-PEC3_B.')/L_PEC3);
Unit_LAT1=simplify((LAT1_A(1:3,4)-LAT1_B.')/L_LAT1);
Unit_LAT2=simplify((LAT2_A(1:3,4)-LAT2_B.')/L_LAT2);
Unit_LAT3=simplify((LAT3_A(1:3,4)-LAT3_B.')/L_LAT3);
Unit_TRIlat=simplify((TRIlat_A(1:3,4)-TRIlat_B(1:3,4))/L_TRIlat);
Unit_TRImed=simplify((TRImed_A(1:3,4)-TRImed_B(1:3,4))/L_TRImed);
Unit_BIClong=simplify((BIClong_A.'-BIClong_B(1:3,4))/L_BIClong);
Unit_BICshort=simplify((BICshort_A.'-BICshort_B(1:3,4))/L_BICshort);
Unit_BRA=simplify((BRA_A(1:3,4)-BRA_B(1:3,4))/L_BRA);
Unit_BRD=simplify((BRD_A(1:3,4)-BRD_B(1:3,4))/L_BRD);
Unit_ANC=simplify((ANC_A(1:3,4)-ANC_B(1:3,4))/L_ANC);
Unit_PT=simplify((PT_A(1:3,4)-PT_B(1:3,4))/L_PT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Optimal lengths:
L_opt_DELT1=0.0976;
L_opt_DELT2=0.1078;
L_opt_DELT3=0.1367;
L_opt_SUPSP=0.0682;
L_opt_INFSP=0.0755;
L_opt_SUBSC=0.0873;
L_opt_TMIN=0.0741;
L_opt_TMAJ=0.1624;
L_opt_CORB=0.0932;
L_opt_TRIlong=0.134;
L_opt_PEC1=0.1442;
L_opt_PEC2=0.1385;
L_opt_PEC3=0.1385;
L_opt_LAT1=0.254;
L_opt_LAT2=0.2324;
L_opt_LAT3=0.2789;
L_opt_TRIlat=0.1138;
L_opt_TRImed=0.1138;
L_opt_BIClong=0.1157;
L_opt_BICshort=0.1321;
L_opt_BRA=0.0858;
L_opt_BRD=0.1726;
L_opt_ANC=0.027;
L_opt_PT=0.0492;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Initial lengths:
q1=0;q2=0;q3=0;
L_0_DELT1=eval(L_DELT1);                                                        
L_0_DELT2=eval(L_DELT2);
L_0_DELT3=eval(L_DELT3);
L_0_SUPSP=eval(L_SUPSP);
L_0_INFSP=eval(L_INFSP);
L_0_SUBSC=eval(L_SUBSC);
L_0_TMIN=eval(L_TMIN);
L_0_TMAJ=eval(L_TMAJ);
L_0_CORB=eval(L_CORB);
L_0_TRIlong=eval(L_TRIlong);
L_0_PEC1=eval(L_PEC1);
L_0_PEC2=eval(L_PEC2);
L_0_PEC3=eval(L_PEC3);
L_0_LAT1=eval(L_LAT1);
L_0_LAT2=eval(L_LAT2);
L_0_LAT3=eval(L_LAT3);
L_0_TRIlat=eval(L_TRIlat);
L_0_TRImed=eval(L_TRImed);
L_0_BIClong=eval(L_BIClong);
L_0_BICshort=eval(L_BICshort);
L_0_BRA=eval(L_BRA);
L_0_BRD=eval(L_BRD);
L_0_ANC=eval(L_ANC);
L_0_PT=eval(L_PT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Dynamic Equation for initial values of q=[0 0 0].
DE_0=simplify(eval(M_q*qddot+C_qqdot*qdot+g_q));