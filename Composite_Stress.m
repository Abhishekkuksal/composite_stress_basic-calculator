%*this program is askes group and

clear all;
clc;
%disp("The program is created to find to find strain in a multigroup laminate /n The below program is strickly made for SYMMETRIC laminate")
%disp("X is the PLY ANGLE in degree, range is from -90 to +90");
%disp("N is the NUMBER, it should be an integer");
%disp("S represents Symmetricity and it must be present the format as this code is for symmetric laminate only");
%disp("Ply groups are seperated by ""|"", you can have N number of groups in a laminate");
%disp("Below is the generalised format for 2 group");
%disp("[X_N|X_N]S");

inputCommandActive= 1; % 1 for active 0 for inactive
disp("The program is created to find to find strain in a multigroup laminate \nThe below program is strickly made for SYMMETRIC laminate")
halfGroup=input('How many SYMMETRIC lamina (not including repeting lamina) does it has?:-');
plyAngleGroup=zeros(1,halfGroup);
for i=1:halfGroup

  if i~=halfGroup
  capture=input(sprintf('Type the %d group PLY laminate starting from top:- ', i));
  plyAngleGroup(1,i)=capture;

  else
  capture=input(sprintf('Type the last group i.e %d of PLY laminate before SYMMETRY:- ', i));
  plyAngleGroup(1,i)=capture;
  end
end
%M=['This is your Ply Group with SYMMETRY [',plyAngleGroup,']'];
disp('This is your Ply Group with SYMMETRY [ %s ]',plyAngleGroup);
%fprintf('This is your Ply Group with SYMMETRY %i ', plyAngleGroup);
disp(plyAngleGroup);
%disp(sprintf('this is your ply angle config %d  S', plyAngleGroup));

%#######__Properties__#######
if (inputCommandActive);
Ex= input('Enter Ex (on-Axis) Value in GPa:- ');
Ey= input('Enter EY (on-Axis) Value in GPa:- ');
Es= input('Enter Es (on-Axis) Value in GPa:- ');
Vx= input('Enter Vx (on-Axis) Value:- ');
h=input('Enter Height of lamina:- ');
else
Ex= 38.6;
Ey= 8.27;
Es= 4.14;
Vx= 0.26;
h= 125e-6;
end
H=(2*halfGroup)*h %Total Height

%#######_____Stress At Offaxis______######
if (inputCommandActive)
stress_1= input('Enter Stress1 (S1) Value:- ');
stress_2=input('Enter Stress2 (S2) Value:- ');
stress_6=input('Enter Stress6 (S6) Value:- ');
else
stress_1= 100;
stress_2=50;
stress_6=20;
end
stress_offAxis=[stress_1;
                stress_2;
                stress_6];

%######__finding Unknown__#########
Vy=Ey*Vx/Ex;
m=(1-Vx*Vy)^-1;
Qxx=m*Ex;
Qxy=m*Vx*Ey;
Qyy=m*Ey;
Qss=Es;
Q_onAxis=[Qxx Qxy 0;
          Qxy Qyy 0;
          0   0  Qss];
%###__compliance matrix__type

%######__Invariant Terms__#######
U1=1/8*(3*Qxx+3*Qyy+2*Qxy+4*Qss);
U2=1/2*(Qxx-Qyy);
U3=1/8*(Qxx+Qyy-2*Qxy-4*Qss);
U4=1/8*(Qxx+Qyy+6*Qxy-4*Qss);
U5=1/8*(Qxx+Qyy-2*Qxy+4*Qss);

%####__V*___########
v1_sum=0;
v2_sum=0;
v3_sum=0;
v4_sum=0;
for i=1:halfGroup
  v1_current=cosd(2*plyAngleGroup(1,i))*h;
  v1_sum=v1_current+v1_sum;

  v2_current=cosd(4*plyAngleGroup(1,i))*h;
  v2_sum=v2_current+v2_sum;

  v3_current=sind(2*plyAngleGroup(1,i))*h;
  v3_sum=v3_current+v3_sum;

  v4_current=sind(4*plyAngleGroup(1,i))*h;
  v4_sum=v4_current+v4_sum;
end

v1=2*v1_sum/H;
v2=2*v2_sum/H;
v3=2*v3_sum/H;
v4=2*v4_sum/H;

%##%laminate stiffness modulus
A11=(U1+(v1*U2)+(v2*U3))*H;
A12=(U4-(v2*U3))*H;
A16=((v3/2)*U2)*H;
A21=(U4-(v2-U3))*H;
A22=(U1-(v1*U2)+(v2*U3))*H;
A26=((v3/2)-(v4*U3))*H;
A66=(U5-(v2*U3))*H;
A21=A12;
A62=A26;
A61=A16;

AH=[A11/H A12/H A16/H;
   A21/H A22/H A26/H;
   A61/H A62/H A66/H]

   A=[A11 A12 A16;
   A21 A22 A26;
   A61 A62 A66]

%laminate compliance stresses
a=[A]^-1
a11=(A11)^-1
a12=(A12)^-1
a21=(A21)^-1
a22=(A22)^-1
a26=(A26)^-1
a66=(A66)^-1
a16=(A16)^-1
a61=a16;
a62=a26;

%stress resultant N1,N2,N6
N1 = 1;
N2 = 0;
N6 = 0;
N = [N1;N2;N6];
% strain resultant EE1,EE2,EE6
A = [A11 A12 A16;A21 A22 A26;A61 A62 A66];%#remove /H
EE = pinv(A)*N
EE = [a]*[N]
%EE = [EE1; EE2; EE6]
EE1 = EE(1);
EE2 = EE(2);
EE6 = EE(3);
% stress value sigma_1,sigma_2,sigma_6
%sigma = [sigma_1;sigma_2;sigma_6]
sigma = (1/H)*[A]*[EE]

