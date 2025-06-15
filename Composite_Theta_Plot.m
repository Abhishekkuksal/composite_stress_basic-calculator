clear all;
clc;
for i_theta=0:1:90

%#######__Properties__#######
Ex= 181e3;
Ey=10.3e3;
Es=7.17e3;
Vx=0.28;
theta=i_theta;

%#######_____Stress At Offaxis______######

stress_1=100;
stress_2=50;
stress_6=20;
stress_offAxis=[stress_1;
                stress_2;
                stress_6];

%######__finding Unknown__#########

Vy=Ey*Vx/Ex;
m=(1-Vx*Vy)^-1;

Qxx=m*Ex
Qxy=m*Vx*Ey;
Qyy=m*Ey;
Qss=Es

Q_onAxis=[Qxx Qxy 0;
          Qxy Qyy 0;
          0   0  Qss]
%###__compliance matrix__type

%######__Invariance Terms__#######
U1=1/8*(3*Qxx+3*Qyy+2*Qxy+4*Qss);
U2=1/2*(Qxx-Qyy);
U3=1/8*(Qxx+Qyy-2*Qxy-4*Qss);
U4=1/8*(Qxx+Qyy+6*Qxy-4*Qss);
U5=1/8*(Qxx+Qyy-2*Qxy+4*Qss);

%######__S__#########
% v_1=2/h*cos[2*theta[groupNo]]*H[groupNo];
% v_2=2/h*cos[4*theta[groupNo]]*H[groupNo];
% v_3=2/h*sin[2*theta[groupNo]]*H[groupNo];
% v_4=2/h*sin[4*theta[groupNo]]*H[groupNo];


%######__Off-Axis Stiffness Matrix__#########
Q11=U1+U2*cosd(2*theta)+U3*cosd(4*theta);
Q12=U4-U3*cosd(4*theta);
Q16=1/2*(U2*sind(2*theta))+U3*sind(4*theta);
Q66=U5-U3*cosd(4*theta);
Q22=U1-U2*cosd(2*theta)+U3*cosd(4*theta);
Q26=1/2*(U2*sind(2*theta))-U3*sind(4*theta);

Q_offAxis=[Q11 Q12 Q16;
           Q12 Q22 Q26;
           Q16 Q26 Q66]

%######__off-axis compliance__#########

det(Q_offAxis)

%##S11=Q22*Q66-Q26^2/delta;
% ##S12=Q16*Q26-Q12*Q66/delta;
% ##S16=Q12*Q26-Q22*Q16/delta;
% ##S22=Q11*Q66-Q16^2/delta;
% ##S26=Q12*Q16-Q11*Q26/delta;
% ##S66=Q11*Q12-Q12^2/delta;

% ##S=[S11 S12 S16;
% ##   S12 S22 S26;
% ##   S16 S26 S66]

S=inv(Q_offAxis)

% ######__off-axis strain__#########

strain_offAxis=[S] * [stress_offAxis]

% ######__graph plot__#########
title('Compliance(S) vs \theta');
xlabel('\Theta');
ylabel('Compliance');
plot(i_theta,S(1,1),"_r",i_theta,S(1,2),"_g",i_theta,S(1,3),"+r",i_theta,S(2,2),"_b",i_theta,S(2,3),"_k",i_theta,S(3,3),"_m");
legend({'S11','S12','S16','S22','S26','S66'},'Location','northeast');
% ##subplot(i_theta,S12,'-');
% ##plot(i_theta,S16,'-');
% ##plot(i_theta,S22,'-');
% ##plot(i_theta,S26,'-');
% ##plot(i_theta,S66,'-');
hold on;

end
