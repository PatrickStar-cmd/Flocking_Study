clear all; close all;
fs = 100;
T = 20;
dt = 1/fs;
N = T * fs;
t = (0 : N-1) * dt;

%% ture angle
true_angle = 30 * sin(2*pi * 0.2*t) * pi / 180;
%% true angle velocity
true_Psai = [0 diff(true_angle)/dt];

%% gyro
gyro_bias = 0.5 * pi / 180;
gyro_noise_std = 0.2 * pi / 180;
Psai = true_Psai + gyro_bias + gyro_noise_std * randn(1, N);

%% acc
acc_noise_std = 0.8 * pi / 180;
acc_angle = true_angle + acc_noise_std * randn(1, N);

%% 互补滤波
angle_cf = zeros(1, N);
angle_cf(1) = acc_angle(1);
alpha = 0.98;
for k = 2 : N
    tmp1=angle_cf(k-1)+Psai(k)*dt;
    angle_cf(k)=alpha*tmp1+(1-alpha)*acc_angle(k);
end

%%acc
angle_acc = zeros(1, N);
angle_acc(1) = acc_angle(1);
for k=2:N
    angle_acc(k)=0.75*angle_acc(k-1)+0.25*acc_angle(k);
end

%%gyro inter
angle_gyro=zeros(1,N);
angle_gyro(1)=acc_angle(1);
for k=2:N
    angle_gyro(k)=angle_gyro(k-1)+Psai(k)*dt;

end

%%kf
xk=[acc_angle(1);0];
Pk=eye(2);
q_theta=acc_noise_std^2;
q_b=(0.01*pi/180)^2;
Qk=diag([q_theta*dt;q_b*dt]);
Rk=acc_noise_std^2;
Fk=[1-dt;0 1];
Bk=[dt;0];
theta_kf=zero(1,N);
bias_kf=zeros(1,N);
for k=1:N
    xk1=Fk*xk+Bk*Psai(k);
    Pk1=Fk*Pk*Fk'+Qk;
    Gk=Pk1*Hk'*inv(Hk*Pk1*Hk'+Pk);
    yk=acc_angle(k);
    xk=xk1+Gk*(yk-Hk*xk1);
    Pk=(eye(2)-Gk*Hk)*Pk1;
    theta_kf(k)=xk(1);
    bias_kf(k)=xk(2);
end




figure,hold on;
plot(t,true_angle,'0');
plot(t,angle_cf,'*');
plot(t,angle_acc,'*');
legend('true','cf','acc');

err_cf=true_angle-angle_cf;
err_acc=true_angle-angle_acc;
err_gyro=true_angle-angle_gyro
err_kf=true_angle-theta_kf;

[rms(err_cf) rms(err_acc) rms(err_gyro)]*180/pi



