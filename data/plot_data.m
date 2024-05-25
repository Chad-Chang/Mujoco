clear all
clc
close all
filename = 'data/DOB.csv';

T = readtable(filename); %check T.Properties
VariableNames = T.Properties.VariableNames;

data_length = 13;
Arr = table2array(T(:,1:data_length));
[m,n] = size(Arr);

% disturbance estimation
% figure(1)
% subplot(2,1,1);
% plot(Arr(:,1),Arr(:,2),'b-');
% hold on 
% plot(Arr(:,1),Arr(:,3),'r-');
% hold on 
% plot(Arr(:,1),Arr(:,8),'k-');
% grid on 
% legend("dist", "e-dist", "gravity term"  );
% title("disturbace and estimated disturbance")
% ylim([-10,10]);
% 
% subplot(2,1,2);
% plot(Arr(:,1),Arr(:,4),'k-');
% title("disturbance error" );
% legend("dist err");
% grid on 
% ylim([-10,10]);
% 
figure(2)
subplot(2,1,1);
plot(Arr(:,1),Arr(:,5),'b.');
hold on 
plot(Arr(:,1),Arr(:,6),'r.');
grid on 
legend("input", "output" );
title("input and output");
% 
% subplot(2,1,2);
% plot(Arr(:,1),Arr(:,7),'k.');
% legend("I/O err" );
% legend("Input output error");
% grid on 
% 
% figure(3)
% subplot(2,1,1)
% % plot(T.t, T.ext_force_x_hat,'b');
% plot(Arr(:,1), Arr(:,8), 'r-');
% hold on
% plot(Arr(:,1), Arr(:,9), 'b-');
% % plot(T.t, T.ext_force_y_hat,'r');
% hold on 
% legend("force_x","force_y");
% grid on 
% title("force observer");
% ylim([-1,1]);
% 
% 
% subplot(2,1,2)
% % plot(T.t, T.dist,'r');
% plot(Arr(:,1), Arr(:,2), 'r-');
% hold on 
% % plot(T.t, T.torque_hat,'b-');
% plot(Arr(:,1), Arr(:,10), 'b-');
% hold on 
% plot(Arr(:,1), Arr(:,4), ' k-');
% legend("disturbance", "estimated torque", "disturbance error");
% % plot(T.t, T.dist-T.torque_hat,'-');
% % hold on 
% % plot(T.t, T.dist_err-T.torque_hat,'k-');
% % legend("torque(disturbance)", "estimated torque","est err")
% grid on
% title("torque estimation")
% ylim([-10,10]);

% figure(4)
% subplot(2,1,1)
% plot(Arr(:,1), Arr(:,6), 'r-');
% hold on 
% plot(Arr(:,1), Arr(:,13), 'b-');
% hold on 
% plot(Arr(:,1), Arr(:,12), 'k-');
% hold on
% grid on 
% legend("q_pos", "delta_x", "reference");
% title("admittance control")





