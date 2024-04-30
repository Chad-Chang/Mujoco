clear all
clc
close all
filename = 'data/DOB.csv';

T = readtable(filename); %check T.Properties
VariableNames = T.Properties.VariableNames;

Arr = table2array(T(:,1:8));
[m,n] = size(Arr);

% disturbance estimation
figure(1)
subplot(2,1,1);
plot(Arr(:,1),Arr(:,2),'b-');
hold on 
plot(Arr(:,1),Arr(:,3),'r-');
hold on 
plot(Arr(:,1),Arr(:,8),'k-');
grid on 
legend("dist", "e-dist", "gravity term"  );
title("disturbace and estimated disturbance")
ylim([-10,10]);

subplot(2,1,2);
plot(Arr(:,1),Arr(:,4),'k-');
title("disturbance error" );
legend("dist err");
grid on 
ylim([-10,10]);

figure(2)
subplot(2,1,1);
plot(Arr(:,1),Arr(:,5),'b.');
hold on 
plot(Arr(:,1),Arr(:,6),'r.');
grid on 
legend("input", "output" );
title("input and output");

subplot(2,1,2);
plot(Arr(:,1),Arr(:,7),'k.');
legend("I/O err" );
legend("Input output error");
grid on 


% 
% for i=2:m
%     figure(i)
%     yy = i;
%     plot(Arr(:,yy),'r');
%     ylabel(cell2mat(VariableNames(yy)))
% 
% end