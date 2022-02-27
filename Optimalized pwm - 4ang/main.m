clear;
options = optimoptions('fmincon');
options = optimoptions(options,'Display', 'off');
options = optimoptions(options,'Algorithm', 'sqp');


Modulation = 0:0.001:4/pi;

THD = zeros(0,length(Modulation ));
THDtemp = zeros(0,length(Modulation ));

alpha1 = zeros(0,length(Modulation ));
alpha11 = zeros(0,length(Modulation ));

alpha2 = zeros(0,length(Modulation ));
alpha22 = zeros(0,length(Modulation ));

alpha3 = zeros(0,length(Modulation ));
alpha33 = zeros(0,length(Modulation ));

alpha4 = zeros(0,length(Modulation ));
alpha44 = zeros(0,length(Modulation ));



global M;
M = 0.01;

alpha00(1,1) = 60.0747*pi/180;%A
alpha00(1,2) = 75.2968*pi/180;
alpha00(1,3) = 75.4056*pi/180;
alpha00(1,4) = 89.9*pi/180;

alpha00(1,5) = 1*pi/180;%B
alpha00(1,6) = 60*pi/180;
alpha00(1,7) = 79*pi/180;
alpha00(1,8) = 81*pi/180;

alpha00(1,9) = 0.84743*pi/180;
alpha00(1,10) = 1.7697*pi/180;
alpha00(1,11) = 63.0388*pi/180;
alpha00(1,12) = 84.4739*pi/180;

alpha00(1,13) = 0.56477*pi/180;
alpha00(1,14) = 1.1352*pi/180;
alpha00(1,15) = 1.6188*pi/180;
alpha00(1,16) = 66*pi/180;

alpha00(1,17) = 7.9646*pi/180;
alpha00(1,18) = 13.848*pi/180;
alpha00(1,19) = 22.7392*pi/180;
alpha00(1,20) = 40.9001*pi/180;
lb = [0,0,0,0];
ub = [pi/2,pi/2,pi/2,pi/2];
A = [1,-1,0,0;
     0,1,-1,0;
     0,0,1,-1];
b = [0;0;0];
for ii = 1:length(Modulation)

 M = Modulation(ii);
if(Modulation(ii)<=0.05)
temp(1) = alpha00(1,1);
temp(2) = alpha00(1,2);
temp(3) = alpha00(1,3);
temp(4) = alpha00(1,4);
else
temp(1) = alpha00(ii,1);
temp(2) = alpha00(ii,2);
temp(3) = alpha00(ii,3);
temp(4) = alpha00(ii,4);
end
[x,fval] = fmincon(@WTHD_A,temp,A,b,[],[],lb,ub,@ncf_A,options);
THDtemp(1,ii) = fval;
alpha00(ii+1,1) = x(1);
alpha00(ii+1,2) = x(2);
alpha00(ii+1,3) = x(3);
alpha00(ii+1,4) = x(4);
alpha11(1,ii) = x(1)*180/pi;
alpha22(1,ii) = x(2)*180/pi;
alpha33(1,ii) = x(3)*180/pi;
alpha44(1,ii) = x(4)*180/pi;
end

for ii = 1:length(Modulation)
M = Modulation(ii);
if(Modulation(ii)<=0.05)
temp(1) = alpha00(1,5);
temp(2) = alpha00(1,6);
temp(3) = alpha00(1,7);
temp(4) = alpha00(1,8);
else
temp(1) = alpha00(ii,5);
temp(2) = alpha00(ii,6);
temp(3) = alpha00(ii,7);
temp(4) = alpha00(ii,8);
end
[x,fval] = fmincon(@WTHD_B,temp,A,b,[],[],lb,ub,@ncf_B,options);
THDtemp(2,ii) = fval;
alpha00(ii+1,5) = x(1);
alpha00(ii+1,6) = x(2);
alpha00(ii+1,7) = x(3);
alpha00(ii+1,8) = x(4);
alpha11(2,ii) = x(1)*180/pi;
alpha22(2,ii) = x(2)*180/pi;
alpha33(2,ii) = x(3)*180/pi;
alpha44(2,ii) = x(4)*180/pi;
end

for ii = 1:length(Modulation)
M = Modulation(ii);
if(Modulation(ii)<=0.22)
temp(1) = alpha00(1,9);
temp(2) = alpha00(1,10);
temp(3) = alpha00(1,11);
temp(4) = alpha00(1,12);
else
temp(1) = alpha00(ii,9);
temp(2) = alpha00(ii,10);
temp(3) = alpha00(ii,11);
temp(4) = alpha00(ii,12);
end
[x,fval] = fmincon(@WTHD_A,temp,A,b,[],[],lb,ub,@ncf_A,options);

if(Modulation(ii)<=0.22)
THDtemp(3,ii) = 1;
else
THDtemp(3,ii) = fval;
end
alpha00(ii+1,9) = x(1);
alpha00(ii+1,10) = x(2);
alpha00(ii+1,11) = x(3);
alpha00(ii+1,12) = x(4);
alpha11(3,ii) = x(1)*180/pi;
alpha22(3,ii) = x(2)*180/pi;
alpha33(3,ii) = x(3)*180/pi;
alpha44(3,ii) = x(4)*180/pi;
end

for ii = 1:length(Modulation)
M = Modulation(ii);
if(Modulation(ii)<=0.17)
temp(1) = alpha00(1,13);
temp(2) = alpha00(1,14);
temp(3) = alpha00(1,15);
temp(4) = alpha00(1,16);
else
temp(1) = alpha00(ii,13);
temp(2) = alpha00(ii,14);
temp(3) = alpha00(ii,15);
temp(4) = alpha00(ii,16);
end
[x,fval] = fmincon(@WTHD_B,temp,A,b,[],[],lb,ub,@ncf_B,options);
THDtemp(4,ii) = fval;
alpha00(ii+1,13) = x(1);
alpha00(ii+1,14) = x(2);
alpha00(ii+1,15) = x(3);
alpha00(ii+1,16) = x(4);
alpha11(4,ii) = x(1)*180/pi;
alpha22(4,ii) = x(2)*180/pi;
alpha33(4,ii) = x(3)*180/pi;
alpha44(4,ii) = x(4)*180/pi;
end

for ii = 1:length(Modulation)
% M = Modulation(ii);
% if(Modulation(ii)<=0.8)
% temp(1) = alpha00(1,17);
% temp(2) = alpha00(1,18);
% temp(3) = alpha00(1,19);
% temp(4) = alpha00(1,20);
% else
% temp(1) = alpha00(ii,17);
% temp(2) = alpha00(ii,18);
% temp(3) = alpha00(ii,19);
% temp(4) = alpha00(ii,20);
% end
% [x,fval] = fmincon(@WTHD_A,temp,A,b,[],[],lb,ub,@ncf_A,options);
% if(Modulation(ii)<=0.8)
THDtemp(5,ii) = 1;
% else
% THDtemp(5,ii) = fval;
% end
% alpha00(ii+1,17) = x(1);
% alpha00(ii+1,18) = x(2);
% alpha00(ii+1,19) = x(3);
% alpha00(ii+1,20) = x(4);
% alpha11(5,ii) = x(1)*180/pi;
% alpha22(5,ii) = x(2)*180/pi;
% alpha33(5,ii) = x(3)*180/pi;
% alpha44(5,ii) = x(4)*180/pi;
end

for ii = 1:length(Modulation)
    if(Modulation(ii)<=0.1)
    alpha1(ii) = alpha11(1,ii);
    alpha2(ii) = alpha22(1,ii);
    alpha3(ii) = alpha33(1,ii);
    alpha4(ii) = alpha44(1,ii);
    THD(ii) = THDtemp(1,ii);
    
    elseif(flag == 1)
    alpha1(ii) = alpha11(5,ii);
    alpha2(ii) = alpha22(5,ii);
    alpha3(ii) = alpha33(5,ii);
    alpha4(ii) = alpha44(5,ii);
    THD(ii) = THDtemp(5,ii);
    
    elseif (THDtemp(1,ii)<=THDtemp(2,ii)&&THDtemp(1,ii)<=THDtemp(3,ii)&&THDtemp(1,ii)<=THDtemp(4,ii)&&THDtemp(1,ii)<=THDtemp(5,ii))
    alpha1(ii) = alpha11(1,ii);
    alpha2(ii) = alpha22(1,ii);
    alpha3(ii) = alpha33(1,ii);
    alpha4(ii) = alpha44(1,ii);
    THD(ii) = THDtemp(1,ii);

    elseif (THDtemp(2,ii)<=THDtemp(1,ii)&&THDtemp(2,ii)<=THDtemp(3,ii)&&THDtemp(2,ii)<=THDtemp(4,ii)&&THDtemp(2,ii)<=THDtemp(5,ii))
    alpha1(ii) = alpha11(2,ii);
    alpha2(ii) = alpha22(2,ii);
    alpha3(ii) = alpha33(2,ii);
    alpha4(ii) = alpha44(2,ii);
    THD(ii) = THDtemp(2,ii);
    
    
    elseif (THDtemp(3,ii)<=THDtemp(1,ii)&&THDtemp(3,ii)<=THDtemp(2,ii)&&THDtemp(3,ii)<=THDtemp(4,ii)&&THDtemp(3,ii)<=THDtemp(5,ii))
    alpha1(ii) = alpha11(3,ii);
    alpha2(ii) = alpha22(3,ii);
    alpha3(ii) = alpha33(3,ii);
    alpha4(ii) = alpha44(3,ii);
    THD(ii) = THDtemp(3,ii);   
   
    elseif (THDtemp(4,ii)<=THDtemp(1,ii)&&THDtemp(4,ii)<=THDtemp(2,ii)&&THDtemp(4,ii)<=THDtemp(3,ii)&&THDtemp(4,ii)<=THDtemp(5,ii))
    alpha1(ii) = alpha11(4,ii);
    alpha2(ii) = alpha22(4,ii);
    alpha3(ii) = alpha33(4,ii);
    alpha4(ii) = alpha44(4,ii);
    THD(ii) = THDtemp(4,ii);   
    
%     elseif (THDtemp(5,ii)<=THDtemp(1,ii)&&THDtemp(5,ii)<=THDtemp(2,ii)&&THDtemp(5,ii)<=THDtemp(3,ii)&&THDtemp(5,ii)<=THDtemp(4,ii))
%     alpha1(ii) = alpha11(5,ii);
%     alpha2(ii) = alpha22(5,ii);
%     alpha3(ii) = alpha33(5,ii);
%     alpha4(ii) = alpha44(5,ii);
%     THD(ii) = THDtemp(5,ii);  
%     flag = 1;%折角标志位
    
    end
end

subplot(2,1,1);
plot(Modulation,THD);
xlim([0 1.27]);
ylim([0 0.11]);
xlabel('Modulation Index(M)');
ylabel('WTHD');
title('谐波畸变率')

subplot(2,1,2);
plot(Modulation,alpha1);
hold on;
plot(Modulation,alpha2);
hold on;
plot(Modulation,alpha3);
hold on;
plot(Modulation,alpha4);
hold on;
xlim([0 1.27]);
ylim([0 90]);
xlabel('Modulation Index(M)');
ylabel('Swiching Angles(deg)');
title('开关角度')
for ii = 1:length(Modulation)
    alpha(ii,1) = alpha1(ii);
    alpha(ii,2) = alpha2(ii);
    alpha(ii,3) = alpha3(ii);
    alpha(ii,4) = alpha4(ii);
end
clear alpha00 alpha00_A alpha00_B alpha11_A A alpha11 alpha22 alpha33 alpha44 b  THDtemp temp alpha11_B alpha1_A alpha1_B fval ii lb M Modulation options THDA THDB ub x







