clear;
options = optimoptions('fmincon');
options = optimoptions(options,'Display', 'off');
options = optimoptions(options,'Algorithm', 'sqp');
Modulation = 0:0.001:4/pi;
THD = zeros(0,length(Modulation ));
THDA = zeros(0,length(Modulation ));
alpha1 = zeros(0,length(Modulation ));
alpha11 = zeros(0,length(Modulation ));
alpha2 = zeros(0,length(Modulation ));
alpha22 = zeros(0,length(Modulation ));
global M;
M = 0.01;

alpha00_A(1,1) = 61*pi/180;
alpha00_A(1,2) = 89*pi/180;

alpha00_A(1,3) = 0.11111*pi/180;
alpha00_A(1,4) =  59.998*pi/180;

alpha00_B(1,1) = 0.1111*pi/180;
alpha00_B(1,2) =  60.1111*pi/180;
lb = [0,0];
ub = [pi/2+0.1,pi/2+0.1];
for ii = 1:length(Modulation)
M = Modulation(ii);
if(Modulation(ii)<=0.01)
temp(1) = alpha00_A(1,1);
temp(2) = alpha00_A(1,2);
else
temp(1) = alpha00_A(ii,1);
temp(2) = alpha00_A(ii,2);   
end
[x,fval] = fmincon(@WTHD_A,temp,[],[],[],[],lb,ub,@ncf_A,options);
THDA(1,ii) = fval;
if(Modulation(ii)<=0.01)
THDA(1,ii) = 1;
else
THDA(1,ii) = fval;  
end
alpha00_A(ii+1,1) = x(1);
alpha00_A(ii+1,2) = x(2);
alpha11(1,ii) = x(1)*180/pi;
alpha22(1,ii) = x(2)*180/pi;
end

for ii = 1:length(Modulation)
M = Modulation(ii);
if(Modulation(ii)<=0.2)
temp(1) = alpha00_A(1,3);
temp(2) = alpha00_A(1,4);
else
temp(1) = alpha00_A(ii,3);
temp(2) = alpha00_A(ii,4);   
end
[x,fval] = fmincon(@WTHD_A,temp,[],[],[],[],lb,ub,@ncf_A,options);
if(Modulation(ii)<=0.01)
THDA(2,ii) = 1;
else
THDA(2,ii) = fval;  
end
alpha00_A(ii+1,3) = x(1);
alpha00_A(ii+1,4) = x(2);
alpha11(2,ii) = x(1)*180/pi;
alpha22(2,ii) = x(2)*180/pi;
end

for ii = 1:length(Modulation)
M = Modulation(ii);
if(Modulation(ii)<=0.2)
temp(1) = alpha00_B(1,1);
temp(2) = alpha00_B(1,2);
else
temp(1) = alpha00_B(ii,1);
temp(2) = alpha00_B(ii,2);   
end
[x,fval] = fmincon(@WTHD_B,temp,[],[],[],[],lb,ub,@ncf_B,options);
if(Modulation(ii)<=0.01)
THDA(3,ii) = 1;
else
THDA(3,ii) = fval;  
end
alpha00_B(ii+1,1) = x(1);
alpha00_B(ii+1,2) = x(2);
alpha11(3,ii) = x(1)*180/pi;
alpha22(3,ii) = x(2)*180/pi;
end

for ii = 1:length(Modulation)
    if(Modulation(ii)<=0.5)
    alpha1(ii) = alpha11(1,ii);
    alpha2(ii) = alpha22(1,ii);
    THD(ii) = THDA(1,ii);
    

    elseif (THDA(1,ii)<=THDA(3,ii))
    alpha1(ii) = alpha11(1,ii);
    alpha2(ii) = alpha22(1,ii);
    THD(ii) = THDA(1,ii);

 
    elseif (THDA(3,ii)<=THDA(1,ii))
    alpha1(ii) = alpha11(3,ii);
    alpha2(ii) = alpha22(3,ii);
    THD(ii) = THDA(3,ii);   
     
    end
end

subplot(2,1,1);
plot(Modulation,THD);
xlim([0 1.27]);
ylim([0 0.20]);
xlabel('Modulation Index(M)');
ylabel('WTHD');
title('谐波畸变率')
% subplot(2,2,2);
% plot(Modulation,THDB);
% xlim([0 1.27]);
% ylim([0 0.35]);


subplot(2,1,2);
plot(Modulation,alpha1);
hold on;
plot(Modulation,alpha2);
hold on;
xlim([0 1.27]);
ylim([0 90]);
xlabel('Modulation Index(M)');
ylabel('Swiching Angles(deg)');
title('开关角度')
% subplot(2,2,4);
% plot(Modulation,alpha1_B);
% hold on;
% xlim([0 1.27]);
% ylim([0 90]);

for ii = 1:length(Modulation)
    alpha(ii,1) = alpha1(ii);
    alpha(ii,2) = alpha2(ii);
end

% clear alpha00 alpha00_A alpha00_B alpha11_A  alpha11 alpha22 temp alpha11_B alpha1_A alpha1_B fval ii lb M Modulation options THDA THDB ub x









