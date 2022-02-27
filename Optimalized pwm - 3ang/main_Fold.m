clear;
options = optimoptions('fmincon');
options = optimoptions(options,'Display', 'off');
options = optimoptions(options,'Algorithm', 'sqp');
lb = [0,0,0];
ub = [pi/2,pi/2,pi/2];
A = [1,-1,0;
     0,1,-1];
b = [0;0];

Modulation = 0:0.001:1.27;
THD = zeros(0,length(Modulation ));
alpha_Fold = zeros(length(Modulation ),3);
alpha00 = alpha_Fold;
global M;
M = 0.01;
alpha00(1,1) = 0.56477*pi/180;
alpha00(1,2)= 1.1352*pi/180;
alpha00(1,3) = 58*pi/180;


for ii = 1:length(Modulation)
M = Modulation(ii);
if(Modulation(ii)<=0.54)
temp(1) = alpha00(1,1);
temp(2) = alpha00(1,2);
temp(3) = alpha00(1,3);
else
temp(1) = alpha00(ii,1);
temp(2) = alpha00(ii,2);
temp(3) = alpha00(ii,3);
end
[x,fval] = fmincon(@WTHD_B,temp,[],[],[],[],lb,ub,@ncf_B,options);
if(Modulation(ii)<=0.54)
THD(ii) = 1;
else
THD(ii) = fval;    
end
alpha00(ii+1,1) = x(1);
alpha00(ii+1,2) = x(2);
alpha00(ii+1,3) = x(3);
alpha_Fold(ii,1) = x(1)*180/pi;
alpha_Fold(ii,2) = x(2)*180/pi;
alpha_Fold(ii,3) = x(3)*180/pi;
end

subplot(2,1,1);
plot(Modulation,THD);
xlim([0 1.27]);
ylim([0 0.16]);
xlabel('Modulation Index(M)');
ylabel('WTHD');
title('谐波畸变率')
subplot(2,1,2);
plot(Modulation,alpha_Fold);
hold on;
xlim([0 1.27]);
ylim([0 90]);
xlabel('Modulation Index(M)');
ylabel('Swiching Angles(deg)');
title('开关角度')
