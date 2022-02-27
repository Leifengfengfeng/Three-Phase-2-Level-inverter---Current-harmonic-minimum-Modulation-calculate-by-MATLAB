clear;
options = optimoptions('fmincon');
options = optimoptions(options,'Display', 'off');
options = optimoptions(options,'Algorithm', 'sqp');
lb = [0,0];
ub = [pi/2+0.1,pi/2+0.1];

global M;
M = 0.01;

Modulation = 0:0.001:1.27;
THD = zeros(0,length(Modulation));
alpha_fold = zeros(length(Modulation),2);
alpha00 = zeros(length(Modulation) + 1,2);

alpha00(1,1) = 0.11111*pi/180;
alpha00(1,2) =  59.998*pi/180;


for ii = 1:length(Modulation)
M = Modulation(ii);
if(Modulation(ii)<=0.01)
temp(1) = alpha00(1,1);
temp(2) = alpha00(1,2);
else
temp(1) = alpha00(ii,1);
temp(2) = alpha00(ii,2);   
end
[x,fval] = fmincon(@WTHD_A,temp,[],[],[],[],lb,ub,@ncf_A,options);
if(Modulation(ii)<=0.01)
THD(ii) = 1;
else
THD(ii) = fval;  
end
alpha00(ii+1,1) = x(1);
alpha00(ii+1,2) = x(2);
alpha_fold(ii,1) = x(1)*180/pi;
alpha_fold(ii,2) = x(2)*180/pi;
end

subplot(2,1,1);
plot(Modulation,THD);
xlim([0 1.27]);
ylim([0 0.20]);
xlabel('Modulation Index(M)');
ylabel('WTHD');
title('谐波畸变率')
subplot(2,1,2);
plot(Modulation,alpha_fold);
hold on;
xlim([0 1.27]);
ylim([0 90]);
xlabel('Modulation Index(M)');
ylabel('Swiching Angles(deg)');
title('开关角度')