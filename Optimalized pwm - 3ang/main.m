clear;
options = optimoptions('fmincon');
options = optimoptions(options,'Display', 'off');
options = optimoptions(options,'Algorithm', 'sqp');


Modulation = 0:0.001:4/pi;
THD = zeros(0,length(Modulation ));
THDtemp  = THD;

alpha1 = THD;
alpha11 = THD;

alpha2 = THD;
alpha22 = THD;

alpha3 = THD;
alpha33 = THD;
alpha_Fold = zeros(length(Modulation ),3);
global M;
M = 0.01;

alpha00(1,1) = 60*pi/180;
alpha00(1,2)= 79*pi/180;
alpha00(1,3) = 81*pi/180;

alpha00(1,4) = 1*pi/180;
alpha00(1,5)= 61*pi/180;
alpha00(1,6) = 89*pi/180;

alpha00(1,7) = 0*pi/180;
alpha00(1,8)=  1*pi/180;
alpha00(1,9) = 61*pi/180;


lb = [0,0,0];
ub = [pi/2,pi/2,pi/2];
A = [1,-1,0;
     0,1,-1];
b = [0;0];

for ii = 1:length(Modulation)
M = Modulation(ii);
if(Modulation(ii)<=0.01)
temp(1) = alpha00(1,1);
temp(2) = alpha00(1,2);
temp(3) = alpha00(1,3);
else
temp(1) = alpha00(ii,1);
temp(2) = alpha00(ii,2);
temp(3) = alpha00(ii,3);
end
[x,fval] = fmincon(@WTHD_A,temp,A,b,[],[],lb,ub,@ncf_A,options);
if(Modulation(ii)<=0.01)
THDtemp(1,ii) = 1;
else
THDtemp(1,ii) = fval;
end

alpha00(ii+1,1) = x(1);
alpha00(ii+1,2) = x(2);
alpha00(ii+1,3) = x(3);
alpha11(1,ii) = x(1)*180/pi;
alpha22(1,ii) = x(2)*180/pi;
alpha33(1,ii) = x(3)*180/pi;
end

for ii = 1:length(Modulation)
M = Modulation(ii);
if(Modulation(ii)<=0.01)
temp(1) = alpha00(1,4);
temp(2) = alpha00(1,5);
temp(3) = alpha00(1,6);
else
temp(1) = alpha00(ii,4);
temp(2) = alpha00(ii,5);
temp(3) = alpha00(ii,6);
end
[x,fval] = fmincon(@WTHD_B,temp,[],[],[],[],lb,ub,@ncf_B,options);
if(Modulation(ii)<=0.01)
THDtemp(2,ii) = 1;
else
THDtemp(2,ii) = fval;
end
alpha00(ii+1,4) = x(1);
alpha00(ii+1,5) = x(2);
alpha00(ii+1,6) = x(3);
alpha11(2,ii) = x(1)*180/pi;
alpha22(2,ii) = x(2)*180/pi;
alpha33(2,ii) = x(3)*180/pi;
end

for ii = 1:length(Modulation)
M = Modulation(ii);
if(Modulation(ii)<=0.01)
temp(1) = alpha00(1,7);
temp(2) = alpha00(1,8);
temp(3) = alpha00(1,9);
else
temp(1) = alpha00(ii,7);
temp(2) = alpha00(ii,8);
temp(3) = alpha00(ii,9);
end
[x,fval] = fmincon(@WTHD_A,temp,[],[],[],[],lb,ub,@ncf_A,options);
if(Modulation(ii)<=0.01)
THDtemp(3,ii) = 1;
else
THDtemp(3,ii) = fval;
end
alpha00(ii+1,7) = x(1);
alpha00(ii+1,8) = x(2);
alpha00(ii+1,9) = x(3);
alpha11(3,ii) = x(1)*180/pi;
alpha22(3,ii) = x(2)*180/pi;
alpha33(3,ii) = x(3)*180/pi;
end


for ii = 1:length(Modulation)
    if(Modulation(ii)<=0.5)
    alpha1(ii) = alpha11(1,ii);
    alpha2(ii) = alpha22(1,ii);
    alpha3(ii) = alpha33(1,ii);
    THD(ii) = THDtemp(1,ii);
    
    
    elseif(flag == 1)
       alpha1(ii) = alpha11(3,ii);
    alpha2(ii) = alpha22(3,ii);
    alpha3(ii) = alpha33(3,ii);
    THD(ii) = THDtemp(3,ii); 
    
    elseif (THDtemp(1,ii)<=THDtemp(2,ii)&&THDtemp(1,ii)<=THDtemp(3,ii))
    alpha1(ii) = alpha11(1,ii);
    alpha2(ii) = alpha22(1,ii);
    alpha3(ii) = alpha33(1,ii);
    THD(ii) = THDtemp(1,ii);

    elseif (THDtemp(3,ii)<=THDtemp(2,ii)&&THDtemp(3,ii)<=THDtemp(1,ii))
    alpha1(ii) = alpha11(3,ii);
    alpha2(ii) = alpha22(3,ii);
    alpha3(ii) = alpha33(3,ii);
    THD(ii) = THDtemp(3,ii);   
   flag = 1;
    elseif (THDtemp(2,ii)<=THDtemp(1,ii)&&THDtemp(2,ii)<=THDtemp(3,ii))
    alpha1(ii) = alpha11(2,ii);
    alpha2(ii) = alpha22(2,ii);
    alpha3(ii) = alpha33(2,ii);
    THD(ii) = THDtemp(2,ii);   
    
    end
end
subplot(2,1,1);
plot(Modulation,THD);
xlim([0 1.27]);
ylim([0 0.16]);
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
xlim([0 1.27]);
ylim([0 90]);
xlabel('Modulation Index(M)');
ylabel('Swiching Angles(deg)');
title('开关角度')






for ii = 1:length(Modulation)
    alpha(ii,1) = alpha1(ii);
    alpha(ii,2) = alpha2(ii);
    alpha(ii,3) = alpha3(ii);
end




% clear alpha00 alpha00_A alpha00_B flag alpha11_A A alpha11 alpha22 alpha33 b  THDtemp temp alpha11_B alpha1_A alpha1_B fval ii lb M Modulation options THDA THDB ub x





