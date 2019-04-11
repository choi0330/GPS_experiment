clear all;
close all;
clc;

load(['EXP_20171011_test4_data']) 

%% time setting
ini_epoch = 5*60;              % 초기 미지정수를 위한 정지구간
n_epoch = size(true.GPSTime,1); % 총 데이터 길이
maskangle = 15;                 
nt = 32;

%% Visible satellite selection
PRN = [1:32];
vis_sat = zeros(1,32);
for ii = 1:32
    
    % mask angle보다 elevation angle이 높고, 측정치가 있는 위성 찾기
    temp_el1 = find(ref.El(:,ii) < maskangle| user.El(:,ii) < maskangle); 
    temp_el2 = find(ref.pr1(:,ii)<10 |user.pr1(:,ii)<10);
    if isempty(temp_el1)
        if isempty(temp_el2)
            vis_sat(ii)=1;
        end
    end
end
SV_vis = find(vis_sat==1);
n_vis = length(SV_vis);

%% Standalone 
 user.pos_standalone = zeros(n_epoch,3);
 user.pos_standalone_error = zeros(n_epoch,3);
 user.pos_standalone_enu = zeros(n_epoch,3);
 R_ur = zeros(n_epoch,3);
 
 for ti = 1: n_epoch
     
     dx = 100;
     % Step 1. 사용자 위치 원점 가정
     R_user = zeros(3,1).'; % 초기 위치 가정
     B_ur = 0;              % 수신기 시계오차
     x_old = [R_user.';B_ur];
     
     while(dx> 10^-4)
         
         H_standalone = zeros(n_vis,4);
         z_standalone =  zeros(n_vis,1);

         for jj = 1: n_vis
             
             % Step 2. e_hat 계산
             R_su = user.svpos(ti,3*SV_vis(jj)-2:3*SV_vis(jj)) - R_user;
             e_hat = R_su/norm(R_su);
             
             % Step 3. 사용자 위치 계산 하기
             % Step 3.1. H matrix 계산
             H_standalone(jj,1:3) = e_hat;
             H_standalone(jj,4) = -1;

             % Step 3.2. measurement vector 계산
             z_standalone(jj,1) = e_hat*user.svpos(ti,3*SV_vis(jj) ...
             -2:3*SV_vis(jj)).'-(user.pr1(ti,SV_vis(jj))+user.b(ti,SV_vis(jj))*c);
         end
         
         % Step 3.3. user position 계산
         x = pinv(H_standalone)*z_standalone;
         R_user = x(1:3).';
         B_ur = x(4);
         
         % Step 4. 사용자 위치 비교 후 update
         dx = norm(x-x_old);
         x_old = x;
         
     end
     
     user.pos_standalone(ti,:) = R_user;
     user.pos_standalone_enu(ti,:) = (Rtran*(user.pos_standalone(ti,:) ...
                                     -ref_xyz).').';
     error.standalone(ti,:) = user.pos_standalone(ti,:)  - true.xyz(ti,:);
     
 end
 
 %% Standalone plot

figure(1);
hold on;
plot(true.enu(ini_epoch:n_epoch,1),true.enu(ini_epoch:n_epoch,2),'r.');
plot(user.pos_standalone_enu(ini_epoch:n_epoch,1), ...
     user.pos_standalone_enu(ini_epoch:n_epoch,2),'k.','Marker' ...
     ,'o','MarkerSize',4);
legend('True','Standalone');
title('Standalone Horizontal Trajectory');
xlabel('East(m)');
ylabel('North(m)');
hold off;
figure(2);
hold on;
plot(true.GPSTime(ini_epoch:n_epoch,1),true.enu(ini_epoch:n_epoch,3),'r.');
plot(user.GPSTime(ini_epoch:n_epoch,1),user.pos_standalone_enu( ...
     ini_epoch:n_epoch,3),'k.','Marker','o','MarkerSize',4);
legend('True','Standalone');
title('Standalone Vertical Trajectory');
xlabel('time(s)');
ylabel('Up(m)');
hold off;

%% RMS
add = 0;
add_temp = 0;

for ti=ini_epoch:n_epoch
    add_temp = (norm(error.standalone(ti,:)))^2;
    add = add+add_temp;
end
mean = add/(n_epoch-ini_epoch+1);
rms = sqrt(mean);