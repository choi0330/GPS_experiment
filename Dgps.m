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

%% DGPS 
 user.pos_dgps = zeros(n_epoch,3);
 user.pos_dgps_error = zeros(n_epoch,3);
 user.pos_dgps_enu = zeros(n_epoch,3);
 R_ur = zeros(n_epoch,3);
 HDOP = zeros(1,n_epoch);
 VDOP = zeros(1,n_epoch);
 PDOP = zeros(1,n_epoch);
 TDOP = zeros(1,n_epoch);
 GDOP = zeros(1,n_epoch);
 
 for ti = 1: n_epoch
     
     dx = 100;
     % Step 1. 사용자 위치 원점 가정
     R_user = zeros(3,1).'; % 초기 위치 가정
     B_del = 0;              % B_u-B_r
     x_old = [R_user.';B_del];
     
     while(dx> 10^-4)
         
         H_dgps = zeros(n_vis,4);
         z_dgps = zeros(n_vis,1);

         for jj = 1: n_vis
             
             % Step 2. e_hat 계산
             R_su = user.svpos(ti,3*SV_vis(jj)-2:3*SV_vis(jj)) - R_user;
             e_hat = R_su/norm(R_su);
             
             % Step 3. 사용자 위치 계산 하기
             % Step 3.1. H matrix 계산
             H_dgps(jj,1:3) = e_hat;
             H_dgps(jj,4) = -1;
             % Step 3.2 d_r_j 계산 (correction)
             d_r_j=norm(ref.svpos(ti,3*SV_vis(jj)-2:3*SV_vis(jj))-ref_xyz);
                                                %기준국에서 위성까지 true거리
             % Step 3.2. measurement vector 계산
             z_dgps(jj,1) = e_hat*user.svpos(ti,3*SV_vis(jj) ...
                            -2:3*SV_vis(jj)).'-(user.pr1(ti,SV_vis(jj)) ...
                            -(ref.pr1(ti,SV_vis(jj))-d_r_j));
        
         end
         
         % Step 3.3. user position 계산
         x = pinv(H_dgps)*z_dgps;
         R_user = x(1:3).';
         B_del = x(4);
         
         % Step 4. 사용자 위치 비교 후 update
         dx = norm(x-x_old);
         x_old = x;
         
     end
     
     % DOP 계산
     H_enu = zeros(n_vis,4);
     for i = 1:n_vis
         H_enu(i,1:3) = (Rtran*(H_dgps(i,1:3)).').';
         H_enu(i,4) = -1;
     end
             
     Sat_geo = inv((H_enu).' * (H_enu));
     HDOP(1,ti) = sqrt((Sat_geo(1,1))+(Sat_geo(2,2)));
     VDOP(1,ti) = sqrt(Sat_geo(3,3));
     PDOP(1,ti) = sqrt((Sat_geo(1,1))+(Sat_geo(2,2))+(Sat_geo(3,3)));
     TDOP(1,ti) = sqrt(Sat_geo(4,4));
     GDOP(1,ti) = sqrt((Sat_geo(1,1))+(Sat_geo(2,2))+(Sat_geo(3,3)) ...
                  +(Sat_geo(4,4)));
    
     % 사용자 위치
     user.pos_dgps(ti,:) = R_user;
     user.pos_dgps_enu(ti,:) = (Rtran*(user.pos_dgps(ti,:)-ref_xyz).').';
                                            %기준국에서의 사용자의 enu좌표
     error.dgps(ti,:) = user.pos_dgps(ti,:)  - true.xyz(ti,:);
     
 end
 
 %% DGPS Plot
figure(1);
hold on;
plot(true.enu(ini_epoch:n_epoch,1),true.enu(ini_epoch:n_epoch,2),'r.');
plot(user.pos_dgps_enu(ini_epoch:n_epoch,1),user.pos_dgps_enu( ...
     ini_epoch:n_epoch,2),'k.','Marker','o','MarkerSize',4);
legend('True','DGPS');
title('DGPS Horizontal Trajectory');
xlabel('East(m)');
ylabel('North(m)');
hold off;
figure(2);
hold on;
plot(true.GPSTime(ini_epoch:n_epoch,1),true.enu(ini_epoch:n_epoch,3),'r.');
plot(user.GPSTime(ini_epoch:n_epoch,1),user.pos_dgps_enu( ...
     ini_epoch:n_epoch,3),'k.','Marker','o','MarkerSize',4);
legend('True','DGPS');
title('DGPS Vertical Trajectory');
xlabel('time(s)');
ylabel('Up(m)');
hold off;

%% RMS
add = 0;
add_temp = 0;
for ti=ini_epoch:n_epoch
    add_temp = (norm(error.dgps(ti,:)))^2;
    add = add+add_temp;
end
mean = add/(n_epoch-ini_epoch+1);
rms = sqrt(mean);

%% DOP 계산 (실제 움직인 시점 부터)
H = 0;
for ti=ini_epoch:n_epoch
    H_temp = HDOP(1,ti);
    H = H + H_temp;
end
H_DOP = H/(n_epoch-ini_epoch+1);

V = 0;
for ti=ini_epoch:n_epoch
    V_temp = VDOP(1,ti);
    V = V + V_temp;
end
V_DOP = V/(n_epoch-ini_epoch+1);

P = 0;
for ti=ini_epoch:n_epoch
    P_temp = PDOP(1,ti);
    P = P + P_temp;
end
P_DOP = P/(n_epoch-ini_epoch+1);

T = 0;
for ti=ini_epoch:n_epoch
    T_temp = TDOP(1,ti);
    T = T + T_temp;
end
T_DOP = T/(n_epoch-ini_epoch+1);

G = 0;
for ti=ini_epoch:n_epoch
    G_temp = GDOP(1,ti);
    G = G + G_temp;
end
G_DOP = G/(n_epoch-ini_epoch+1);