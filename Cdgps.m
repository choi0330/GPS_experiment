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
HEl_sat = 0;
HEl = 0;
for ii = 1:32
    
    % mask angle보다 elevation angle이 높고, 측정치가 있는 위성 찾기
    temp_el1 = find(ref.El(:,ii) < maskangle| user.El(:,ii) < maskangle); 
    temp_el2 = find(ref.pr1(:,ii)<10 |user.pr1(:,ii)<10);
    if isempty(temp_el1)
        if isempty(temp_el2)
            vis_sat(ii)=1;
        end
    end
    % 현재 관측가능한 위성과 그 위성의 Elevation평균 저장
    temp_sat1 = find(vis_sat==1);
    temp_hel1 = mean(user.El(:,temp_sat1));
end
SV_vis = find(vis_sat==1);
n_vis = length(SV_vis);
% Elevation평균이 제일 큰 위성 찾기
HEl = max(temp_hel1);
for jj = 1:length(SV_vis)
    if HEl == temp_hel1(jj)
        HEl_sat = temp_sat1(jj);
    end
end
% 기준 위성과 다른 위성 구분
SV_vis2 = zeros(1,length(SV_vis)-1);
j = 1;
for ii = 1:n_vis
    if SV_vis(1,ii) ~= HEl_sat
        SV_vis2(j) = SV_vis(1,ii);
        j = j+1;        
    end
end

    

%% CDGPS 
%  user.pos_cdgps = zeros(n_epoch,3);
 user.pos_cdgps_error = zeros(n_epoch,3);
 user.pos_cdgps_enu = zeros(n_epoch,3);
 R_ur = zeros(n_epoch,3);
 HDOP = zeros(1,n_epoch);
 VDOP = zeros(1,n_epoch);
 PDOP = zeros(1,n_epoch);
 TDOP = zeros(1,n_epoch);
 GDOP = zeros(1,n_epoch);
 c = 299792458;    % m/s
 f = 1575420000;   % Hz(1/s)
 Wave_length = c/f;% m
 cp_rju = zeros(n_epoch,n_vis-1);
 cp_rku = zeros(n_epoch,1);
 d_rju = zeros(n_epoch,n_vis-1);
 d_rku = zeros(n_epoch,1);
 cp_rjku = zeros(n_epoch,n_vis-1);
 d_rjku = zeros(n_epoch,n_vis-1);
 N_rjku = zeros(n_epoch,n_vis-1);
 N = zeros(1,n_vis-1);
 
 % Step 1. 미지정수 계산
for ti = 1:n_epoch
    cp_rku(ti,1) = ref.cp1(ti,HEl_sat)-user.cp1(ti,HEl_sat);
    d_rku(ti,1) = norm(ref.svpos(ti,3*HEl_sat-2:3*HEl_sat)-ref_xyz) ...
                 -norm(user.svpos(ti,3*HEl_sat-2:3*HEl_sat)-user_xyz0);
    for jj = 1:n_vis-1
     cp_rju(ti,jj) = ref.cp1(ti,SV_vis2(jj))-user.cp1(ti,SV_vis2(jj));
     d_rju(ti,jj) = norm(ref.svpos(ti,3*SV_vis2(jj)-2:3*SV_vis2(jj)) ...
     -ref_xyz)-norm(user.svpos(ti,3*SV_vis2(jj)-2:3*SV_vis2(jj))-user_xyz0);
     cp_rjku(ti,jj) = cp_rju(ti,jj)-cp_rku(ti,1);
     d_rjku(ti,jj) = d_rju(ti,jj)-d_rku(ti,1);
     N_rjku(ti,jj) = round((cp_rjku(ti,jj)-d_rjku(ti,jj))/(Wave_length));
     end 
end
% N 계산
for ii=1:n_vis-1
    N(1,ii) = round(mean(N_rjku(1:ini_epoch,ii)));
end
     
     
 for ti = 1: n_epoch
     
     dx = 100;
     % Step 2. 사용자 위치 원점 가정
     R_user = zeros(3,1).'; % 초기 위치 가정
     R_j = zeros(3,1).';
     R_k = zeros(3,1).';
     x_old = R_user.';
     
     while(dx> 10^-4)
         
         H_cdgps = zeros(n_vis-1,3);
         z_cdgps =  zeros(n_vis-1,1);
         
         R_k = user.svpos(ti,3*HEl_sat-2:3*HEl_sat);
         R_uk = R_k - R_user;
         e_hatk = R_uk/norm(R_uk);
         R_rk = ref.svpos(ti,3*HEl_sat-2:3*HEl_sat) - ref_xyz;
         e_hatrk = R_rk/norm(R_rk);

         
         for jj = 1: n_vis-1
             
             % Step 3.1 e_hat 계산
             R_j = user.svpos(ti,3*SV_vis2(jj)-2:3*SV_vis2(jj));
             R_su = R_j - R_user;
             e_hat = R_su/norm(R_su);
             R_sr = ref.svpos(ti,3*SV_vis2(jj)-2:3*SV_vis2(jj)) - ref_xyz;
             e_hatr = R_sr/norm(R_sr);
             

             % Step 3.2. H matrix 계산
             H_cdgps(jj,1:3) = e_hat-e_hatk;
             % Step 3.3. z matrix 계산
             z_cdgps(jj,1) = cp_rjku(ti,jj)-N(jj)*Wave_length-((e_hatr ...
             -e_hat)*R_j.'-(e_hatrk-e_hatk)*R_k.'-(e_hatr-e_hatrk)*ref_xyz.');
         end
         
         % Step 3.5. user position 계산
         x = pinv(H_cdgps)*z_cdgps;
         R_user = x(1:3).';
         
         % Step 4. 사용자 위치 비교 후 update
         dx = norm(x-x_old);
         x_old = x;
         
     end
          
     % DOP 계산
     H_enu = zeros(n_vis-1,3);
     for i = 1:n_vis-1
         H_enu(i,1:3) = (Rtran*(H_cdgps(i,1:3)).').';
     end
             
     Sat_geo = inv((H_enu).' * (H_enu));
     HDOP(1,ti) = sqrt((Sat_geo(1,1))+(Sat_geo(2,2)));
     VDOP(1,ti) = sqrt(Sat_geo(3,3));
     PDOP(1,ti) = sqrt((Sat_geo(1,1))+(Sat_geo(2,2))+(Sat_geo(3,3)));
    
     % 사용자 위치
     user.pos_cdgps(ti,:) = R_user;
     user.pos_cdgps_enu(ti,:) = (Rtran*(user.pos_cdgps(ti,:)-ref_xyz).').';
     error.cdgps(ti,:) = user.pos_cdgps(ti,:)  - true.xyz(ti,:);
     
 end
 
 %% CDGPS Plot
figure(1);
hold on;
plot(true.enu(ini_epoch:n_epoch,1),true.enu(ini_epoch:n_epoch,2),'r.');
plot(user.pos_cdgps_enu(ini_epoch:n_epoch,1),user.pos_cdgps_enu( ...
     ini_epoch:n_epoch,2),'k.','Marker','o','MarkerSize',4);
legend('True','CDGPS');
title('CDGPS Horizontal Trajectory');
xlabel('East(m)');
ylabel('North(m)');
hold off;
figure(2);
hold on;
plot(true.GPSTime(ini_epoch:n_epoch,1),true.enu(ini_epoch:n_epoch,3),'r.');
plot(user.GPSTime(ini_epoch:n_epoch,1),user.pos_cdgps_enu( ...
     ini_epoch:n_epoch,3),'k.','Marker','o','MarkerSize',4);
legend('True','CDGPS');
title('CDGPS Vertical Trajectory');
xlabel('time(s)');
ylabel('Up(m)');
hold off;

%% RMS
add = 0;
add_temp = 0;
for ti=ini_epoch:n_epoch
    add_temp = (norm(error.cdgps(ti,:)))^2;
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