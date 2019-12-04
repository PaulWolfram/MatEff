% This script contains all calculations to produce the results shown in 
% "Material efficiency for immediate climate change mitigation of passenger vehicles",
% Journal of Industrial Ecology

% Paul Wolfram, 9/11/2019 
% Last updated 12/4/2019

tic
clear all
close all
clc
cd 'C:/Users/pw379/Box Sync/_New_projects/G7 RECC priv/___JIE SI'
% specify folder that contains data file and script


%% Read in data 
file = 'ME_data_Nov13.xlsx';
% specify data file (needs to be in the same folder as the script)

% Material composition of baseline and LW vehicles:
[x_mat,hdrs_x_mat,all_x_mat] = xlsread(file,'Mat_com','A3:I195');
x_mat_rel = xlsread(file,'Mat_com_%','B4:I51');

% Vehicle weight by component
x_comp = xlsread(file,'Comp','B4:J51');

% Material emission factors:
f_mat_vir = xlsread(file,'Mat_EF','D4:E11'); % virgin
f_mat_rec = xlsread(file,'Mat_EF_rec','D4:E11'); % recycled

% Energy carrier emission factors:
f_en = xlsread(file,'En_EF','C4:D9');

% Assembly energy use factors (kWh/kg vehicle):
x_en_ass = xlsread(file,'Ass_en','B4:G27');
x_en_ass_reu = xlsread(file,'Reu_en','B4:G27');

% Drive-cycle energy use factors by fuel:
[x_en_dc,hdrs_x_en,all_x_en] = xlsread(file,'Use_en','A3:G51');

% Drive-cycle energy use factors total:
x_en_tot = xlsread(file,'Use_en','H4:H51');

% Material recycling rates:
rec_mat = xlsread(file,'Rec','B4:I27');

% Component Remanufacturing rates:
reu_mat = xlsread(file,'Reu','B4:I27');
reu_mat = reu_mat ./ 1.02;
% (corrected by the fact that about 2% of new materials are needed during Remanufacturing)

% Define other parameters:
para_occ = 1.5; % Vehicle occupancies
para_occ_MIU = 2.0;
para_ltm = 150000; % Lifetime mileage (km)


%% Add material composition data of other MESs for virign material production stage:
x_mat(isnan(x_mat)) = 0;
x_mat_prd = [x_mat; zeros(size(x_mat,1)/2*6,size(x_mat,2))]; % initialize as zeros

% Rec:
x_mat_prd(2*24+1:3*24,:) = x_mat(1:24,:) .* (1 - rec_mat);
% pick baseline and apply recycling rates

% Rem:
x_mat_prd(3*24+1:4*24,:) = x_mat(1:24,:) .* (1 - reu_mat) ;
% pick baseline and apply Remanufacturing rates 

% DWN:
% switch vehicle segment rows sothat ... 
for d = 0:5;
    x_mat_prd(4*(24+d)+1,:) = x_mat(4*d+1,:); % micro remains micro because DWN is not possible for micro
    x_mat_prd(4*(24+d)+2,:) = x_mat(4*d+1,:); % ... PC becomes micro
    x_mat_prd(4*(24+d)+3,:) = x_mat(4*d+2,:); % ... SUV becomes PC
    x_mat_prd(4*(24+d)+4,:) = x_mat(4*d+3,:); % ... LT becomes SUV
end

% MIU:
x_mat_prd(5*24+1:6*24,:) = x_mat(1:24,:) .* (para_occ / para_occ_MIU);
% pick baseline vehicle and adjust occupancy factor
% (we haven't YET converted to t CO2e/p!)

% All:
% switch vehicle segment rows sothat ...
for d = 0:5;
    x_mat_prd(6*24+1+(d*4),:) = x_mat(4*d+1+24,:) .* (1 - rec_mat(1,:)) .* (1 - reu_mat(d*4+1,:)) .* (para_occ / para_occ_MIU); 
    % ... micro becomes LW micro, then apply REC, ReU and MIU
    % we can always pick first recycling row because same for all powertrains and segments 
    % but we pick powertrain-specific reman. rate (same for all segments)
    x_mat_prd(6*24+2+(d*4),:) = x_mat(4*d+1+24,:) .* (1 - rec_mat(1,:)) .* (1 - reu_mat(d*4+1,:)) .* (para_occ / para_occ_MIU); 
    % ... PC becomes LW micro, then apply REC, ReU and MIU
    x_mat_prd(6*24+3+(d*4),:) = x_mat(4*d+2+24,:) .* (1 - rec_mat(1,:)) .* (1 - reu_mat(d*4+1,:)) .* (para_occ / para_occ_MIU); 
    % ... SUV becomes LW PC, then apply REC, ReU and MIU
    x_mat_prd(6*24+4+(d*4),:) = x_mat(4*d+3+24,:) .* (1 - rec_mat(1,:)) .* (1 - reu_mat(d*4+1,:)) .* (para_occ / para_occ_MIU); 
    % ... LT becomes LW SUV, then apply REC, ReU and MIU
end

% All but LW:
x_mat_prd(7*24+1:8*24,:) = x_mat_prd(4*24+1:5*24,:) .* (1 - rec_mat) .* (1 - reu_mat) .* (para_occ / para_occ_MIU);
% pick DWN and apply REC, ReU and MIU factors 


%% Add material composition data for material recycling stage:
x_mat_rec = zeros(192,8);

% Rec:
x_mat_rec(2*24+1:3*24,:) = x_mat(1:24,:) .* rec_mat;
% pick baseline and apply recycling rates

% All:
% switch vehicle segment rows sothat ...
for d = 0:5;
    x_mat_rec(6*24+1+(d*4),:) = x_mat(4*d+1+24,:) .* rec_mat(1,:) .* (1 - reu_mat(d*4+1,:)) .* (para_occ / para_occ_MIU); 
    % ... micro becomes LW micro, then apply REC, ReU and MIU
    x_mat_rec(6*24+2+(d*4),:) = x_mat(4*d+1+24,:) .* rec_mat(1,:) .* (1 - reu_mat(d*4+1,:)) .* (para_occ / para_occ_MIU); 
    % ... PC becomes LW micro, then apply REC, ReU and MIU
    x_mat_rec(6*24+3+(d*4),:) = x_mat(4*d+2+24,:) .* rec_mat(1,:) .* (1 - reu_mat(d*4+1,:)) .* (para_occ / para_occ_MIU); 
    % ... SUV becomes LW PC, then apply REC, ReU and MIU
    x_mat_rec(6*24+4+(d*4),:) = x_mat(4*d+3+24,:) .* rec_mat(1,:) .* (1 - reu_mat(d*4+1,:)) .* (para_occ / para_occ_MIU); 
    % ... LT becomes LW SUV, then apply REC, ReU and MIU
end

% All but LW:
x_mat_rec(7*24+1:8*24,:) = x_mat(1:24,:) .* rec_mat .* (1 - reu_mat) .* (para_occ / para_occ_MIU);
% pick DWN and apply REC, ReU and MIU factors 


%% Caculate total mass that enters vehicle assembly stage:
% same for most strategies but for REC, 'all' and 'all but LW' we need 
% to add mass of recycled materials because these in fact do enter assembly 
% stage (but did not enter virgin material stage)
 
x_veh_ass = x_mat_prd;

%REC:
x_veh_ass(2*24+1:3*24,:) = x_mat(1:24,:);
% pick baseline, do not apply REC factor

% All:
% switch vehicle segment rows sothat ...
for d = 0:5;
    x_veh_ass(6*24+1+(d*4),:) = x_mat(4*d+1+24,:) .* (1 - reu_mat(d*4+1,:)) .* (para_occ / para_occ_MIU); 
    % ... micro becomes LW micro, then apply ReU and MIU, not REC
    x_veh_ass(6*24+2+(d*4),:) = x_mat(4*d+1+24,:) .* (1 - reu_mat(d*4+1,:)) .* (para_occ / para_occ_MIU); 
    % ... PC becomes LW micro, then apply ReU and MIU, not REC
    x_veh_ass(6*24+3+(d*4),:) = x_mat(4*d+2+24,:) .* (1 - reu_mat(d*4+1,:)) .* (para_occ / para_occ_MIU); 
    % ... SUV becomes LW PCV, then apply ReU and MIU, not REC
    x_veh_ass(6*24+4+(d*4),:) = x_mat(4*d+3+24,:) .* (1 - reu_mat(d*4+1,:)) .* (para_occ / para_occ_MIU); 
    % ... LT becomes LW SUV, then apply ReU and MIU, not REC
end

% All but LW:
x_veh_ass(7*24+1:8*24,:) = x_mat_prd(4*24+1:5*24,:) .* (1 - reu_mat) .* (para_occ / para_occ_MIU);
% pick DWN and apply ReU and MIU factors but not REC

% Calculate totals:
x_veh_ass_tot = sum(x_veh_ass,2);


%% Add remanufacturing layer
% initilize as zeros
x_veh_ass_reu = zeros(192,8);

% Rem:
x_veh_ass_reu(73:96,:) = x_mat(1:24,:) .* reu_mat;
% pick baseline and apply reman. rates

% All:
% switch vehicle segment rows sothat ...
for d = 0:5;
    x_veh_ass_reu(6*24+1+(d*4),:) = x_mat(4*d+1+24,:) .* reu_mat(d*4+1,:) .* (para_occ / para_occ_MIU); 
    % ... micro becomes LW micro, then apply ReU and MIU
    x_veh_ass_reu(6*24+2+(d*4),:) = x_mat(4*d+1+24,:) .* reu_mat(d*4+1,:) .* (para_occ / para_occ_MIU); 
    % ... PC becomes LW micro, then apply ReU and MIU
    x_veh_ass_reu(6*24+3+(d*4),:) = x_mat(4*d+2+24,:) .* reu_mat(d*4+1,:) .* (para_occ / para_occ_MIU); 
    % ... SUV becomes LW PC, then apply ReU and MIU
    x_veh_ass_reu(6*24+4+(d*4),:) = x_mat(4*d+3+24,:) .* reu_mat(d*4+1,:) .* (para_occ / para_occ_MIU); 
    % ... LT becomes LW SUV, then apply ReU and MIU
end

% All but LW:
x_veh_ass_reu(7*24+1:8*24,:) = x_mat(1:24,:) .* reu_mat .* (para_occ / para_occ_MIU);
% pick DWN and apply ReU and MIU factors 

% Total:
x_veh_ass_reu_tot = sum(x_veh_ass_reu,2);


%% Plot fuel consumption, weight and material composition (Figure 1)
font = 7;
labels = hdrs_x_mat(2:193,1);
x = 1:24;

figure
% Subplot a: energy consumption 
subplot(2,2,1)
scatter(x,x_en_tot(1:24),'MarkerFaceColor',[.5 .5 .5],'MarkerEdgeColor',[0 0 0]) % fuel consumption of baseline vehicle
hold on
scatter(x,x_en_tot(25:48),'MarkerFaceColor',[1 .8 .4],'MarkerEdgeColor',[0 0 0]) % LW vehicle
hold off
title('a) Vehicle energy consumption')
ylabel('kWh (100 km)^{-1}')
xticks(1:1:24)
set(gca, 'XTickLabel', {[]},'FontSize',font,'YGrid','on','XGrid','on')
legend({'Baseline', 'Lightweight'},'Location','south','Orientation','vertical')
alpha(.75)


% Subplot b: weight by component
subplot(2,2,3)
bar0 = bar(x_comp(1:24,1:8),'stacked'); % baseline vehicle weight
hold on 
scatter(x,x_comp(25:48,9),'MarkerEdgeColor',[.25 .25 .25],'MarkerFaceColor',[.25 .25 .25]) % LW vehicle weight
hold off
title('b) Weight by component')
ylabel('kg')
xticks(1:1:24)
set(gca,'xticklabel', labels(1:1:24),'FontSize',font,'YGrid','on','XGrid','on','XTickLabelRotation',45)
legend2 = legend({'Engine', 'Motor','Fuel cell','Transmission','Battery','Tank','Chassis','Body','Total weight, lightweight'},'Location','best','Orientation','vertical');
legend2.NumColumns = 3;
alpha(.75)

bar0(1).FaceColor = [.145 .702 .890]; % light blue
bar0(2).FaceColor = [.098 .584 .757]; % medium blue
bar0(3).FaceColor = [.071 .435 .565]; % dark blue
bar0(4).FaceColor = [.337 .824 .533]; % light green
bar0(5).FaceColor = [.988 .439 .439]; % light red
bar0(6).FaceColor = [1 .078 .078]; % red
bar0(7).FaceColor = [.5 .5 .5]; % grey
bar0(8).FaceColor = [.25 .25 .25]; % darker grey 

for i =1:8
    bar0(i).EdgeAlpha = 0;
end


% Subplot c: material composition, baseline
subplot(2,2,2)
barx = bar(x_mat_rel(1:24,:),'stacked');
title('c) Material composition, baseline')
ylim([0 1])
xticks(1:1:24)
set(gca,'xticklabel',{[]},'FontSize',font,'YGrid','on','XGrid','on')
legend3 = legend({'Autom. steel', 'Stainl. steel','Cast iron','Wrought Al','Cast Al','Copper','Plastics','Other'},'Location','south','Orientation','horizontal');
legend3.NumColumns = 3;
alpha(.75)

for i =1:8
     barx(i).EdgeAlpha = 0;
end

 barx(1).FaceColor = [.098 .584 .757]; % light blue
 barx(2).FaceColor = [0.145 0.702 0.89]; % medium blue
 barx(3).FaceColor = [0.071 0.435 0.565]; % medium blue
 barx(4).FaceColor = [0.4 0.4 0.4]; % dark grey 
 barx(5).FaceColor = [0.6 0.6 0.6]; % grey 
 barx(6).FaceColor = [0.84 .37 0.3]; % orange
 barx(7).FaceColor = [1 .904 .402]; % light yellow
 barx(8).FaceColor = [.949 .745 0]; % medium yellow


% Subplot d: material composition, lightweight:
subplot(2,2,4)
bary = bar(x_mat_rel(25:48,:),'stacked');
title('d) Material composition, lightweight')
ylim([0 1])
xticks(1:1:24)
set(gca, 'XTickLabel', labels(1:1:24),'XTickLabelRotation',45,'FontSize',font,'YGrid','on','XGrid','on')
alpha(.75)

 for i =1:8
     bary(i).EdgeAlpha = 0;
 end

 bary(1).FaceColor = [.098 .584 .757]; %light blue
 bary(2).FaceColor = [0.145 0.702 0.89]; % medium blue
 bary(3).FaceColor = [0.071 0.435 0.565]; % dark blue
 bary(4).FaceColor = [0.4 0.4 0.4]; % dark grey 
 bary(5).FaceColor = [0.6 0.6 0.6]; % grey 
 bary(6).FaceColor = [0.84 .37 0.3]; % orange
 bary(7).FaceColor = [1 .904 .402]; % light yellow
 bary(8).FaceColor = [.949 .745 0]; % medium yellow


%% Calculate material production and recycling footprint:
F_mat_vir = (x_mat_prd * f_mat_vir) ./ para_occ; % 1st col. = conv., 2nd col. = low-carbon
F_mat_rec = (x_mat_rec * f_mat_rec) ./ para_occ;
F_mat = F_mat_vir + F_mat_rec;


%% Calculate vehicle assembly and reman. footprint:
f_ass = x_en_ass * f_en / 1000; % assembly emission factors (kg CO2e/kg vehicle)
f_ass = repmat(f_ass,8,1);

f_ass_reu = x_en_ass_reu * f_en / 1000; % remanufacturing emission factors
f_ass_reu = repmat(f_ass_reu,8,1);
% 1st col. = conv., 2nd col. = low-carbon/CCS

F_ass = (f_ass .* [x_veh_ass_tot x_veh_ass_tot]) ./ para_occ;

F_ass_reu = (f_ass_reu .* [x_veh_ass_reu_tot x_veh_ass_reu_tot]) ./ para_occ

F_ass = F_ass + F_ass_reu;


%% Fuel cycle (well-to-wheel):
% wtw emission factors (g CO2e/kWh * kWh/100km / 100 = g CO2e/km)
f_wtw = zeros(192,2); % initialize as zeros

% baseline and LW:
f_wtw(0*24+1:2*24,:) = x_en_dc * f_en / 100;  

% REC = baseline:
f_wtw(2*24+1:3*24,:) = f_wtw(1:24,:); 

% ReU = baseline:
f_wtw(3*24+1:4*24,:) = f_wtw(1:24,:); 

% DWN:
% switch vehicle segment rows sothat ... 
for d = 0:5;
    f_wtw(4*(24+d)+1,:) = f_wtw(4*d+1,:); % micro remains micro because DWN is not possible for micro
    f_wtw(4*(24+d)+2,:) = f_wtw(4*d+1,:); % ... PC becomes micro
    f_wtw(4*(24+d)+3,:) = f_wtw(4*d+2,:); % ... SUV becomes PC
    f_wtw(4*(24+d)+4,:) = f_wtw(4*d+3,:); % ... LT becomes SUV
end

% MIU:
f_wtw(5*24+1:6*24,:) = f_wtw(1:24,:) .* (para_occ / para_occ_MIU);
% pick baseline and scale with occupancy factor

% All but LW:
f_wtw(7*24+1:8*24,:) = f_wtw(4*24+1:5*24,:) .* (para_occ / para_occ_MIU);
% pick DWN and apply MIU factors 

% All:
% switch vehicle segment rows sothat ...
for d = 0:5;
    f_wtw(6*24+1+(d*4),:) = f_wtw(4*d+1+24,:) .* (para_occ / para_occ_MIU); 
    % ... micro becomes LW micro, then apply MIU
    f_wtw(6*24+2+(d*4),:) = f_wtw(4*d+1+24,:) .* (para_occ / para_occ_MIU); 
    % ... PC becomes LW micro, then apply MIU
    f_wtw(6*24+3+(d*4),:) = f_wtw(4*d+2+24,:) .* (para_occ / para_occ_MIU); 
    % ... SUV becomes LW PC, then apply MIU
    f_wtw(6*24+4+(d*4),:) = f_wtw(4*d+3+24,:) .* (para_occ / para_occ_MIU); 
    % ... LT becomes LW SUV, then apply MIU
end

F_wtw = (f_wtw .* para_ltm / 1000) ./ para_occ;


%% Total:
F_tot = F_mat + F_ass + F_wtw;


%% Print highest and lowest totals to command window:
% Current:
[min_cur,idx_min_cur] = min(F_tot(:,1)); % find smallest value and its index
all_x_mat(idx_min_cur+1,1); % find corresponding label
fprintf('\nCurrent:\n')
fprintf('- With %.2f t CO2e, %s has the smallest life-cycle GHG footprint per person.\n',min_cur/1000,all_x_mat{idx_min_cur+1,1})
fprintf('   - Materials contribute %.0f%%. ',F_mat(idx_min_cur,1)/min_cur*100)
fprintf('Vehicle assembly contributes %.0f%%. ',F_ass(idx_min_cur,1)/min_cur*100)
fprintf('The fuel cycle contributes %.0f%%.',F_wtw(idx_min_cur,1)/min_cur*100)

[max_cur,idx_max_cur] = max(F_tot(:,1)); % find largest value and its index
all_x_mat(idx_max_cur+1,1); % find corresponding label
fprintf('\n- With %.2f t CO2e, %s has the largest life-cycle GHG footprint per person.\n',max_cur/1000,all_x_mat{idx_max_cur+1,1})
fprintf('   - Materials contribute %.0f%%. ',F_mat(idx_max_cur,1)/max_cur*100)
fprintf('Vehicle assembly contributes %.0f%%. ',F_ass(idx_max_cur,1)/max_cur*100)
fprintf('The fuel cycle contributes %.0f%%.\n',F_wtw(idx_max_cur,1)/max_cur*100)

% Low-carbon:
[min_lo,idx_min_lo] = min(F_tot(:,2)); % find smallest value and its index
all_x_mat(idx_min_lo+1,1); % find corresponding label
fprintf('\nLow-carbon:\n');
fprintf('- With %.2f t CO2e, %s has the smallest life-cycle GHG footprint per person.\n',min_lo/1000,all_x_mat{idx_min_lo+1,1})
fprintf('   - Materials contribute %.0f%%. ',F_mat(idx_min_lo,2)/min_lo*100)
fprintf('Vehicle assembly contributes %.0f%%. ',F_ass(idx_min_lo,2)/min_lo*100)
fprintf('The fuel cycle contributes %.0f%%.',F_wtw(idx_min_lo,2)/min_lo*100)

[max_lo,idx_max_lo] = max(F_tot(:,2)); % find largest value and its index
all_x_mat(idx_max_lo+1,1); % find corresponding label
fprintf('\n- With %.2f t CO2e, %s has the largest life-cycle GHG footprint per person.\n',max_lo/1000,all_x_mat{idx_max_lo+1,1})
fprintf('   - Materials contribute %.0f%%. ',F_mat(idx_max_lo,2)/max_lo*100)
fprintf('Vehicle assembly contributes %.0f%%. ',F_ass(idx_max_lo,2)/max_lo*100)
fprintf('The fuel cycle contributes %.0f%%.\n',F_wtw(idx_max_lo,1)/max_lo*100)


%% Plot all, disaggregated by life cycle phase (Figure 2):
wtw_col = [.00 .62 .45]; % green
ass_col = [1 .078 .078]; % red
mat_col = [.988 .439 .439]; % light red 
skip = 4;

plot_data=[(F_tot(1:192,1)/1000)]; 
[plot_data_sort,idx_sort]=sort(plot_data,'ascend'); % sort totals and extract index

figure
% current:
subplot(1,2,1)
plot_stack_data=[F_mat(idx_sort,1)/1000, F_ass(idx_sort,1)/1000, F_wtw(idx_sort,1)/1000]; % use same index to sort life cycle phases
bar1 = barh(plot_stack_data,'stacked');
bar1(1).FaceColor = mat_col;
bar1(2).FaceColor = ass_col;
bar1(3).FaceColor = wtw_col;
bar1(1).EdgeAlpha = 0;
bar1(2).EdgeAlpha = 0;
bar1(3).EdgeAlpha = 0;

labels = labels(idx_sort);
labels_skip = labels(1:skip:192);

set(gca, 'YTickLabel', labels_skip,'YTickLabelRotation',0,'FontSize',font,'YGrid','on','XGrid','on')
yticks(1:skip:192)
xlabel('t CO_{2}e person^{-1}')
title('a) Current')
alpha(.7)
xlim([0 45])

% low-carbon:
plot_data_lo=[(F_tot(:,2)/1000)]; 
[plot_data_lo_sort,idx_sort_lo]=sort(plot_data_lo,'ascend');

subplot(1,2,2)
plot_stack_data_lo=[F_mat(idx_sort_lo,2)/1000, F_ass(idx_sort_lo,2)/1000, F_wtw(idx_sort_lo,2)/1000];
bar2 = barh(plot_stack_data_lo,'stacked');
bar2(1).FaceColor = mat_col;
bar2(2).FaceColor = ass_col;
bar2(3).FaceColor = wtw_col;
bar2(1).EdgeAlpha = 0;
bar2(2).EdgeAlpha = 0;
bar2(3).EdgeAlpha = 0;

labels_lo = [all_x_mat(2:193,1)];
labels_lo = labels_lo(idx_sort_lo);
labels_lo_skip = labels_lo(1:4:192);
set(gca, 'YTickLabel', labels_lo_skip,'YTickLabelRotation',0,'FontSize',font,'YGrid','on','XGrid','on')
yticks(1:skip:192)
xlabel('t CO_{2}e person^{-1}')
title('b) Low-carbon') % Carbon footprint per person over vehicle lifetime
legend({'Material production', 'Vehicle assembly', 'Energy cycle'},'Location','bestoutside','Orientation','horizontal')
alpha(.7)
xlim([0 45])


%% Plot lowest and highest only (figure not used in manuscript)
figure
% current:
subplot(2,1,1)
bar([plot_data_sort(1:30); 0; 0; 0; plot_data_sort(163:192)])
labels=[all_x_mat(2:193,1)];
labels=labels(idx_sort);
labels_maxmin=[labels;hdrs_x_mat(1,1);hdrs_x_en(1,1)]; % add an empty label (first cell of hdrs_x_mat) to labels
labels_maxmin=labels_maxmin([1:30, 194, 193, 194, 163:192]);
set(gca, 'XTickLabel', labels_maxmin,'XTickLabelRotation',45,'FontSize',9,'YGrid','on','XGrid','on')
xticks(1:1:63)
ylabel('t CO_{2}e person^{-1}')
title('a) Carbon footprint per person over vehicle lifetime, current')
alpha(.7)

% low-carbon:
subplot(2,1,2)
bar([plot_data_lo_sort(1:30); 0; 0; 0; plot_data_lo_sort(163:192)])
labels_lo=[all_x_mat(2:193,1)];
labels_lo=labels_lo(idx_sort_lo);
labels_lo_maxmin=[labels_lo;hdrs_x_mat(1,1);hdrs_x_en(1,1)]; % add an empty label (first cell of hdrs_x_mat) to labels
labels_lo_maxmin=labels_lo_maxmin([1:30, 194, 193, 194, 163:192]);
set(gca, 'XTickLabel', labels_lo_maxmin,'XTickLabelRotation',45,'FontSize',9,'YGrid','on','XGrid','on')
xticks(1:1:63)
ylabel('t CO_{2}e person^{-1}')
title('b) Carbon footprint per person over vehicle lifetime, low-carbon')
alpha(.7)


%% Plot highest and lowest, disaggregated by life cycle phase (figure not used in manuscript):

plot_stack_data = zeros(51,3); % initialize data matrix

plot_stack_data(1:30, 1) = F_mat(idx_sort(1:30), 1)/1000;
plot_stack_data(31, 1) = 0;
plot_stack_data(32:end, 1) = F_mat(idx_sort(192-20+1:end), 1)/1000;

plot_stack_data(1:30, 2) = F_ass(idx_sort(1:30), 1)/1000;
plot_stack_data(31, 2) = 0;
plot_stack_data(32:end, 2) = F_ass(idx_sort(192-20+1:end), 1)/1000;

plot_stack_data(1:30, 3) = F_wtw(idx_sort(1:30), 1)/1000;
plot_stack_data(31, 3) = 0;
plot_stack_data(32:end, 3) = F_wtw(idx_sort(192-20+1:end), 1)/1000;

figure
% current:
subplot(1,2,1)
bar3 = barh(plot_stack_data,'stacked');
bar3(1).FaceColor = mat_col;
bar3(2).FaceColor = ass_col;
bar3(3).FaceColor = wtw_col;
bar3(1).EdgeAlpha = 0;
bar3(2).EdgeAlpha = 0;
bar3(3).EdgeAlpha = 0;

labels=[all_x_mat(2:193,1)];
labels_cur=labels(idx_sort);

labels_maxmin_dis_cur=[labels_cur;hdrs_x_mat(1,1);hdrs_x_en(1,1)];
labels_maxmin_dis_cur=labels_maxmin_dis_cur([1:30, 193, 173:end]);
set(gca, 'YTickLabel', labels_maxmin_dis_cur,'YTickLabelRotation',0,'FontSize',7,'XGrid','on','YGrid','on')
yticks(1:1:51)
xlabel('t CO_{2}e person^{-1}')
title('a) Current')
alpha(.7)

% low-carbon:
plot_stack_data_lo = zeros(51,3); % initialize data matrix

plot_stack_data_lo(1:30, 1) = F_mat(idx_sort_lo(1:30), 2)/1000;
plot_stack_data_lo(31, 1) = 0;
plot_stack_data_lo(32:end, 1) = F_mat(idx_sort_lo(192-20+1:end), 2)/1000;

plot_stack_data_lo(1:30, 2) = F_ass(idx_sort_lo(1:30), 2)/1000;
plot_stack_data_lo(31, 2) = 0;
plot_stack_data_lo(32:end, 2) = F_ass(idx_sort_lo(192-20+1:end), 2)/1000;

plot_stack_data_lo(1:30, 3) = F_wtw(idx_sort_lo(1:30), 2)/1000;
plot_stack_data_lo(31, 3) = 0;
plot_stack_data_lo(32:end, 3) = F_wtw(idx_sort_lo(192-20+1:end), 2)/1000;

subplot(1,2,2)
bar3 = barh(plot_stack_data_lo,'stacked');
bar3(1).FaceColor = mat_col;
bar3(2).FaceColor = ass_col;
bar3(3).FaceColor = wtw_col;
bar3(1).EdgeAlpha = 0;
bar3(2).EdgeAlpha = 0;
bar3(3).EdgeAlpha = 0;

labels_lo=[all_x_mat(2:193,1)];
labels_lo=labels_lo(idx_sort_lo);

labels_maxmin_dis_lo=[labels_lo;hdrs_x_mat(1,1);hdrs_x_en(1,1)];
labels_maxmin_dis_lo=labels_maxmin_dis_lo([1:30, 193, 173:end]);
set(gca, 'YTickLabel', labels_maxmin_dis_lo,'YTickLabelRotation',0,'FontSize',7,'XGrid','on','YGrid','on')
legend({'Material production', 'Vehicle assembly', 'Energy cycle'},'Location','bestoutside','Orientation','horizontal')
yticks(1:1:51)
xlabel('t CO_{2}e person^{-1}')
title('b) Low-carbon')
alpha(.7)
%close


%% Plot share wtw/total (lowest and highest) (figure not used in manuscript)
print_range_lo=1:20;
print_range_hi=192-20+1:192; % range of lowest and highest that we wish to plot

% current:
figure
subplot(1,2,1)
F_wtw_shr = F_wtw ./ F_tot; % WTW emissions share
F_wtw_shr_cur = F_wtw_shr(:,1); % current vector
[F_wtw_shr_cur_sort,idx_wtw_shr_cur_sort]=sort(F_wtw_shr_cur,'ascend'); % sort from low to high
F_wtw_shr_cur_sort_plot = [F_wtw_shr_cur_sort(print_range_lo); 0; F_wtw_shr_cur_sort(print_range_hi)]; 
F_veh_prd_shr_cur_sort = 1-F_wtw_shr_cur_sort; % vehicle prod. emissions = 1 - WTW emissions
F_veh_prd_shr_cur_sort_plot = [F_veh_prd_shr_cur_sort(print_range_lo); 0; F_veh_prd_shr_cur_sort(print_range_hi)]; % range to be plotted
bar5 = barh([F_wtw_shr_cur_sort_plot F_veh_prd_shr_cur_sort_plot],'stacked'); % bar chart stacked
bar5(1).FaceColor = wtw_col;
bar5(2).FaceColor = ass_col;
bar5(1).EdgeAlpha = 0;
bar5(2).EdgeAlpha = 0;

labels=[all_x_mat(2:193,1)]; % original labels 
labels_wtw_cur_sort=labels(idx_wtw_shr_cur_sort); % sort labels 
labels_wtw_cur_sort_rng=[labels_wtw_cur_sort(print_range_lo); all_x_mat(1,1); labels_wtw_cur_sort(print_range_hi)]; % range of labels to print
set(gca, 'YTickLabel', labels_wtw_cur_sort_rng,'YTickLabelRotation',0,'FontSize',9,'XGrid','on','YGrid','on')
yticks(1:1:41)
title('a) Current')
alpha(.7)

% low-carbon:
subplot(1,2,2)
F_wtw_shr = F_wtw ./ F_tot; % WTW emissions share
F_wtw_shr_lo = F_wtw_shr(:,2); % low-carbon vector
[F_wtw_shr_lo_sort,idx_wtw_shr_lo_sort]=sort(F_wtw_shr_lo,'ascend'); % sort from low to high
F_wtw_shr_lo_sort_plot = [F_wtw_shr_lo_sort(print_range_lo); 0; F_wtw_shr_lo_sort(print_range_hi)]; % range of lowest and highest that we wish to plot
F_veh_prd_shr_lo_sort = 1-F_wtw_shr_lo_sort; % vehicle prod. emissions = 1 - WTW emissions
F_veh_prd_shr_lo_sort_plot = [F_veh_prd_shr_lo_sort(print_range_lo); 0; F_veh_prd_shr_lo_sort(print_range_hi)]; % range to be plotted
bar6 = barh([F_wtw_shr_lo_sort_plot F_veh_prd_shr_lo_sort_plot],'stacked'); % bar chart stacked
bar6(1).FaceColor = wtw_col;
bar6(2).FaceColor = ass_col;
bar6(1).EdgeAlpha = 0;
bar6(2).EdgeAlpha = 0;

labels=[all_x_mat(2:193,1)]; % original labels 
labels_wtw_lo_sort=labels(idx_wtw_shr_lo_sort); % sort labels 
labels_wtw_lo_sort_rng=[labels_wtw_lo_sort(print_range_lo); all_x_mat(1,1); labels_wtw_lo_sort(print_range_hi)]; % range of labels to print
set(gca, 'YTickLabel', labels_wtw_lo_sort_rng,'YTickLabelRotation',0,'FontSize',9,'XGrid','on','YGrid','on')
yticks(1:1:41)
title('b) Low-carbon')
legend({'Energy cycle','Vehicle supply chain'},'Location','bestoutside', 'Orientation','horizontal','Position',[0.5 0.95 0.3 0.02])
alpha(.7)
%close


%% Plot all WTW shares (Figure 3)
% calculate averages for each powertrain/strategy bin 
wtw_shr_avg = zeros(48,2);
% current:
for i = 1:48;
        wtw_shr_avg(i,1) = mean(F_wtw_shr_cur(1+(4*(i-1)):4*i));
end
% low carbon:
for i = 1:48;
        wtw_shr_avg(i,2) = mean(F_wtw_shr_lo(1+(4*(i-1)):4*i));
end

% caculate min
wtw_shr_min = zeros(48,2);
% current:
for i = 1:48;
        wtw_shr_min(i,1) = min(F_wtw_shr_cur(1+(4*(i-1)):4*i));
end
% low-carbon:
for i = 1:48;
        wtw_shr_min(i,2) = min(F_wtw_shr_lo(1+(4*(i-1)):4*i));
end

% caculate max
% current:
wtw_shr_max = zeros(48,2);
for i = 1:48;
        wtw_shr_max(i,1) = max(F_wtw_shr_cur(1+(4*(i-1)):4*i));
end
% low-carbon:
for i = 1:48;
        wtw_shr_max(i,2) = max(F_wtw_shr_lo(1+(4*(i-1)):4*i));
end


% switch columns to arrange average WTW shares by powertrain and MES
% current:
wtw_shr_avg_pwrtrn = zeros(48,2); 
for i = 1:8;
    for j = 0:5;    
        wtw_shr_avg_pwrtrn(i+(j*8),1) = wtw_shr_avg(1+6*(i-1)+j,1);
    end
end

% low-carbon:
for i = 1:8;
    for j = 0:5;
        wtw_shr_avg_pwrtrn(i+(j*8),2) = wtw_shr_avg(1+6*(i-1)+j,2);
    end
end


% same for min
% current:
wtw_shr_min_pwrtrn = zeros(48,2); 
for i = 1:8;
    for j = 0:5;    
        wtw_shr_min_pwrtrn(i+(j*8),1) = wtw_shr_min(1+6*(i-1)+j,1);
    end
end

% low-carbon:
for i = 1:8;
    for j = 0:5;
        wtw_shr_min_pwrtrn(i+(j*8),2) = wtw_shr_min(1+6*(i-1)+j,2);
    end
end


% same for max
% current:
wtw_shr_max_pwrtrn = zeros(48,2); 
for i = 1:8;
    for j = 0:5;    
        wtw_shr_max_pwrtrn(i+(j*8),1) = wtw_shr_max(1+6*(i-1)+j,1);
    end
end

% low-carbon:
for i = 1:8;
    for j = 0:5;
        wtw_shr_max_pwrtrn(i+(j*8),2) = wtw_shr_max(1+6*(i-1)+j,2);
    end
end


figure
% Plot means as bar
shrlabels = ({'No strategy','Lightweighting','Recycling','Remanufacturing','Downsizing','More intensive use','All strategies','All but lightweighting'});
y = 1:8;

subplot(6,2,1)
barshr1 = barh([1-wtw_shr_avg_pwrtrn(1:8,1) wtw_shr_avg_pwrtrn(1:8,1)], 'stacked');
title('a) ICEV-g, current')
set(gca,'yticklabel', shrlabels, 'YGrid','on','XGrid','on')
barshr1(1).FaceColor = ass_col;
barshr1(2).FaceColor =  wtw_col;
barshr1(1).EdgeAlpha = 0;
barshr1(2).EdgeAlpha = 0;
set(gca, 'YDir', 'reverse')
alpha(.7)

% Plot min and max as whiskers around mean
hold on
err_neg = [wtw_shr_max_pwrtrn(1:8,1)-wtw_shr_avg_pwrtrn(1:8,1)]; % x and y error, negative and positive 
err_pos = [wtw_shr_avg_pwrtrn(1:8,1)-wtw_shr_min_pwrtrn(1:8,1)];
err = errorbar(1-wtw_shr_avg_pwrtrn(1:8,1), y, err_neg, err_pos,'horizontal');
err.LineStyle = 'none';  
err.Color = [0 0 0];
err.CapSize = 4;


subplot(6,2,3)
barshr1 = barh([1-wtw_shr_avg_pwrtrn(9:16,1) wtw_shr_avg_pwrtrn(9:16,1)], 'stacked');
title('b) ICEV-d, current')
set(gca,'yticklabel', shrlabels, 'YGrid','on','XGrid','on')
barshr1(1).FaceColor = ass_col;
barshr1(2).FaceColor =  wtw_col;
barshr1(1).EdgeAlpha = 0;
barshr1(2).EdgeAlpha = 0;
set(gca, 'YDir', 'reverse')
alpha(.7)

hold on
err_neg = [wtw_shr_max_pwrtrn(9:16,1)-wtw_shr_avg_pwrtrn(9:16,1)]; % x and y error, negative and positive 
err_pos = [wtw_shr_avg_pwrtrn(9:16,1)-wtw_shr_min_pwrtrn(9:16,1)];
err = errorbar(1-wtw_shr_avg_pwrtrn(9:16,1), y, err_neg, err_pos,'horizontal');
err.LineStyle = 'none';  
err.Color = [0 0 0];
err.CapSize = 4;


subplot(6,2,5)
barshr1 = barh([1-wtw_shr_avg_pwrtrn(17:24,1) wtw_shr_avg_pwrtrn(17:24,1)], 'stacked');
title('c) HEV, current')
set(gca,'yticklabel', shrlabels, 'YGrid','on','XGrid','on')
barshr1(1).FaceColor = ass_col;
barshr1(2).FaceColor =  wtw_col;
barshr1(1).EdgeAlpha = 0;
barshr1(2).EdgeAlpha = 0;
set(gca, 'YDir', 'reverse')
alpha(.7)

hold on
err_neg = [wtw_shr_max_pwrtrn(17:24,1)-wtw_shr_avg_pwrtrn(17:24,1)]; % x and y error, negative and positive 
err_pos = [wtw_shr_avg_pwrtrn(17:24,1)-wtw_shr_min_pwrtrn(17:24,1)];
err = errorbar(1-wtw_shr_avg_pwrtrn(17:24,1), y, err_neg, err_pos,'horizontal');
err.LineStyle = 'none';  
err.Color = [0 0 0];
err.CapSize = 4;


subplot(6,2,7)
barshr1 = barh([1-wtw_shr_avg_pwrtrn(25:32,1) wtw_shr_avg_pwrtrn(25:32,1)], 'stacked');
title('d) PHEV, current')
set(gca,'yticklabel', shrlabels, 'YGrid','on','XGrid','on')
barshr1(1).FaceColor = ass_col;
barshr1(2).FaceColor =  wtw_col;
barshr1(1).EdgeAlpha = 0;
barshr1(2).EdgeAlpha = 0;
set(gca, 'YDir', 'reverse')
alpha(.7)

hold on
err_neg = [wtw_shr_max_pwrtrn(25:32,1)-wtw_shr_avg_pwrtrn(25:32,1)]; % x and y error, negative and positive 
err_pos = [wtw_shr_avg_pwrtrn(25:32,1)-wtw_shr_min_pwrtrn(25:32,1)];
err = errorbar(1-wtw_shr_avg_pwrtrn(25:32,1), y, err_neg, err_pos,'horizontal');
err.LineStyle = 'none';  
err.Color = [0 0 0];
err.CapSize = 4;


subplot(6,2,9)
barshr1 = barh([1-wtw_shr_avg_pwrtrn(33:40,1) wtw_shr_avg_pwrtrn(33:40,1)], 'stacked');
title('e) BEV, current')
set(gca,'yticklabel', shrlabels, 'YGrid','on','XGrid','on')
barshr1(1).FaceColor = ass_col;
barshr1(2).FaceColor =  wtw_col;
barshr1(1).EdgeAlpha = 0;
barshr1(2).EdgeAlpha = 0;
set(gca, 'YDir', 'reverse')
alpha(.7)

hold on
err_neg = [wtw_shr_max_pwrtrn(33:40,1)-wtw_shr_avg_pwrtrn(33:40,1)]; % x and y error, negative and positive 
err_pos = [wtw_shr_avg_pwrtrn(33:40,1)-wtw_shr_min_pwrtrn(33:40,1)];
err = errorbar(1-wtw_shr_avg_pwrtrn(33:40,1), y, err_neg, err_pos,'horizontal');
err.LineStyle = 'none';  
err.Color = [0 0 0];
err.CapSize = 4;


subplot(6,2,11)
barshr1 = barh([1-wtw_shr_avg_pwrtrn(41:48,1) wtw_shr_avg_pwrtrn(41:48,1)], 'stacked');
title('f) HFCEV, current')
set(gca,'yticklabel', shrlabels, 'YGrid','on','XGrid','on')
barshr1(1).FaceColor = ass_col;
barshr1(2).FaceColor =  wtw_col;
barshr1(1).EdgeAlpha = 0;
barshr1(2).EdgeAlpha = 0;
set(gca, 'YDir', 'reverse')
alpha(.7)

hold on
err_neg = [wtw_shr_max_pwrtrn(41:48,1)-wtw_shr_avg_pwrtrn(41:48,1)]; % x and y error, negative and positive 
err_pos = [wtw_shr_avg_pwrtrn(41:48,1)-wtw_shr_min_pwrtrn(41:48,1)];
err = errorbar(1-wtw_shr_avg_pwrtrn(41:48,1), y, err_neg, err_pos,'horizontal');
err.LineStyle = 'none';  
err.Color = [0 0 0];
err.CapSize = 4;


subplot(6,2,2)
barshr1 = barh([1-wtw_shr_avg_pwrtrn(1:8,2) wtw_shr_avg_pwrtrn(1:8,2)], 'stacked');
title('g) ICEV-g, low-carbon')
set(gca,'yticklabel', [{}], 'YGrid','on','XGrid','on')
barshr1(1).FaceColor = ass_col;
barshr1(2).FaceColor =  wtw_col;
barshr1(1).EdgeAlpha = 0;
barshr1(2).EdgeAlpha = 0;
set(gca, 'YDir', 'reverse')
alpha(.7)

hold on
err_neg = [wtw_shr_max_pwrtrn(1:8,2)-wtw_shr_avg_pwrtrn(1:8,2)]; % x and y error, negative and positive 
err_pos = [wtw_shr_avg_pwrtrn(1:8,2)-wtw_shr_min_pwrtrn(1:8,2)];
err = errorbar(1-wtw_shr_avg_pwrtrn(1:8,2), y, err_neg, err_pos,'horizontal');
err.LineStyle = 'none';  
err.Color = [0 0 0];
err.CapSize = 4;


subplot(6,2,4)
barshr1 = barh([1-wtw_shr_avg_pwrtrn(9:16,2) wtw_shr_avg_pwrtrn(9:16,2)], 'stacked');
title('h) ICEV-d, low-carbon')
set(gca,'yticklabel', [{}], 'YGrid','on','XGrid','on')
barshr1(1).FaceColor = ass_col;
barshr1(2).FaceColor =  wtw_col;
barshr1(1).EdgeAlpha = 0;
barshr1(2).EdgeAlpha = 0;
set(gca, 'YDir', 'reverse')
alpha(.7)

hold on
err_neg = [wtw_shr_max_pwrtrn(9:16,2)-wtw_shr_avg_pwrtrn(9:16,2)]; % x and y error, negative and positive 
err_pos = [wtw_shr_avg_pwrtrn(9:16,2)-wtw_shr_min_pwrtrn(9:16,2)];
err = errorbar(1-wtw_shr_avg_pwrtrn(9:16,2), y, err_neg, err_pos,'horizontal');
err.LineStyle = 'none';  
err.Color = [0 0 0];
err.CapSize = 4;


subplot(6,2,6)
barshr1 = barh([1-wtw_shr_avg_pwrtrn(17:24,2) wtw_shr_avg_pwrtrn(17:24,2)], 'stacked');
title('i) HEV, low-carbon')
set(gca,'yticklabel', [{}], 'YGrid','on','XGrid','on')
barshr1(1).FaceColor = ass_col;
barshr1(2).FaceColor =  wtw_col;
barshr1(1).EdgeAlpha = 0;
barshr1(2).EdgeAlpha = 0;
set(gca, 'YDir', 'reverse')
alpha(.7)

hold on
err_neg = [wtw_shr_max_pwrtrn(17:24,2)-wtw_shr_avg_pwrtrn(17:24,2)]; % x and y error, negative and positive 
err_pos = [wtw_shr_avg_pwrtrn(17:24,2)-wtw_shr_min_pwrtrn(17:24,2)];
err = errorbar(1-wtw_shr_avg_pwrtrn(17:24,2), y, err_neg, err_pos,'horizontal');
err.LineStyle = 'none';  
err.Color = [0 0 0];
err.CapSize = 4;


subplot(6,2,8)
barshr1 = barh([1-wtw_shr_avg_pwrtrn(25:32,2) wtw_shr_avg_pwrtrn(25:32,2)], 'stacked');
title('j) PHEV, low-carbon')
set(gca,'yticklabel', [{}], 'YGrid','on','XGrid','on')
barshr1(1).FaceColor = ass_col;
barshr1(2).FaceColor =  wtw_col;
barshr1(1).EdgeAlpha = 0;
barshr1(2).EdgeAlpha = 0;
set(gca, 'YDir', 'reverse')
alpha(.7)

hold on
err_neg = [wtw_shr_max_pwrtrn(25:32,2)-wtw_shr_avg_pwrtrn(25:32,2)]; % x and y error, negative and positive 
err_pos = [wtw_shr_avg_pwrtrn(25:32,2)-wtw_shr_min_pwrtrn(25:32,2)];
err = errorbar(1-wtw_shr_avg_pwrtrn(25:32,2), y, err_neg, err_pos,'horizontal');
err.LineStyle = 'none';  
err.Color = [0 0 0];
err.CapSize = 4;


subplot(6,2,10)
barshr1 = barh([1-wtw_shr_avg_pwrtrn(33:40,2) wtw_shr_avg_pwrtrn(33:40,2)], 'stacked');
title('k) BEV, low-carbon')
set(gca,'yticklabel', [{}], 'YGrid','on','XGrid','on')
barshr1(1).FaceColor = ass_col;
barshr1(2).FaceColor =  wtw_col;
barshr1(1).EdgeAlpha = 0;
barshr1(2).EdgeAlpha = 0;
set(gca, 'YDir', 'reverse')
alpha(.7)

hold on
err_neg = [wtw_shr_max_pwrtrn(33:40,2)-wtw_shr_avg_pwrtrn(33:40,2)]; % x and y error, negative and positive 
err_pos = [wtw_shr_avg_pwrtrn(33:40,2)-wtw_shr_min_pwrtrn(33:40,2)];
err = errorbar(1-wtw_shr_avg_pwrtrn(33:40,2), y, err_neg, err_pos,'horizontal');
err.LineStyle = 'none';  
err.Color = [0 0 0];
err.CapSize = 4;


subplot(6,2,12)
barshr1 = barh([1-wtw_shr_avg_pwrtrn(41:48,2) wtw_shr_avg_pwrtrn(41:48,2)], 'stacked');
title('l) HFCEV, low-carbon')
set(gca,'yticklabel', [{}], 'YGrid','on','XGrid','on')
barshr1(1).FaceColor = ass_col;
barshr1(2).FaceColor =  wtw_col;
barshr1(1).EdgeAlpha = 0;
barshr1(2).EdgeAlpha = 0;
set(gca, 'YDir', 'reverse')
alpha(.7)
legend({'Vehicle supply chain', 'Energy cycle'}, 'Location', 'best', 'Orientation', 'horizontal','AutoUpdate','off')

hold on
err_neg = [wtw_shr_max_pwrtrn(41:48,2)-wtw_shr_avg_pwrtrn(41:48,2)]; % x and y error, negative and positive 
err_pos = [wtw_shr_avg_pwrtrn(41:48,2)-wtw_shr_min_pwrtrn(41:48,2)];
err = errorbar(1-wtw_shr_avg_pwrtrn(41:48,2), y, err_neg, err_pos,'horizontal');
err.LineStyle = 'none';  
err.Color = [0 0 0];
err.CapSize = 4;


%% Print highest and lowest WTW share to command window:
% Current:
[min_wtw_shr_cur,idx_min_wtw_shr_cur] = min(F_wtw_shr_cur); % find smallest value and its index
all_x_mat(idx_min_wtw_shr_cur+1,1); % find corresponding label
fprintf('\nCurrent:\n')
fprintf('- With %.0f%%, %s has the highest Vehicle supply chain share.\n',100-min_wtw_shr_cur*100,all_x_mat{idx_min_wtw_shr_cur+1,1})

[max_wtw_shr_cur,idx_max_wtw_shr_cur] = max(F_wtw_shr_cur);
all_x_mat(idx_max_wtw_shr_cur+1,1);
fprintf('- With %.0f%%, %s has the lowest Vehicle supply chain share.\n',100-max_wtw_shr_cur*100,all_x_mat{idx_max_wtw_shr_cur+1,1})

% Low-carbon:
[min_wtw_shr_lo,idx_min_wtw_shr_lo] = min(F_wtw_shr_lo); 
all_x_mat(idx_min_wtw_shr_lo+1,1); 
fprintf('\nLow-carbon:\n')
fprintf('- With %.0f%%, %s has the highest Vehicle supply chain share.\n',100-min_wtw_shr_lo*100,all_x_mat{idx_min_wtw_shr_lo+1,1})

[max_wtw_shr_lo,idx_max_wtw_shr_lo] = max(F_wtw_shr_lo);
all_x_mat(idx_max_wtw_shr_lo+1,1);
fprintf('- With %.0f%%, %s has the lowest Vehicle supply chain share.\n',100-max_wtw_shr_lo*100,all_x_mat{idx_max_wtw_shr_lo+1,1})


%% Plot CEROI (Figure 5) 
% calculate CEROI (kg CO2 saved in Fuel cycle per 1
% additional kg CO2 during vehicle production):

F_wtw_cur = F_wtw(:,1);
F_wtw_lo = F_wtw(:,2);
F_veh_prd = F_mat + F_ass;
F_veh_prd_cur = F_mat(:,1) + F_ass(:,1);
F_veh_prd_lo = F_mat(:,2) + F_ass(:,2);
F_tot_cur = F_tot(:,1);
F_tot_lo = F_tot(:,2);

% Initialize deltas as zeros
% for 'LW only'
F_wtw_cur_dlt = zeros(1,24);
F_wtw_lo_dlt = zeros(1,24);
F_veh_prd_cur_dlt = zeros(1,24);
F_veh_prd_lo_dlt = zeros(1,24);
F_tot_cur_dlt = zeros(1,24);
F_tot_lo_dlt = zeros(1,24);

% for 'all MESs'
F_wtw_cur_dlt_abLW = zeros(1,24);
F_wtw_lo_dlt_abLW = zeros(1,24);
F_veh_prd_cur_dlt_abLW = zeros(1,24);
F_veh_prd_lo_dlt_abLW = zeros(1,24);
F_tot_cur_dlt_abLW = zeros(1,24);
F_tot_lo_dlt_abLW = zeros(1,24);

% Fill deltas for 'LW only'
F_wtw_cur_dlt = F_wtw_cur(1:24) - F_wtw_cur(25:48);
F_wtw_lo_dlt = F_wtw_lo(1:24) - F_wtw_lo(25:48);
F_veh_prd_cur_dlt = F_veh_prd_cur(1:24) - F_veh_prd_cur(25:48);
F_veh_prd_lo_dlt = F_veh_prd_lo(1:24) - F_veh_prd_lo(25:48);

% Fill deltas for 'all MESs'
F_wtw_cur_dlt_abLW = F_wtw_cur(169:192) - F_wtw_cur(145:168);
F_wtw_lo_dlt_abLW = F_wtw_lo(169:192) - F_wtw_lo(145:168);
F_veh_prd_cur_dlt_abLW = F_veh_prd_cur(169:192) - F_veh_prd_cur(145:168);
F_veh_prd_lo_dlt_abLW = F_veh_prd_lo(169:192) - F_veh_prd_lo(145:168);

labels_CEROI = [{'ICEV-g', 'ICEV-d', 'HEV', 'PHEV', 'BEV', 'HFCEV'}];

CEROI_cur = -F_wtw_cur_dlt ./ F_veh_prd_cur_dlt; % calculate CEROI
CEROI_lo = -F_wtw_lo_dlt ./ F_veh_prd_lo_dlt;
CEROI_cur_abLW = -F_wtw_cur_dlt_abLW ./ F_veh_prd_cur_dlt_abLW; 
CEROI_lo_abLW = -F_wtw_lo_dlt_abLW ./ F_veh_prd_lo_dlt_abLW;

unfav_color = [.84 .37 .00]; % red color 
fav_color = [.00 .62 .45]; % green color

figure
subplot(1,2,1)
boxplot([CEROI_cur(1:4),CEROI_cur(5:8),CEROI_cur(9:12),CEROI_cur(13:16), CEROI_cur(17:20),CEROI_cur(21:24)],'orientation','horizontal');
title('a) Current')%, LW only')
line([1 1],[25 0],'Color',unfav_color);
set(gca, 'YDir', 'reverse')
set(gca,'YTickLabel',labels_CEROI,'YTickLabelRotation',0,'FontSize',10,'YGrid','on','XGrid','on')
yticks(1:1:length(labels_CEROI))
alpha(.7)
xlim([0 16.5])


subplot(1,2,2)
boxplot([CEROI_lo(1:4),CEROI_lo(5:8),CEROI_lo(9:12),CEROI_lo(13:16), CEROI_lo(17:20),CEROI_lo(21:24)],'orientation','horizontal');
title('b) Low-carbon')%, LW only')
line([1 1],[25 0],'Color',unfav_color);
set(gca, 'YDir', 'reverse')
set(gca,'YTickLabel',[{ }],'YTickLabelRotation',0,'FontSize',10,'YGrid','on','XGrid','on')
yticks(1:1:length(labels_CEROI))
alpha(.7)
xlim([0 16.5])
text(1.25, 5.9,'CEROI = 1:1','Color',unfav_color,'FontSize',9)


%% Scatter plot, all, current, absolute (Figure S2)
% loop to calcuate differences between baseline and MES emissions 
% for each power train all four segments for base vehicle and 7 strategies
F_wtw_cur_scatter = zeros(8,24);
F_veh_prd_cur_scatter = zeros(8,24);

for j=1:8
    for k=1:24
        F_wtw_cur_scatter(j,k) = F_wtw_cur([(j-1) * 24 + 1 + (k-1)]);
        F_veh_prd_cur_scatter(j,k) = F_veh_prd_cur([(j-1) * 24 + 1 + (k-1)]);
    end
end

% scatter plot labels
labels_swap = labels(1:24);

% scatter plot colors
scattercol = {[0 0 0]; [.94 .89 .26]; [.00 .45 .70]; [.00 .62 .45]; [.84 .37 .00]; [.80 .47 .65]; [.75 .75 .75]; [.4 .4 .4 ]};
scattercol = repmat(scattercol,1,24);

marker_MES = 30;
sz = {45; marker_MES; marker_MES; marker_MES; marker_MES; marker_MES; marker_MES; marker_MES}; 
% marker size, baseline has a larger marker
sz = repmat(sz,1,24);

figure
for j = 1:8;
    for k = 1:24;
        subplot_scatter(k) = subplot(6,4,k);
        scatter1(j) = scatter(F_veh_prd_cur_scatter(j,k)/1000, F_wtw_cur_scatter(j,k)/1000, sz{j,k}, scattercol{j,k},'filled');
        title(labels{k});
        hold on     
        alpha(.7)
        grid on    
    end
end


legend1 = legend(subplot(6,4,3),{'No strategy','Lightweighting','Recycling','Remanufacturing','Downsizing','More intensive use','All strategies','All but lightweighting'},'Location','north', 'Orientation','horizontal','AutoUpdate','off');
legend1.NumColumns = 4;
hold off

for j = 1:8;
    for k = 1:24;
        subplot(6,4,k)
        line([F_veh_prd_cur_scatter(1,k)/1000 F_veh_prd_cur_scatter(j,k)/1000],[F_wtw_cur_scatter(1,k)/1000 F_wtw_cur_scatter(j,k)/1000],'Color',[.5 .5 .5]); % 'LineStyle','--',
        hold on
    end
end

axes1 = axes('visible','off');
xlabel('Vehicle supply chain emissions, t CO_{2}e person^{-1}');
ylabel('Energy cycle emissions, t CO_{2}e person^{-1}');
axes1(1).YLabel.Visible='on';
axes1(1).XLabel.Visible='on';


%% Scatter plot, all, low-carbon, absolute (Figure S3) 
F_wtw_lo_scatter = zeros(8,24);
F_veh_prd_lo_scatter = zeros(8,24);

for j=1:8
    for k=1:24
        F_wtw_lo_scatter(j,k) = F_wtw_lo([(j-1) * 24 + 1 + (k-1)]);
        F_veh_prd_lo_scatter(j,k) = F_veh_prd_lo([(j-1) * 24 + 1 + (k-1)]);
    end
end

% scatter plot labels
labels_swap = labels(1:24);

% scatter plot colors
scattercol = {[0 0 0]; [.94 .89 .26]; [.00 .45 .70]; [.00 .62 .45]; [.84 .37 .00]; [.80 .47 .65]; [.75 .75 .75]; [.4 .4 .4 ]};
scattercol = repmat(scattercol,1,24);

marker_MES = 30;
sz = {45; marker_MES; marker_MES; marker_MES; marker_MES; marker_MES; marker_MES; marker_MES}; % marker size
sz = repmat(sz,1,24);

figure
for j = 1:8;
    for k = 1:24;
        subplot_scatter(k) = subplot(6,4,k);
        scatter1(j) = scatter(F_veh_prd_lo_scatter(j,k)/1000, F_wtw_lo_scatter(j,k)/1000, sz{j,k}, scattercol{j,k},'filled');
        title(labels_swap{k});
        hold on     
        alpha(.7)
        grid on    
    end
end

legend1 = legend(subplot(6,4,3),{'No strategy','Lightweighting','Recycling','Remanufacturing','Downsizing','More intensive use','All strategies','All but lightweighting'},'Location','north', 'Orientation','horizontal','AutoUpdate','off');
legend1.NumColumns = 4;

hold off


for j = 1:8;
    for k = 1:24;
        subplot(6,4,k)
        line([F_veh_prd_lo_scatter(1,k)/1000 F_veh_prd_lo_scatter(j,k)/1000],[F_wtw_lo_scatter(1,k)/1000 F_wtw_lo_scatter(j,k)/1000],'Color',[.5 .5 .5]); % 'LineStyle','--',
        hold on 
    end
end

axes1 = axes('visible','off');
xlabel('Vehicle supply chain emissions, t CO_{2}e person^{-1}');
ylabel('Energy cycle emissions, t CO_{2}e person^{-1}');
axes1(1).YLabel.Visible='on';
axes1(1).XLabel.Visible='on';


%% Scatter plot, all, current, relative (Figure S1)
F_wtw_cur_scatter = zeros(8,24);
F_veh_prd_cur_scatter = zeros(8,24);

for j=1:8
    for k=1:24
        F_wtw_cur_scatter(j,k) = F_wtw_cur([(j-1) * 24 + 1 + (k-1)]);
        F_veh_prd_cur_scatter(j,k) = F_veh_prd_cur([(j-1) * 24 + 1 + (k-1)]);
    end
end

% scatter plot labels
labels_swap = labels(1:24);

% scatter plot colors
scattercol = {[0 0 0]; [.94 .89 .26]; [.00 .45 .70]; [.00 .62 .45]; [.84 .37 .00]; [.80 .47 .65]; [.75 .75 .75]; [.4 .4 .4 ]};
scattercol = repmat(scattercol,1,24);

marker_MES = 30;
sz = {45; marker_MES; marker_MES; marker_MES; marker_MES; marker_MES; marker_MES; marker_MES}; 
% marker size, baseline has a larger marker
sz = repmat(sz,1,24);

figure
for j = 1:8;
    for k = 1:24;
        subplot_scatter(k) = subplot(6,4,k);
        scatter1(j) = scatter(F_veh_prd_cur_scatter(j,k)/F_veh_prd_cur_scatter(1,k), F_wtw_cur_scatter(j,k)/F_wtw_cur_scatter(1,k), sz{j,k}, scattercol{j,k},'filled');
        title(labels{k});
        hold on     
        alpha(.7)
        grid on
        xlim([.25 1.35])
        xticks([.25 .5 .75 1 1.25])
        ylim([.4 1.1])
        yticks([.5 .75 1])
    end
end


legend1 = legend(subplot(6,4,3),{'No strategy','Lightweighting','Recycling','Remanufacturing','Downsizing','More intensive use','All strategies','All but lightweighting'},'Location','north', 'Orientation','horizontal','AutoUpdate','off');
legend1.NumColumns = 4;
hold off

for j = 1:8;
    for k = 1:24;
        subplot(6,4,k)
        line([1 F_veh_prd_cur_scatter(j,k)/F_veh_prd_cur_scatter(1,k)],[1 F_wtw_cur_scatter(j,k)/F_wtw_cur_scatter(1,k)],'Color',[.5 .5 .5]); % 'LineStyle','--',
        hold on
    end
end

axes1 = axes('visible','off');
xlabel('Relative vehicle supply chain emissions');
ylabel('Relative energy cycle emissions');
axes1(1).YLabel.Visible='on';
axes1(1).XLabel.Visible='on';


%% Scatter plot, all, low-carbon, relative (fig 9) 
F_wtw_lo_scatter = zeros(8,24);
F_veh_prd_lo_scatter = zeros(8,24);

for j=1:8
    for k=1:24
        F_wtw_lo_scatter(j,k) = F_wtw_lo([(j-1) * 24 + 1 + (k-1)]);
        F_veh_prd_lo_scatter(j,k) = F_veh_prd_lo([(j-1) * 24 + 1 + (k-1)]);
    end
end

% scatter plot labels
labels_swap = labels(1:24);

% scatter plot colors
scattercol = {[0 0 0]; [.94 .89 .26]; [.00 .45 .70]; [.00 .62 .45]; [.84 .37 .00]; [.80 .47 .65]; [.75 .75 .75]; [.4 .4 .4 ]};
scattercol = repmat(scattercol,1,24);

marker_MES = 30;
sz = {45; marker_MES; marker_MES; marker_MES; marker_MES; marker_MES; marker_MES; marker_MES}; % marker size
sz = repmat(sz,1,24);

figure
for j = 1:8;
    for k = 1:24;
        subplot_scatter(k) = subplot(6,4,k);
        scatter1(j) = scatter(F_veh_prd_lo_scatter(j,k)/F_veh_prd_lo_scatter(1,k), F_wtw_lo_scatter(j,k)/F_wtw_lo_scatter(1,k), sz{j,k}, scattercol{j,k},'filled');
        title(labels_swap{k});
        hold on     
        alpha(.7)
        grid on
        xlim([.25 1.35])
        xticks([.25 .5 .75 1 1.25])
        ylim([.4 1.1])
        yticks([.5 .75 1])
    end
end

legend1 = legend(subplot(6,4,3),{'No strategy','Lightweighting','Recycling','Remanufacturing','Downsizing','More intensive use','All strategies','All but lightweighting'},'Location','north', 'Orientation','horizontal','AutoUpdate','off');
legend1.NumColumns = 4;

hold off


for j = 1:8;
    for k = 1:24;
        subplot(6,4,k)
        line([1 F_veh_prd_lo_scatter(j,k)/F_veh_prd_lo_scatter(1,k)],[1 F_wtw_lo_scatter(j,k)/F_wtw_lo_scatter(1,k)],'Color',[.5 .5 .5]); % 'LineStyle','--',
        hold on
        line([],[])
    end
end

axes1 = axes('visible','off');
xlabel('Relative vehicle supply chain emissions');
ylabel('Relative energy cycle emissions');
axes1(1).YLabel.Visible='on';
axes1(1).XLabel.Visible='on';


%% Print relative fuel cycle, vehicle cycle, and total emission reductions
% Initialize as zeros
F_wtw_cur_red = zeros(168,1);
F_wtw_lo_red = zeros(168,1);
F_veh_prd_cur_red = zeros(168,1);
F_veh_prd_lo_red = zeros(168,1);
F_tot_cur_red = zeros(168,1);
F_tot_lo_red = zeros(168,1);

% calculate relative reduction potentials of all MESs
for a=0:6;
    F_wtw_cur_red(a*24+1:(a+1)*24) = (F_wtw_cur(1:24) - F_wtw_cur((a+1)*24+1:(a+2)*24)) ./ F_wtw_cur(1:24);
    F_wtw_lo_red(a*24+1:(a+1)*24) = (F_wtw_lo(1:24) - F_wtw_lo((a+1)*24+1:(a+2)*24)) ./ F_wtw_lo(1:24);
    F_veh_prd_cur_red(a*24+1:(a+1)*24) = (F_veh_prd_cur(1:24) - F_veh_prd_cur((a+1)*24+1:(a+2)*24)) ./ F_veh_prd_cur(1:24);
    F_veh_prd_lo_red(a*24+1:(a+1)*24) = (F_veh_prd_lo(1:24) - F_veh_prd_lo((a+1)*24+1:(a+2)*24)) ./ F_veh_prd_lo(1:24);
    F_tot_cur_red(a*24+1:(a+1)*24) = (F_tot_cur(1:24) - F_tot_cur((a+1)*24+1:(a+2)*24)) ./ F_tot_cur(1:24);
    F_tot_lo_red(a*24+1:(a+1)*24) = (F_tot_lo(1:24) - F_tot_lo((a+1)*24+1:(a+2)*24)) ./ F_tot_lo(1:24);
end

% Sort current reduction potentials:
% Total:
[F_tot_cur_red_sort,idx_tot_cur_red_sort]=sort(F_tot_cur_red,'ascend'); 
labels_tot_cur_red=[all_x_mat(26:193,1)];
labels_tot_cur_red=labels_tot_cur_red(idx_tot_cur_red_sort);

% Fuel cycle:
[F_wtw_cur_red_sort,idx_wtw_cur_red_sort]=sort(F_wtw_cur_red,'ascend'); % sort current reduction potentials, and extract index
labels_wtw_cur_red=[all_x_mat(26:193,1)];
labels_wtw_cur_red=labels_wtw_cur_red(idx_wtw_cur_red_sort);

% Vehicle cycle:
[F_veh_prd_cur_red_sort,idx_veh_prd_cur_red_sort]=sort(F_veh_prd_cur_red,'ascend'); % sort current reduction potentials, and extract index
labels_veh_prd_cur_red=[all_x_mat(26:193,1)];
labels_veh_prd_cur_red=labels_veh_prd_cur_red(idx_veh_prd_cur_red_sort);

% Sort low-carbon reduction potentials:
% Total:
[F_tot_lo_red_sort,idx_tot_lo_red_sort]=sort(F_tot_lo_red,'ascend'); 
labels_tot_lo_red=[all_x_mat(26:193,1)];
labels_tot_lo_red=labels_tot_lo_red(idx_tot_lo_red_sort);

% Fuel cycle:
[F_wtw_lo_red_sort,idx_wtw_lo_red_sort]=sort(F_wtw_lo_red,'ascend'); % sort current reduction potentials, and extract index
labels_wtw_lo_red=[all_x_mat(26:193,1)];
labels_wtw_lo_red=labels_wtw_lo_red(idx_wtw_lo_red_sort);

% Vehicle cycle:
[F_veh_prd_lo_red_sort,idx_veh_prd_lo_red_sort]=sort(F_veh_prd_lo_red,'ascend'); % sort current reduction potentials, and extract index
labels_veh_prd_lo_red=[all_x_mat(26:193,1)];
labels_veh_prd_lo_red=labels_veh_prd_lo_red(idx_veh_prd_lo_red_sort);

toc
