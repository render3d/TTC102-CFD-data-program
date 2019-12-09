%% CFD Data Processing
% %% Header
% tic
% if loc == 'u'
%     cd 'C:\Users\ttvr3\OneDrive - Loughborough University\1. University\4. Part C\Modules\CFD\CW 1200 09.12.2019'
% elseif loc == 'h'
%     cd 'C:\Users\Vince\OneDrive - Loughborough University\1. University\4. Part C\Modules\CFD\CW 1200 09.12.2019'
% elseif loc == 'm'
%     cd 'C:\Users\vjren\OneDrive - Loughborough University\1. University\4. Part C\Modules\CFD\CW 1200 09.12.2019'
% end
close all

%% SETUP
%% MESH 1
% Data import
% H Probe
hztl_probe_m1 = importdata('H Probe1.csv') ;
Hztl_Probe_M1 = struct2cell(hztl_probe_m1) ;
h_M1 = Hztl_Probe_M1{1,1} ;
H_M1 = sortrows(h_M1,6,'ascend') ;

% V Probe
vert_probe_m1 = importdata('V Probe1.csv') ;
Vert_Probe_M1 = struct2cell(vert_probe_m1) ;
v_M1 = Vert_Probe_M1{1,1} ;
V_M1 = sortrows(v_M1,7,'ascend') ;

% C Probe
cntl_probe_m1 = importdata('C Probe1.csv') ;
Cntl_Probe_M1 = struct2cell(cntl_probe_m1) ;
c_M1 = Cntl_Probe_M1{1,1} ;
C_M1 = sortrows(c_M1,8,'ascend') ;

% Data partitioning
xx_H1 = H_M1(:,6) ;
pm_H1 = H_M1(:,5) ;
vm_H1 = H_M1(:,4) ;
vi_H1 = H_M1(:,1) ;
vj_H1 = H_M1(:,2) ;
% vk_H1 = H_M1(:,3) ;


xy_V1 = V_M1(:,7) ;
pm_V1 = V_M1(:,1) ;
vm_V1 = V_M1(:,2) ;
vi_V1 = V_M1(:,3) ;
vj_V1 = V_M1(:,4) ;
% vk_V2 = V_M1(:,5) ;

xz_C1 = C_M1(:,8) ;
pm_C1 = C_M1(:,5) ;
vm_C1 = C_M1(:,1) ;
vi_C1 = C_M1(:,2) ;
vj_C1 = C_M1(:,3) ;
% vk_C1 = C_M1(:,4) ;

%% MESH 2
% Data import
% H Probe
hztl_probe_m2 = importdata('H Probe2.csv') ;
Hztl_Probe_M2 = struct2cell(hztl_probe_m2) ;
h_M2 = Hztl_Probe_M2{1,1} ;
H_M2 = sortrows(h_M2,6,'ascend') ;

% V Probe
vert_probe_m2 = importdata('V Probe2.csv') ;
Vert_Probe_M2 = struct2cell(vert_probe_m2) ;
v_M2 = Vert_Probe_M2{1,1} ;
V_M2 = sortrows(v_M2,7,'ascend') ;

% C Probe
cntl_probe_m2 = importdata('C Probe2.csv') ;
Cntl_Probe_M2 = struct2cell(cntl_probe_m2) ;
c_M2 = Cntl_Probe_M2{1,1} ;
C_M2 = sortrows(c_M2,8,'ascend') ;

% Data partitioning
xx_H2 = H_M2(:,6) ;
pm_H2 = H_M2(:,5) ;
vm_H2 = H_M2(:,4) ;
vi_H2 = H_M2(:,1) ;
vj_H2 = H_M2(:,2) ;
% vk_H2 = H_M2(:,3) ;

xy_V2 = V_M2(:,7) ;
pm_V2 = V_M2(:,1) ;
vm_V2 = V_M2(:,2) ;
vi_V2 = V_M2(:,3) ;
vj_V2 = V_M2(:,4) ;
% vk_V2 = V_M2(:,5) ;

xz_C2 = C_M2(:,8) ;
pm_C2 = C_M2(:,5) ;
vm_C2 = C_M2(:,1) ;
vi_C2 = C_M2(:,2) ;
vj_C2 = C_M2(:,3) ;
% vk_C2 = C_M2(:,4) ;

%% MESH 3
% Data import
% H Probe
hztl_probe_m3 = importdata('H Probe3.csv') ;
Hztl_Probe_M3 = struct2cell(hztl_probe_m3) ;
h_M3 = Hztl_Probe_M3{1,1} ;
H_M3 = sortrows(h_M3,6,'ascend') ;

% V Probe
vert_probe_m3 = importdata('V Probe3.csv') ;
Vert_Probe_M3 = struct2cell(vert_probe_m3) ;
v_M3 = Vert_Probe_M3{1,1} ;
V_M3 = sortrows(v_M3,7,'ascend') ;

% C Probe
cntl_probe_m3 = importdata('C Probe3.csv') ;
Cntl_Probe_M3 = struct2cell(cntl_probe_m3) ;
c_M3 = Cntl_Probe_M3{1,1} ;
C_M3 = sortrows(c_M3,8,'ascend') ;

% Data partitioning
xx_H3 = H_M3(:,6) ;
pm_H3 = H_M3(:,5) ;
vm_H3 = H_M3(:,4) ;
vi_H3 = H_M3(:,1) ;
vj_H3 = H_M3(:,2) ;
% vk_H3 = H_M3(:,3) ;

xy_V3 = V_M3(:,7) ;
pm_V3 = V_M3(:,1) ;
vm_V3 = V_M3(:,2) ;
vi_V3 = V_M3(:,3) ;
vj_V3 = V_M3(:,4) ;
% vk_V3 = V_M3(:,5) ;

xz_C3 = C_M3(:,8) ;
pm_C3 = C_M3(:,5) ;
vm_C3 = C_M3(:,1) ;
vi_C3 = C_M3(:,2) ;
vj_C3 = C_M3(:,3) ;
% vk_C3 = C_M3(:,4) ;

%% MESH 4
% Data import
% H Probe
hztl_probe_m4 = importdata('H Probe4.csv') ;
Hztl_Probe_M4 = struct2cell(hztl_probe_m4) ;
h_M4 = Hztl_Probe_M4{1,1} ;
H_M4 = sortrows(h_M4,6,'ascend') ;

% V Probe
vert_probe_m4 = importdata('V Probe4.csv') ;
Vert_Probe_M4 = struct2cell(vert_probe_m4) ;
v_M4 = Vert_Probe_M4{1,1} ;
V_M4 = sortrows(v_M4,7,'ascend') ;

% C Probe
cntl_probe_m4 = importdata('C Probe4.csv') ;
Cntl_Probe_M4 = struct2cell(cntl_probe_m4) ;
c_M4 = Cntl_Probe_M4{1,1} ;
C_M4 = sortrows(c_M4,8,'ascend') ;

% Data partitioning
xx_H4 = H_M4(:,6) ;
pm_H4 = H_M4(:,5) ;
vm_H4 = H_M4(:,4) ;
vi_H4 = H_M4(:,1) ;
vj_H4 = H_M4(:,2) ;
% vk_H6 = H_M4(:,3) ;

xy_V4 = V_M4(:,7) ;
pm_V4 = V_M4(:,1) ;
vm_V4 = V_M4(:,2) ;
vi_V4 = V_M4(:,3) ;
vj_V4 = V_M4(:,4) ;
% vk_V6 = V_M4(:,5) ;

xz_C4 = C_M4(:,8) ;
pm_C4 = C_M4(:,5) ;
vm_C4 = C_M4(:,4) ;
vi_C4 = C_M4(:,1) ;
vj_C4 = C_M4(:,2) ;
% vk_C4 = C_M4(:,3) ;

%% PLOTS
%% H Plots
figure(1) 
% % Pressure
% subplot(2,2,1)
% hold on
% grid on
% plot(xx_H4,pm_H4,'r+-')
% plot(xx_H3,pm_H3,'g+-')
% plot(xx_H2,pm_H2,'b+-')
% plot(xx_H1,pm_H1,'k+-')
% legend('Mesh 4','Mesh 3','Mesh 2','Mesh 1','location','best')
% title('Pressure - H Probe')
% ylabel('Pressure (Pa)')
% xlabel('Flow Field Co-ordinate (x)')
% xlim([-0.1,0.2])

% Velocity Magnitude
% subplot(2,2,2)
hold on
grid on
plot(xx_H4,vm_H4,'r+-')
plot(xx_H3,vm_H3,'g+-')
plot(xx_H2,vm_H2,'b+-')
plot(xx_H1,vm_H1,'k+-')
legend('Mesh 4','Mesh 3','Mesh 2','Mesh 1','location','best')
title('Velocity Magnitude - H Probe')
ylabel('Velocity (m/s)')
xlabel('Flow Field Co-ordinate (x)')
xlim([-0.1,0.2])

% % Velocity (i)
% subplot(2,2,3)
% hold on
% grid on
% plot(xx_H4,vi_H4,'r+-')
% plot(xx_H3,vi_H3,'g+-')
% plot(xx_H2,vi_H2,'b+-')
% plot(xx_H1,vi_H1,'k+-')
% legend('Mesh 4','Mesh 3','Mesh 2','Mesh 1','location','best')
% title('Velocity (i) - H Probe')
% ylabel('Velocity (m/s)')
% xlabel('Flow Field Co-ordinate (x)')
% xlim([-0.1,0.2])
% 
% % Velocity (j)
% subplot(2,2,4)
% hold on
% grid on
% plot(xx_H4,vj_H4,'r+-')
% plot(xx_H3,vj_H3,'g+-')
% plot(xx_H2,vj_H2,'b+-')
% plot(xx_H1,vj_H1,'k+-')
% legend('Mesh 4','Mesh 3','Mesh 2','Mesh 1','location','best')
% title('Velocity (j) - H Probe')
% ylabel('Velocity (m/s)')
% xlabel('Flow Field Co-ordinate (x)')
% xlim([-0.1,0.2])

% % Velocity (k)
% subplot(3,2,5)
% hold on
% grid on
% plot(xx_H6,vk_H6,'r+-')
% plot(xx_H3,vk_H3,'g+-')
% plot(xx_H2,vk_H2,'b+-')
% plot(xx_H1,vk_H1,'k+-')
% legend('Mesh 4','Mesh 3','Mesh 2','Mesh 1','location','best')
% title('Velocity (k) - H Probe')
% ylabel('Velocity (m/s)')
% xlabel('Flow Field Co-ordinate (x)')

filename1 = 'Fig 1 Horizontal Probe Sensitivity' ;
savefig(filename1)

%% V Plots
figure(2)
% % Pressure
% subplot(2,2,1)
% hold on
% grid on
% plot(xy_V4,pm_V4,'r+-')
% plot(xy_V3,pm_V3,'g+-')
% plot(xy_V2,pm_V2,'b+-')
% plot(xy_V1,pm_V1,'k+-')
% legend('Mesh 4','Mesh 3','Mesh 2','Mesh 1','location','best')
% title('Pressure - V Probe')
% ylabel('Pressure (Pa)')
% xlabel('Flow Field Co-ordinate (y)')

% % Velocity Magnitude
% subplot(2,2,2)
hold on
grid on
plot(xy_V4,vm_V4,'r+-')
plot(xy_V3,vm_V3,'g+-')
plot(xy_V2,vm_V2,'b+-')
plot(xy_V1,vm_V1,'k+-')
legend('Mesh 4','Mesh 3','Mesh 2','Mesh 1','location','southeast')
title('Velocity Magnitude - V Probe')
ylabel('Velocity (m/s)')
xlabel('Flow Field Co-ordinate (y)')

% % Velocity (i)
% subplot(2,2,3)
% hold on
% grid on
% plot(xy_V4,vi_V4,'r+-')
% plot(xy_V3,vi_V3,'g+-')
% plot(xy_V2,vi_V2,'b+-')
% plot(xy_V1,vi_V1,'k+-')
% legend('Mesh 4','Mesh 3','Mesh 2','Mesh 1','location','southeast')
% title('Velocity (i) - V Probe')
% ylabel('Velocity (m/s)')
% xlabel('Flow Field Co-ordinate (y)')

% % Velocity (j)
% subplot(2,2,4)
% hold on
% grid on
% plot(xy_V4,vj_V4,'r+-')
% plot(xy_V3,vj_V3,'g+-')
% plot(xy_V2,vj_V2,'b+-')
% plot(xy_V1,vj_V1,'k+-')
% legend('Mesh 4','Mesh 3','Mesh 2','Mesh 1','location','best')
% title('Velocity (j) - V Probe')
% ylabel('Velocity (m/s)')
% xlabel('Flow Field Co-ordinate (y)')

% % Velocity (k)
% subplot(3,2,5)
% hold on
% grid on
% plot(xy_V6,vk_V6,'r+-')
% plot(xy_V3,vk_V3,'g+-')
% plot(xy_V2,vk_V2,'b+-')
% plot(xy_V1,vk_V1,'k+-')
% legend('Mesh 4','Mesh 3','Mesh 2','Mesh 1','location','best')
% title('Velocity (k) - V Probe')
% ylabel('Velocity (m/s)')
% xlabel('Flow Field Co-ordinate (y)')

filename2 = 'Fig 2 Vertical Probe Sensitivity' ;
savefig(filename2)

%% C Plots
figure(3)
% % Pressure
% subplot(2,2,1)
% hold on
% grid on
% plot(xz_C4,pm_C4,'r+-')
% plot(xz_C3,pm_C3,'g+-')
% plot(xz_C2,pm_C2,'b+-')
% plot(xz_C1,pm_C1,'k+-')
% legend('Mesh 4','Mesh 3','Mesh 2','Mesh 1','location','best')
% title('Pressure - C Probe')
% ylabel('Pressure (Pa)')
% xlabel('Flow Field Co-ordinate (z)')

% % Velocity Magnitude
% subplot(2,2,2)
hold on
grid on
plot(xz_C4,vm_C4,'r+-')
plot(xz_C3,vm_C3,'g+-')
plot(xz_C2,vm_C2,'b+-')
plot(xz_C1,vm_C1,'k+-')
legend('Mesh 4','Mesh 3','Mesh 2','Mesh 1','location','best')
title('Velocity Magnitude - C Probe')
ylabel('Velocity (m/s)')
xlabel('Flow Field Co-ordinate (z)')
ylim([38,42])
xlim([0,0.0215])

% % Velocity (i)
% subplot(2,2,3)
% hold on
% grid on
% plot(xz_C4,vi_C4,'r+-')
% plot(xz_C3,vi_C3,'g+-')
% plot(xz_C2,vi_C2,'b+-')
% plot(xz_C1,vi_C1,'k+-')
% legend('Mesh 4','Mesh 3','Mesh 2','Mesh 1','location','southeast')
% title('Velocity (i) - C Probe')
% ylabel('Velocity (m/s)')
% xlabel('Flow Field Co-ordinate (z)')

% % Velocity (j)
% subplot(2,2,4)
% hold on
% grid on
% plot(xz_C4,vj_C4,'r+-')
% plot(xz_C3,vj_C3,'g+-')
% plot(xz_C2,vj_C2,'b+-')
% plot(xz_C1,vj_C1,'k+-')
% legend('Mesh 4','Mesh 3','Mesh 2','Mesh 1','location','best')
% title('Velocity (j) - C Probe')
% ylabel('Velocity (m/s)')
% xlabel('Flow Field Co-ordinate (z)')

% % Velocity (k)
% subplot(3,2,5)
% hold on
% grid on
% plot(xy_V6,vk_V6,'r+-')
% plot(xy_V3,vk_V3,'g+-')
% plot(xy_V2,vk_V2,'b+-')
% plot(xy_V1,vk_V1,'k+-')
% legend('Mesh 4','Mesh 3','Mesh 2','Mesh 1','location','best')
% title('Velocity (k) - S Probe')
% ylabel('Velocity (m/s)')
% xlabel('Flow Field Co-ordinate (z)')

filename3 = 'Fig 3 Span Probe Sensitivity' ;
savefig(filename3)

%% UNSTEADY-EXPERIMENTAL DATA COMPARISON
clear
W = 0.02 ;

expm_raw = importdata('TTC102_Exp_profile_csv.csv') ;
expm_flt = struct2cell(expm_raw) ;
expm_dat = expm_flt{1,1} ;
pos_yWE = expm_dat(:,1) ;
pos_yE = pos_yWE*W ;
UvsUin = expm_dat(:,2) ;
Um_exp = 40*UvsUin ;

% edit data import here 
cfd_raw = importdata('E Probe.csv') ;
cfd_flt = struct2cell(cfd_raw) ;
cfd_dat = cfd_flt{1,1} ;
pos_yH = cfd_dat(:,8) ;
pos_yHS = pos_yH/W ;
Vm_cfd = cfd_dat(:,5) ;
VvsVin = Vm_cfd.*(1/40) ;

% close all
figure
subplot(2,1,1)
title('Comparison of Relative Velocity Magnitude Profiles of CFD & Experimental Data')
hold on
grid on
% yyaxis left
plot(VvsVin,pos_yHS,'b-')
% ylim([-0.045,0.045])
ylabel('Relative Flow Field Position y / W')
% yyaxis right
plot(UvsUin,pos_yWE,'g-')
ylim([-2.25,2.25])
% ylabel('Simulation Flow Field Position y / H')
legend('CFD Relative Velicity Magnitude, V/Vin','Experimental Velicity Magnitude, U/Uin','location','northwest')
xlim([0.5,1.5])
xlabel('Flow Velocity Magnitude Relative to Inlet')
set(gca,'XAxisLocation','origin')

% figure
subplot(2,1,2)
title('Comparison of Absolute Velocity Magnitude Profiles of CFD & Experimental Data')
hold on
grid on
plot(Vm_cfd,pos_yH,'b-')
plot(Um_exp,pos_yE,'g-')
legend('CFD Velicity Magnitude, Vm','Experimental Velicity Magnitude, Um','location','southwest')
ylim([-0.045,0.045])
ylabel('Flow Field Co-ordinate (y)')
% xlim([30,50])
xlabel('Absolute Flow Velocity Magnitude (m/s)')
ylabel('Absolute Flow Field Co-ordinate, y (m)')
set(gca,'XAxisLocation','origin')

filename4 = 'Fig 4 Sim vs Expm Flow Velocity Plots' ;
savefig(filename4)

%% MATLAB Workspace Data
filename5 = 'simdata' ;
save(filename5)

%% Footer
disp(' ')
toc