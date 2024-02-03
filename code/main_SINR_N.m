% This Matlab script can be used to generate Fig. 2 and Fig. 3 in the paper:
% R. Liu, M. Li, H. Luo, Q. Liu, and A. L. Swindlehurst, “Integrated sensing and communication with reconfigurable intelligent surfaces: Opportunities, applications, and future directions,” IEEE Wireless Commun., vol. 30, no. 1, pp. 50-57, Feb. 2023.
% Download this paper at: https://ieeexplore.ieee.org/document/10077119
% Last edited by Rang Liu (rangl2@uci.edu) in 2024-02-01

clear
clc
%%%% system settings
Prms.M = 16;  M = Prms.M; %%% number of transmit/receive antennas
Locations_clutter = [35 15]; %%% range-angle of clutter
Locations_target = [35 7];
Locations_RIS = [30,0];
Locations_user = [25 -7; 35 -7];
Locations = [Locations_clutter;Locations_RIS;Locations_user;Locations_target;0 0];
Prms.K = length(Locations_user(:,1));  K = Prms.K; %%% number of users
Prms.L = 20; L = Prms.L; %%% number of samples   
Prms.Q = length(Locations_clutter(:,1)); Q = Prms.Q; %%% number of clutter patches
Prms.sigmar2 = 10^(-11); sigmar2 = Prms.sigmar2; %%% radar noise
Prms.sigmac2 = 10^(-11); sigmac2 = Prms.sigmac2; %%% communication noise
Prms.Phi = pi/4; Phi = Prms.Phi;  %%% QPSK modulated
Prms.Nmax = 1000; Nmax = Prms.Nmax; %%% maximum iteration
Prms.res_th = 5e-4; %%%convergence tolerance
Prms.P = L*100; P = Prms.P; %%% total transmit power
Prms.sigma2 = 1; sigma2 = Prms.sigma2; %%% RCS
%%% channel settings
%%%% distances
dg = Locations_RIS(1);
drt = norm(Locations_target-Locations_RIS,2);
dt = norm(Locations_target,2);
theta_t = atan(Locations_target(2)/Locations_target(1));
if Locations_target(1) < Locations_RIS(1)
    theta_rt = pi - atan(Locations_target(2)/(Locations_RIS(1)-Locations_target(1)));
else
    theta_rt = atan(Locations_target(2)/(Locations_target(1)-Locations_RIS(1)));
end
drk = zeros(1,K);
dk = zeros(1,K);
theta_k = zeros(1,K);
theta_rk = zeros(1,K);
for k = 1:1:K
    drk(k) = norm(Locations_user(k,:)-Locations_RIS,2);
    dk(k) = norm(Locations_user(k,:),2);
    theta_k(k) = atan(Locations_user(k,2)/Locations_user(k,1));
    if Locations_user(k,1) < Locations_RIS(1)
        theta_rk(k) = -pi + atan(abs(Locations_user(k,2))/(Locations_RIS(1)-Locations_user(k,1)));
    else
        theta_rk(k) = atan(Locations_user(k,2)/(Locations_user(k,1)-Locations_RIS(1)));
    end
end
dc = zeros(1,Q);
drc = zeros(1,Q);
theta_c = zeros(1,Q);
theta_rc = zeros(1,Q);
for q = 1:1:Q
    dc(q) = norm(Locations_clutter(q,:),2);
    drc(q) = norm(Locations_clutter(q,:)-Locations_RIS,2);
    theta_c(q) = atan(Locations_clutter(q,2)/Locations_clutter(q,1));
    if Locations_clutter(q,1) < Locations_RIS(1)
        theta_rc(q) = pi - atan((Locations_clutter(q,2))/(Locations_RIS(1)-Locations_clutter(q,1)));
    else
        theta_rc(q) = atan((Locations_clutter(q,2))/(Locations_clutter(q,1)-Locations_RIS(1)));
    end
end
theta_RIS = 0;
%%%% path-loss
alpha_t = 3.2;
alpha_rt = 2.4;
alpha_k = 3.2;
alpha_rk = 2.4;
alpha_g = 2.2;

Nx = (4:2:16);
N_range = Nx.^2;

Channel.ht = sqrt(10^(-3)*dt^(-alpha_t))*exp(-1j*(0:1:M-1)*pi*sin(theta_t));%/sqrt(M);%
ht = Channel.ht;
hrt0 = sqrt(10^(-3)*drt^(-alpha_rt))*exp(-1j*(0:1:N_range(end)-1)*pi*sin(theta_rt));%
Channel.Hc = zeros(Q,M);
Hrc0 = zeros(Q,N_range(end));
for q = 1:1:Q
    Channel.Hc(q,:) = sqrt(10^(-3)*dc(q)^(-alpha_t))*exp(-1j*(0:1:M-1)*pi*sin(theta_c(q))); %
    Hrc0(q,:) = sqrt(10^(-3)*drc(q)^(-alpha_rt))*exp(-1j*(0:1:N_range(end)-1)*pi*sin(theta_rc(q))); %
end
Hc = Channel.Hc;
SNR = 10*ones(1,K);
Prms.gamma = sqrt(sigmac2*10.^(0.1*SNR));

SINR_my_CI = zeros(1,length(N_range));
SINR_wo_RIS_CI = zeros(1,length(N_range));

CG_my_t1 = zeros(1,length(N_range));
CG_my_t2 = zeros(1,length(N_range));
CG_my_t3 = zeros(1,length(N_range));
CG_my_t4 = zeros(1,length(N_range));
CG_my_sum = zeros(1,length(N_range));
CG_woRIS_t = zeros(1,length(N_range));

G0 = sqrt(10^(-3)*dg^(-alpha_g))*exp(-1j*(0:1:N_range(end)-1)'*pi*sin(pi))*exp(-1j*(0:1:M-1)*pi*sin(pi));
Hu = zeros(K,M);
Hru0 = zeros(K,N_range(end));
for k = 1:1:K
    Hu(k,:) = sqrt(10^(-3)*dk(k)^(-alpha_k))*exp(-1j*(0:1:M-1)*pi*sin(theta_k(k)));
    Hru0(k,:) = sqrt(10^(-3)*drk(k)^(-alpha_rk))*exp(-1j*(0:1:N_range(end)-1)*pi*sin(theta_rk(k)));%
end
Channel.Hu = Hu;

S = exp(1i*(Phi+2*Phi*randi([0,pi/Phi-1],K,L)));  

[x_wo_RIS,VSINR_wo_RIS] = get_x_woRIS_CI(Prms,Channel,S);
    SINR_wo_RIS_CI = SINR_wo_RIS_CI + VSINR_wo_RIS(end)*ones(1,length(N_range));
%%% BS-target/clutter-BS
Rt = ht'*ht*reshape(x_wo_RIS,M,L);
CG_woRIS_t = CG_woRIS_t + 10*log10(norm(Rt,'fro')^2/sigmar2); 

for N_index = 1:1:length(N_range)

    N = N_range(N_index)
    Prms.N = N;
    Channel.hrt = hrt0(1:N);
    hrt = Channel.hrt;
    Channel.Hrc = Hrc0(:,1:N);
    Hrc = Channel.Hrc;
    Channel.G = G0(1:N,:);
    G = Channel.G;
    Channel.Hru = Hru0(:,1:N);

    Ni = 10;
    XX_my = zeros(M*L,Ni);
    PPhi_my = zeros(N,Ni);
    Vsinr = zeros(1,Ni);
    for ii = 1:1:Ni
        [x_my,phi_my,VSINR_my] = get_x_phi_CI(Prms,Channel,S);
        XX_my(:,ii) = x_my;
        PPhi_my(:,ii) = phi_my;
        Vsinr(ii) = VSINR_my(end);
    end
    [~,ind] = find(Vsinr == max(Vsinr));
    x_my = XX_my(:,ind);
    phi_my = PPhi_my(:,ind);
    vs = Vsinr(ind);

    SINR_my_CI(N_index) = SINR_my_CI(N_index) + vs;
    X_my = reshape(x_my,M,L);
    %%% BS-target/clutter-BS
    Rt1 = ht'*ht*X_my;
    CG_my_t1(N_index) = CG_my_t1(N_index) + 10*log10(norm(Rt1,'fro')^2/sigmar2); 
    %%% BS-RIS-target-BS
    Rt2 = ht'*(hrt*diag(phi_my)*G)*X_my;
    CG_my_t2(N_index) = CG_my_t2(N_index) + 10*log10(norm(Rt2,'fro')^2/sigmar2); 
    %%% BS-target-RIS-BS
    Rt3 = (hrt*diag(phi_my)*G)'*ht*X_my;
    CG_my_t3(N_index) = CG_my_t3(N_index) + 10*log10(norm(Rt3,'fro')^2/sigmar2); 
    %%% BS-RIS-target-RIS-BS
    Rt4 = (hrt*diag(phi_my)*G)'*(hrt*diag(phi_my)*G)*X_my;
    CG_my_t4(N_index) = CG_my_t4(N_index) + 10*log10(norm(Rt4,'fro')^2/sigmar2); 

    CG_my_sum(N_index) = CG_my_sum(N_index) + norm(Rt1,'fro')^2/sigmar2 + norm(Rt2,'fro')^2/sigmar2 + norm(Rt3,'fro')^2/sigmar2 + norm(Rt4,'fro')^2/sigmar2;
 
end

CG_my_sum = 10*log10(CG_my_sum);

figure
plot(Nx,CG_my_t1,'-d','color',[0,0.35,0.7],'LineWidth',1.5)
hold on
plot(Nx,CG_my_t2,'-s','color',[0,0.35,0.7],'LineWidth',1.5)
plot(Nx,CG_my_t3,'-^','color',[0,0.35,0.7],'LineWidth',1.5)
plot(Nx,CG_my_t4,'-x','color',[0,0.35,0.7],'LineWidth',1.5)
plot(Nx,CG_my_sum,'-o','color',[0.8,0,0],'LineWidth',1.5)
plot(Nx,CG_woRIS_t,'-*','color',[0,0,0],'LineWidth',1.5)
hold off
grid on
xlabel('$$\sqrt{N}$$','Interpreter','latex');
ylabel('SNR gain (dB)');
legend('Path 1: BS-target-BS','Path 2: BS-RIS-target-BS','Path 3: BS-target-RIS-BS',...
    'Path 4: BS-RIS-target-RIS-BS','Sum of path 1-4', 'Benchmark, w/o RIS');
xlabel('$$\sqrt{N}$$','Interpreter','latex');

dd = 10;
xind = -dd+min(Locations(:,1)):0.2:max(Locations(:,1))+dd;
yind = -dd+min(Locations(:,2)):0.2:max(Locations(:,2))+dd;
Bp = zeros(length(xind),length(yind));
for i = 1:1:length(xind)
    x = xind(i);
    for j = 1:1:length(yind)
        y = yind(j);
        if x < 0 && y > 0
            theta1 = atan(-x/y)+pi/2;
            theta2 = pi-atan(y/(dg-x));
        elseif x >= 0 && y > 0
            theta1 = atan(y/x);
            if x < dg
                theta2 = pi-atan(y/(dg-x));
            else
                theta2 = atan(y/(x-dg));
            end
        elseif x < 0 && y <= 0
            theta1 = -(atan(x/y)+pi/2);
            theta2 = atan(-y/(dg-x))-pi;
        else
            theta1 = atan(y/x);
            if x < dg
                theta2 = atan(-y/(dg-x))-pi;
            else
                theta2 = atan(y/(x-dg));
            end
        end
        d1 = sqrt(x^2+y^2);
        d2 = sqrt((x-dg)^2+y^2);
        h1 = sqrt(10^(-3)*d1^(-alpha_t))*exp(-1j*(0:1:M-1)*pi*sin(theta1));
        h2 = sqrt(10^(-3)*d2^(-alpha_rt))*exp(-1j*(0:1:N-1)*pi*sin(theta2));%
        Bp(i,j) = 10*log10( norm((h1 + h2*diag(phi_my)*G)*reshape(x_my,M,L),2)^2 );
    end
end

figure;
pcolor(xind,yind,Bp.');
shading interp;
colorbar;colorbar;
colormap(jet);
hold on;
plot(0,0,'kd','markersize',9.5,'linewidth',1.5)
plot(dg,0,'ks','markersize',9.5,'linewidth',1.5)
plot(Locations_target(1),Locations_target(2),'kp','markersize',9.5,'linewidth',1.5)
for k = 1:1:K
    plot(Locations_user(k,1),Locations_user(k,2),'ko','markersize',9.5,'linewidth',1.5)
end
for q = 1:1:Q
    plot(Locations_clutter(q,1),Locations_clutter(q,2),'k^','markersize',9.5,'linewidth',1.5)
end
grid on
hold off
xlabel('{\it x}-axis (m)')
ylabel('{\it y}-axis (m)')

