% clear
% clc
%% Institude: HuaZhong University of Science and Technology
% 机构： 华中科技大学电气与电子工程学院
%% Written by Xu Shouyu
% 作者： 徐首彧
% 指导老师： 叶才勇
%% 1.参数

Rsb=0.16;                  % 定子槽底部半径
Rt=0.115;                  % 定子槽顶部半径
Rs=0.113;                  % 定子内半径
Rm=0.105;                  % 磁极半径
Rr=0.092;                  % 转子半径
boa=0.02655;                 % boa是槽开口角
bsa=0.0607;                  % bsa是槽宽角
ap=0.82;                     % 极弧系数

p=2;                     % 极对数
Ns=48;                   % 齿槽个数
theta_s=2*pi/Ns;         % 定子齿槽极距，两个齿槽之间的间距，与选取半径无关，有槽数有关
                         % 第一个定子槽中心角的位置设为theta_s/2
                         
a0=0;

u0=4*pi*1e-7;  % 真空磁导率
ur=1.05;       % PM的相对磁导率，解析模型的假设
Br=1;          % PM的剩磁/T，有限元中使用钕铁硼材料的剩磁密度

rpm=10000;     % rpm为电机转速
f=rpm*p/60;    % f为反电动势频率
T=1/f;         % T为转子转一周的周期
wm=2*pi*f;     % wm为机械角频率
we=wm*p;       % we为电角频率

t_accuracy=1000;
t_region=linspace(0,T,t_accuracy);           % 转子转动的时间范围，为转子转一圈的时间，即一个周期

% 在这组参数下，最大谐波次数取值。减小Rr,Rm,Rs半径间的差别，影响不大
N=100;         % 气隙及磁极谐波次
M=100;         % 齿槽谐波次数取值数取值。在这组参数下，最大谐波次数取值。减小Rr,Rm,Rs半径间的差别，影响不大
ii=100;
K=(2*ii-1)*p;

r=0.1094;          % r是所取的气隙半径，计算的是半径为r处的径向磁密
accuracy=1000;
alpha=linspace(0,pi,accuracy);
B2r_t=zeros(accuracy,t_accuracy);
B2a_t=zeros(accuracy,t_accuracy);

%% 循环开始
for t1 = 1:t_accuracy
    t=t_region(1,t1);                       % 转子转动的时刻
    %% 选择充磁方式 1.径向充磁 2.平行充磁

    Mrck=zeros(K,1);
    Mrsk=zeros(K,1);
    Mack=zeros(K,1);
    Mask=zeros(K,1);
    magnetization_method=1;

    switch(magnetization_method)
        case 1
            % radial magnetization 径向磁化
            for i = 1:ii
                n=2*i-1;
                k=n*p;
                Mrk=(4*p*Br/(k*pi*u0))*sin(k*pi*ap/(2*p));
                Mak=0;

                Mrck(k,1)=Mrk*cos(k*wm*t+k*a0);
                Mrsk(k,1)=Mrk*sin(k*wm*t+k*a0);
                Mack(k,1)=-Mak*sin(k*wm*t+k*a0);
                Mask(k,1)=Mak*cos(k*wm*t+k*a0);
            end

        case 2
            % parallel magnetization 平行磁化 k=1时A2k的分母为0要注意
            for i = 1:ii
                n=2*i-1;
                k=n*p;

                A1k=sin((k+1)*ap*(pi/(2*p)))/((k+1)*ap*(pi/(2*p)));
                if k == 1
                    A2k=1;
                else
                    A2k=sin((k-1)*ap*(pi/(2*p)))/((k-1)*ap*(pi/(2*p)));
                end

                Mrk=(Br/u0)*ap*(A1k+A2k);
                Mak=(Br/u0)*ap*(A1k-A2k);

                Mrck(k,1)=Mrk*cos(k*wm*t+k*a0);
                Mrsk(k,1)=Mrk*sin(k*wm*t+k*a0);
                Mack(k,1)=-Mak*sin(k*wm*t+k*a0);
                Mask(k,1)=Mak*cos(k*wm*t+k*a0);
            end
    end

    %% permanent magnet and air-gap 永磁体与气隙之间
    Ik=eye(K,K);
    K_matrix=eye(K,K);
    K_special=eye(K,K);
    for i = 1:K
        if i == 1
            K_special(i,i)=sqrt(2);
        else
            K_matrix(i,i)=i;
            K_special(i,i)=i;
        end
    end

    g1=Rr/Rm;
    g2=Rm/Rs;
    G1=eye(K,K);
    G2=eye(K,K);
    for i = 1:K
        G1(i,i)=g1^i;
        G2(i,i)=g2^i;
    end
    K11=Ik+G1^2;
    K22=Ik+G1^2;
    K13=-G2;
    K25=-G2;
    K14=-Ik;
    K26=-Ik;

    K31=Ik-G1^2;
    K42=Ik-G1^2;
    K33=-ur*G2;
    K45=-ur*G2;
    K34=ur*Ik;
    K46=ur*Ik;

    %% 定义Y向量

    Y1=-u0*(K_special^2-Ik)^(-1)*((Rr*K*G1+Rm*Ik)*Mack-(Rr*G1+Rm*K_matrix)*Mrsk);
    Y2=-u0*(K_special^2-Ik)^(-1)*((Rr*K*G1+Rm*Ik)*Mask+(Rr*G1+Rm*K_matrix)*Mrck);

    Y3=-u0*(K_special^2-Ik)^(-1)*(K_matrix*(Rm*Ik-Rr*G1)*Mack-(Rm*Ik-Rr*G1)*Mrsk);
    Y4=-u0*(K_special^2-Ik)^(-1)*(K_matrix*(Rm*Ik-Rr*G1)*Mask+(Rm*Ik-Rr*G1)*Mrck);

    %% slot opening and slot 槽开口与槽之间
    Fm=eye(M,M);
    En=eye(N,N);
    G3=eye(N,N);
    G4=eye(M,M);
    IN=eye(N,N);
    gamma=eye(M,N);

    for i = 1:M
        Fm(i,i)=i*pi/boa;
        G4(i,i)=(Rs/Rt)^(i*pi/boa);
    end
    for i = 1:N
        En(i,i)=i*pi/bsa;
        G3(i,i)=(Rt/Rsb)^(i*pi/bsa);
    end
    for i = 1:M
        for o = 1:N
            en=o*pi/bsa;
            fm=i*pi/boa;
            gamma(i,o)=-(2/bsa)*(en/(fm^2-en^2))*(cos(i*pi)*sin(en*(bsa+boa)/2)-sin(en*(bsa-boa)/2));
        end
    end
    K97=gamma'*Fm;
    K98=-gamma'*Fm*G4;
    K99=-En*(G3^2-IN);
    K97i=eye(N*Ns,M*Ns);
    K98i=eye(N*Ns,M*Ns);
    K99i=eye(N*Ns,N*Ns);
    Ftm=eye(M*Ns,M*Ns);%%
    G4t=eye(M*Ns,M*Ns);%%
    for i = 1:N
        for o = 1:M
            for n = 0:Ns-1
                K97i(i+N*n,o+M*n)=K97(i,o);
                K98i(i+N*n,o+M*n)=K98(i,o);
            end
        end
    end
    for i = 1:N
        for o = 1:N
            for n = 0:Ns-1
                K99i(i+N*n,o+N*n)=K99(i,o);
                Ftm(i+N*n,o+N*n)=Fm(i,o);
                G4t(i+N*n,o+N*n)=G4(i,o);
            end
        end
    end

    zeta=(bsa/boa)*gamma;
    IM=eye(M,M);
    K87i=eye(M*Ns,M*Ns);
    K88i=eye(M*Ns,M*Ns);
    K89i=eye(N*Ns,N*Ns);
    for i = 1:M
        for o = 1:M
            for n = 0:Ns-1
                K88i(i+M*n,o+M*n)=G4(i,o);
            end
        end
    end
    K89=-zeta*(G3^2+IN);
    for i = 1:N
        for o = 1:N
            for n = 0:Ns-1
                K89i(i+N*n,o+N*n)=K89(i,o);
            end
        end
    end

    %% slot opening and air gap 槽开口与气隙之间

    K53=-K_matrix;
    K65=-K_matrix;
    K54=G2*K_matrix;
    K66=G2*K_matrix;

    yita=zeros(K,M*Ns);
    kesei=zeros(K,M*Ns);

    for i = 1:M
        for o = 1:K
            for n = 0:Ns-1
                fm=i*pi/boa;
                ai=theta_s/2+theta_s*n;
                yita(o,i+M*n)=-(1/pi)*(o/(fm^2-o^2))*(cos(i*pi)*sin(o*ai+o*boa/2)-sin(o*ai-o*boa/2));
                kesei(o,i+M*n)=(1/pi)*(o/(fm^2-o^2))*(cos(i*pi)*cos(o*ai+o*boa/2)-cos(o*ai-o*boa/2));
            end
        end
    end

    K57=yita*Ftm*G4t;
    K58=-yita*Ftm;
    K67=kesei*Ftm*G4t;
    K68=-kesei*Ftm;

    sigma=(2*pi/boa)*yita';
    tao=(2*pi/boa)*kesei';
    K73=sigma;
    K74=sigma*G2;
    K75=tao;
    K76=tao*G2;
    K77=-G4t;
    K78=-eye(M*Ns,M*Ns);

    %% 拼凑系数矩阵

    Z1=zeros(K,K);
    Z2=zeros(K,M*Ns);
    Z3=zeros(M*Ns,K);
    Z4=zeros(M*Ns,M*Ns);

    ratio=[K11,Z1,K13,K14,Z1,Z1,Z2,Z2,Z2;
        Z1,K22,Z1,Z1,K25,K26,Z2,Z2,Z2;
        K31,Z1,K33,K34,Z1,Z1,Z2,Z2,Z2;
        Z1,K42,Z1,Z1,K45,K46,Z2,Z2,Z2;
        Z1,Z1,K53,K54,Z1,Z1,K57,K58,Z2;
        Z1,Z1,Z1,Z1,K65,K66,K67,K68,Z2;
        Z3,Z3,K73,K74,K75,K76,K77,K78,Z4;
        Z3,Z3,Z3,Z3,Z3,Z3,K87i,K88i,K89i;
        Z3,Z3,Z3,Z3,Z3,Z3,K97i,K98i,K99i];

    Y=[Y1;Y2;Y3;Y4;zeros(K*2+M*Ns*3,1)];

    %% 求出系数向量

    R=(ratio^(-1))*Y;

    A1=R(1:K,1);
    C1=R(K+1:2*K,1);
    A2=R(2*K+1:3*K,1);
    B2=R(3*K+1:4*K,1);
    C2=R(4*K+1:5*K,1);
    D2=R(5*K+1:6*K,1);
    C4t=R(6*K+1:6*K+M*Ns,1);
    D4t=R(6*K+M*Ns+1:6*K+M*Ns*2,1);
    D3t=R(6*K+M*Ns*2+1:6*K+M*Ns*3,1);

    %% 求气隙磁密分布

    for i = 1:accuracy
        level1=0;
        level2=0;
        level3=0;
        level4=0;
        for j = 1:ii
            n=2*j-1;
            o=p*n;
            level1=level1+o*((A2(o)/Rs)*(r/Rs)^(o-1)+(B2(o)/Rm)*(r/Rm)^(-o-1))*sin(p*o*alpha(i));
            level2=level2+o*((C2(o)/Rs)*(r/Rs)^(o-1)+(D2(o)/Rm)*(r/Rm)^(-o-1))*cos(p*o*alpha(i));
            level3=level3+o*((A2(o)/Rs)*(r/Rs)^(o-1)-(B2(o)/Rm)*(r/Rm)^(-o-1))*cos(p*o*alpha(i));
            level4=level4+o*((C2(o)/Rs)*(r/Rs)^(o-1)-(D2(o)/Rm)*(r/Rm)^(-o-1))*sin(p*o*alpha(i));
        end
        B2r_t(i,t1)=-level1+level2;
        B2a_t(i,t1)=-level3-level4;
    end
    disp("现在是第")
    disp(t1)
    disp("次循环")
end

%% 绘制气隙磁密波形

figure
plot(180*p*alpha/pi,B2r_t(:,1),'-k','DisplayName',"径向气隙磁密")
hold on
plot(180*p*alpha/pi,B2a_t(:,1),'-r','DisplayName',"周向气隙磁密")
hold on

xlim([0 360])
title("气隙磁密波形")
legend
%% 绘制气隙磁密波

figure
subplot(2,1,1)
waterfall(p*t_region,180*p*alpha./pi,B2r_t);
xlabel("时间 / s")
ylabel("转子位置角 / Electrical Degree")
zlabel("径向气隙磁密 / Tesla")
ylim([0,360]);
xlim([0,0.003]);
title("径向气隙磁密波")
colormap("jet")
colorbar
hold on

subplot(2,1,2)
waterfall(p*t_region,180*p*alpha./pi,B2a_t);
xlabel("时间 / s")
ylabel("转子位置角 / Electrical Degree")
zlabel("周向气隙磁密 / Tesla")
ylim([0,360]);
xlim([0,0.003]);
title("周向气隙磁密波")
colormap("jet")
colorbar
hold on

%% 绘制气隙磁密在时间维度上的分布

figure

subplot(2,1,1)
plot(t_region,B2r_t(1,:),'-k',"DisplayName","周向气隙磁密")
hold on

xlabel("时间 / s")
ylabel("气隙磁密")
title("气隙磁密在时间维度上的分布")
xlim([0,T])
legend

subplot(2,1,2)
plot(t_region,B2a_t(1,:),'-r',"DisplayName","径向气隙磁密")
hold on

xlabel("时间 / s")
ylabel("气隙磁密")
title("气隙磁密在时间维度上的分布")
xlim([0,T])
legend

%% 绘制反电动势波形

wm=2*pi*f;     % wm为机械角频率
we=wm*p;       % we为电角频率
v=wm*r;        % v为绕组的相对线速度

figure

induced_voltage=B2r_t.*L*v;
plot(t_region,16*induced_voltage(1,:),'-k',"DisplayName","线反电动势")

xlabel("时间 / s")
ylabel("反电动势 / V")
title("线反电动势波形")
xlim([0,T])
legend
%% 绘制电磁力分布图

u0=4*pi*1e-7;  % 真空磁导率
fr=(B2r_t.^2-B2a_t.^2)/(2*u0);
fa=(B2r_t.*B2a_t)/u0;

figure
mesh(t_region,180*p*alpha./pi,fr);
xl=xlabel('时间（s）');
yl=ylabel('空间机械角度（°）');
zlabel('径向电磁力密度（N/m^2）');
set(xl,'Rotation',15);
set(yl,'Rotation',-25);
set(gca,'FontSize',14);
title("径向电磁力波")
colormap("jet")
colorbar

figure
mesh(t_region,180*p*alpha./pi,fa);
xl=xlabel('时间（s）');
yl=ylabel('空间机械角度（°）');
zlabel('周向电磁力密度（N/m^2）');
set(xl,'Rotation',15);
set(yl,'Rotation',-25);
set(gca,'FontSize',14);
title("切向电磁力波")
colormap("jet")
colorbar


%% 对径向/周向电磁力密度波形进行空间维度上的一维FFT分析
figure

subplot(2,2,1)
plot(180*p*alpha/pi,fr(:,1),'-k','DisplayName',"径向电磁力");
hold on
legend
set(gca,'FontSize',14)
xlim([0 360])
xlabel("电角度 / deg")
ylabel("电磁力密度 / N/m^2")
title("径向电磁力密度")

subplot(2,2,2)
plot(180*p*alpha/pi,fa(:,1),'-r','DisplayName',"周向电磁力")
hold on
legend
set(gca,'FontSize',14)
xlim([0 360])
xlabel("电角度 / deg")
ylabel("电磁力密度 / N/m^2")
title("周向电磁力密度")

% 对径向电磁力密度进行FFT分析

subplot(2,2,3)
FFT_fr=fft(fr(:,1));
bar(0:1:31,abs(FFT_fr(1:2:64,1))/2002,'FaceColor','r','BarWidth',0.2)
set(gca,'XTick',0:1:32);
xlabel("空间阶数")
ylabel("谐波幅值 / N/m^2")
set(gca,'FontSize',14)
title("径向电磁力密度谐波分布")

% 对周向电磁力密度进行FFT密度

subplot(2,2,4)
FFT_fa=fft(fa(:,1));
bar(0:1:31,abs(FFT_fa(1:2:64,1))/2002,'FaceColor','r','BarWidth',0.2)
set(gca,'XTick',0:1:32);
xlabel("空间阶数")
ylabel("谐波幅值 / N/m^2")
set(gca,'FontSize',14)
title("周向电磁力密度谐波分布")

%% 对径向/周向电磁力密度波形进行时间维度上的一维FFT分析

% 对径向电磁力密度进行FFT分析
figure

subplot(2,2,1)
plot(t_region,fr(1,:),'-k','DisplayName',"径向电磁力");
hold on
legend
set(gca,'FontSize',14)
xlabel("时间 / s")
ylabel("电磁力密度 / N/m^2")
title("径向电磁力密度")

subplot(2,2,2)
plot(t_region,fa(1,:),'-r','DisplayName',"周向电磁力")
hold on
legend
set(gca,'FontSize',14)
xlabel("时间 / s")
ylabel("电磁力密度 / N/m^2")
title("周向电磁力密度")


subplot(2,2,3)
FFT_fr=fft(fr(1,:));
bar(0:1:31,abs(FFT_fr(1,1:2:64))/2002,'FaceColor','r','BarWidth',0.2)
set(gca,'XTick',0:1:32);
xlabel("时间次数")
ylabel("谐波幅值 / N/m^2")
set(gca,'FontSize',14)
title("径向电磁力密度谐波分布")

% 对周向电磁力密度进行FFT密度

subplot(2,2,4)
FFT_fa=fft(fa(1,:));
bar(0:1:31,abs(FFT_fa(1,1:2:64))/2002,'FaceColor','r','BarWidth',0.2)
set(gca,'XTick',0:1:32);
xlabel("时间次数")
ylabel("谐波幅值 / N/m^2")
set(gca,'FontSize',14)
title("周向电磁力密度谐波分布")

%%


FFT_fr=zeros(accuracy,t_accuracy);
FFT_fa=zeros(accuracy,t_accuracy);

for i = 1:t_accuracy
    FFT_fr(:,i)=fft(fr(:,i));
    FFT_fa(:,i)=fft(fa(:,i));
end

FFT_fr=abs(FFT_fr);
FFT_fa=abs(FFT_fa);

figure

bar3(0:1:31,FFT_fr(1:2:64,1:2:64))

set(gca,'XTick',0:1:32);
set(gca,'YTick',0:1:32);
xlabel("空间阶数")
ylabel("时间次数")
zlabel("谐波幅值 / N/m^2")
set(gca,'FontSize',14)
title("周向电磁力密度谐波分布")
%% 电磁力 2D FFT

Y = fft2(fr);
N_space = size(fr,1);
N_time = size(fr,2);
Ayy = abs(Y)/(N_space*N_time/2);
Ayy(1) = Ayy(1)/2;

figure
y=linspace(-64,64,129);
bar3(y,Ayy(1:129,1:65));

set(gca,'XTick',1:2:64);
set(gca,'XTickLabel',0:2:64,'FontSize',7);
set(gca,'YTick',-64:4:64)
set(gca,'YTickLabel',-64:4:64,'FontSize',7)
xlabel("时间次数",'FontSize',14)
ylabel("空间阶数",'FontSize',14)
zlabel("径向电磁力密度 N/m^2",'FontSize',14)
title("径向电磁力密度二维FFT分解",'FontSize',20)
xlim([0 65]);
ylim([-65 65]);
colormap([spring;summer;autumn;winter])
colorbar

%%
X = fft2(fa);
N_space = size(fa,1);
N_time = size(fa,2);
Ayy = abs(X)/(N_space*N_time/2);
Ayy(1) = Ayy(1)/2;

figure
y=linspace(-64,64,129);
bar3(y,Ayy(1:129,1:65));

set(gca,'XTick',1:2:64);
set(gca,'XTickLabel',0:2:64,'FontSize',7);
set(gca,'YTick',-64:4:64)
set(gca,'YTickLabel',-64:4:64,'FontSize',7)
xlabel("时间次数",'FontSize',14)
ylabel("空间阶数",'FontSize',14)
zlabel("切向电磁力密度 N/m^2",'FontSize',14)
title("切向电磁力密度二维FFT分解",'FontSize',20)
xlim([0 65]);
ylim([-65 65]);
colormap([spring;summer;autumn;winter])
colorbar


%% 计算电磁转矩

L=0.46;

dx=2*pi/accuracy;

Tcog=zeros(t_accuracy,1);
for t = 1:t_accuracy
    for i = 1:accuracy-1
        Tcog(t,1)=Tcog(t,1)+(r^2)*L*(fa(i,t)+fa(i+1,t))*dx/2;
    end
end

figure
plot(360*t_region/T,Tcog,'-k');

xlim([0 7.5])