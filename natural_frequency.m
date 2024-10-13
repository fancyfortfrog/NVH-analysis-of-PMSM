% clear
% clc
%% Institude: HuaZhong University of Science and Technology
% 机构： 华中科技大学电气与电子工程学院
%% Written by Xu Shouyu
% 作者： 徐首彧
% 指导老师： 叶才勇
%% 1MW大电机参数

% Rout=0.184;                % 定子铁芯外半径
% Rsb=0.16;                  % 定子槽底部半径
% Rt=0.115;                  % 定子槽顶部半径
% Rs=0.113;                  % 定子内半径
% Rm=0.105;                  % 磁极半径
% Rr=0.092;                  % 转子半径
% Df=0.4;                    % 机座直径
% boa=0.02655;               % boa是槽开口角
% bsa=0.0607;                % bsa是槽宽角
%
% % 定子铁芯材料为B20AT1500硅钢片
% vc=0.27;                   % 定子铁芯泊松比
% Rho_c=7650;                % 定子铁芯密度
% Ec=211e9;                  % 定子铁芯弹性模量
%
% % 绕组材料为铜
% vw=0.35;                   % 绕组铜泊松比
% Rho_w=8890;                % 绕组铜密度
% Ew=115e9;                  % 绕组铜弹性模量

%% 《多相电机噪声书上的参数》
Rout=0.233/2;                 % 定子外半径
hc=8.42e-3;                   % 定子铁轭厚度
Rsb=Rout-hc;                  % 定子槽底半径

s1=36;                        % 定子槽数
ht=0.028;                     % 定子齿长度
ct=5.87e-3;                   % 定子齿宽
Rs=0.16/2;                    % 定子内半径
Lc=0.1975;                    % 定子有效长度
hov=0.048;                    % 绕组端部接线长度

Df=0.246;                     % 机座直径
Lf=0.359;                     % 机座长度
vf=0.33;                      % 机座泊松比
Ef=71e9;                      % 机座弹性模量
Rho_f=2700;                   % 机座密度

vc=0.3;                       % 定子铁芯泊松比
Rho_c=7700;                   % 定子铁芯密度
Ec=200e9;                     % 定子铁心弹性模量

vw=0.35;                      % 绕组泊松比
Rho_w=8890;                   % 绕组密度
Ew=9.4e9;                     % 绕组弹性模量


m_region=linspace(0,8,9);  % 周向模态阶次m的取值范围，m取0,1,2,3,4,5,6,7,8


%% 计算定子铁心的固有频率与集中刚度

hc=Rout-Rsb;               % 定子铁轭厚度，hc=Rout-Rsb
Rc=(2*Rout-hc)/2;          % 定子铁轭平均半径Rc=(D1out-hc)/2
k2_stator=(hc^2)/(12*Rc^2);           % 参数k的计算

ohm_stator=zeros(length(m_region),2); % 分配定子铁心参数ohm_m的存储空间
f_stator=zeros(length(m_region),2);   % 定子铁芯的固有频率
K_stator=zeros(length(m_region),1);   % 定子铁芯的集中刚度

for i = 1:9
    m=i-1;
    if m == 0
        ohm_stator(i,1)=1;
        ohm_stator(i,2)=1;
    else
        ohm_stator(i,1)=0.5*sqrt((1+m^2+k2_stator*m^4)-sqrt((1+m^2+k2_stator*m^4)^2-4*k2_stator*m^6));
        ohm_stator(i,2)=0.5*sqrt((1+m^2+k2_stator*m^4)+sqrt((1+m^2+k2_stator*m^4)^2-4*k2_stator*m^6));
    end
    f_stator(i,1)=ohm_stator(i,1)*sqrt(Ec/(Rho_c/(1-vc^2)))/(pi*2*Rc);  
    f_stator(i,2)=ohm_stator(i,2)*sqrt(Ec/(Rho_c/(1-vc^2)))/(pi*2*Rc);
    K_stator(i,1)=(4*ohm_stator(i,1)^2)/(2*Rc)*(pi*Lc*hc*Ec)/(1-vc^2);   % 定子铁芯的集中刚度采用(5-25)式计算
end

%% 计算绕组的固有频率与集中刚度

Rw=(2*Rs+ht)/2;            % 定子绕组的平均半径
k2_wingdings=(ht^2)/(12*Rw^2);       % 参数k的计算

ohm_windings=zeros(length(m_region),2); % 分配绕组参数ohm_m的存储空间
f_windings=zeros(length(m_region),2);   % 绕组的固有频率
K_windings=zeros(length(m_region),1);   % 绕组的集中刚度

for i = 1:9
    m=i-1;
    if m == 0
        ohm_windings(i,1)=1;
        ohm_windings(i,2)=1;
    else
        ohm_windings(i,1)=0.5*sqrt((1+m^2+k2_wingdings*m^4)-sqrt((1+m^2+k2_wingdings*m^4)^2-4*k2_wingdings*m^6));
        ohm_windings(i,2)=0.5*sqrt((1+m^2+k2_wingdings*m^4)+sqrt((1+m^2+k2_wingdings*m^4)^2-4*k2_wingdings*m^6));
    end
    f_windings(i,1)=ohm_windings(i,1)*sqrt(Ew/(Rho_w/(1-vw^2)))/(pi*2*Rw);
    f_windings(i,2)=ohm_windings(i,2)*sqrt(Ew/(Rho_w/(1-vw^2)))/(pi*2*Rw);
    K_windings(i,1)=(ohm_windings(i,1)^2)*Ew*vw/(Rw^2)/(1-vw^2);         % 绕组的集中刚度采用(5-39)式计算
end

%% 计算机座的固有频率与集中刚度

Rf=(Df+2*Rout)/4;          % 机座平均半径
hf=(Df-2*Rout)/2;          % 机座厚度
k2_saddle=(hf^2)/(12*Rf^2);          % 参数k的计算

lamda_e=zeros(3,1);        % 含端盖的圆柱形机座可视为两端固定壳体，等效波长为lamda_e
ohm_frame=zeros(3,9);         % 存储机座的最小多根Rho_mn
f_frame=zeros(3,9);        % 机座的固有频率
K_frame=zeros(3,9);        % 机座的集中刚度

% 计算机座的等效波长
for n = 1:3
    L0=Lf*(0.3/(0.3+n));
    lamda_e(n,1)=n*pi*Rf/(Lf-L0);
end

% 计算机座的最小多根Rho_mn
for i = 1:9
    m=i-1;
    for n = 1:3
        % 根据Donnell-Mushtari理论计算
        delta=m^2+lamda_e(n,1)^2;

        C2=1+0.5*(3-vf)*delta+k2_saddle*delta^2;
        C1=0.5*(1-vf)*((3+2*vf)*lamda_e(n,1)^2+m^2+delta^2+(3-vf)/(1-vf)*k2_saddle*delta^2);
        C0=0.5*(1-vf)*((1-vf^2)*lamda_e(n,1)^4+k2_saddle*delta^4);
        equation=[1,0,-C2,0,C1,0,-C0];

        roots_equation=roots(equation);
        for o = 1:length(roots_equation)
            % 去除虚数根，虚数根值变为0
            if isreal(roots_equation(o))
            else
                roots_equation=0;
            end
        end
        % 去除值为0和小于0的多根
        roots_equation(roots_equation==0)=[];
        roots_equation(roots_equation < 0)=[];
        % 找到最小多根
        [roots_min,index]=min(roots_equation);
        % 存储最小多根
        ohm_frame(n,i)=roots_min;
        f_frame(n,i)=ohm_frame(n,i)*sqrt(Ef/Rho_f/(1-vf^2))/(2*pi*Rf);      % 机座的固有频率采用(5-33)式计算
        K_frame(n,i)=2*(ohm_frame(n,i)^2)*pi*Lf*hf*Ef/Rf/(1-vf^2);          % 机座的集中刚度采用(5-34)式计算
    end
end


%% 计算定子系统的固有频率

V_stator=Lc*pi*(Rout^2-Rsb^2)+Lc*ct*ht*s1;           % 含齿铁芯的体积；
V_windings=(Lc+2*hov)*pi*(Rsb^2-Rs^2)-Lc*ct*ht*s1;   % 绕组的体积
V_frame=Lf*pi*((Df/2)^2-Rout^2);                     % 机座的体积

Mc=Rho_c*V_stator;                                   % 定子铁芯的质量
Mw=Rho_w*V_windings;                                 % 绕组的质量
Mf=Rho_f*V_frame;                                    % 机座的质量

Mc=18.85;
Mw=12.02;
Mf=6.24;

f_mn=zeros(3,9);                      % 带绕组的定子-机座系统的固有频率

for i = 1:9
    m=i-1;
    for n = 1:3
        f_mn(n,i)=sqrt((K_stator(i,1)+K_frame(n,i)+K_windings(i,1))/(Mc+Mf+Mw))/(2*pi);
    end
end

%%

figure
plotmatrix(f_stator())
