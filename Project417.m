clc;clear all;close all;
%%%%
R_n=0.02; %k.ohm.cm^2 nodal resistance area 
R_m=100;%k.ohm.cm^2 myelin resistance area 
R_a=0.14;%k.ohm.cm^2 axoplasm resistivity 
C_n=5.0 ;%uF/cm^2 nodal capacitance per unit area 
C_m=5*10^-3;%uF/cm^-2 myelin capacitance per unit area 
%%%%myelin timw constant and space constant
lambda_m=0.433; %cm
t_m=500; %us
%%%%nodl time constant and space constant
lambda_n=0.0061; %cm
t_n=100; %us
%%%
delta=1e-4;%cm width of node of Ranvier 
d=0.7*1.5*1e-3;%cm inner diameter ofaxonal membrane 
%%%
alpha_d_start = 0.1;              
alpha_d_end = 0.9;              
alpha_d_delta = 1/20;              
alpha_d_grid = (alpha_d_end-alpha_d_start)/alpha_d_delta; % number of steps
alpha_d = linspace(alpha_d_start,alpha_d_end,alpha_d_grid);
alpha_normal=0.7;%0.77,0.6 place
%%%
D=d./alpha_d;%outer diameter of myelin sheath 
L=100.*D;
c_n=C_n*pi*d;%uF/cm nodal capacitance per unit length 
c_m=C_m*pi*d;%uF/cm myelin capacitance per unit length 
r_n=R_n/(pi*d);%k.ohm.cm nodal resistance length 
r_m=R_m/(pi*d);%k.ohm.cm myelin resistance length 
r_a=R_a/(pi*d*d);%k.ohm/cm axoplasm resistance 
%%%intergated lambda and t
lambda_T=sqrt((1-delta./L).*(r_a/r_m)+delta./L.*(r_a/r_n));
t_T=lambda_T.^2.*((1-delta./L)*r_a*c_m+r_a*c_n*delta./L);
%%%%%%%%
%%%set the time interval
t_start = 0;               % ms
t_end = 60;                 % ms
t_on= 3;
t_off=4;
t_delta = 1/20;              % ms,use approate interval 
t_grid = 1+(t_end-t_start)/t_delta; % number of steps
time = linspace(t_start,t_end,t_grid);% ms
%%% plot the stimulus pulse, generate the pulse signal
[n,m]=size(time);
pulse_signal=zeros(1,m);
for i=1:m
pulse_signal(1,i)=pulseside(time(1,i)-t_on)*pulseside(t_off-time(1,i));
end
plot(time,pulse_signal)

%%%%v_m, use t=1

%%%set the space interval
x_start = -50;               % mm
x_end = 50;                 % mm
x_delta = 1/20;              % mm,use approate interval 
x_grid = 1+(x_end-x_start)/x_delta; % number of steps
space = linspace(x_start,x_end,x_grid);% mm
parameter=lambda_T/(2*sqrt(pi*t_T));
for j=1:t_grid
    dt=(j-1)*t_delta;
    i=16;
%     for i=1:16
%         V_m{i*j}=@(x)exp(-dt./t_T(i)-t_T(i).*(x./(2.*lambda_T(i))).^2./dt);%t=1
%         q(i*j)=integral(V_m{i*j},-Inf,Inf);
%         V_m2{i*j}=@(t)q(i*j)./sqrt(t);
%         p(i*j)=integral(V_m2{i*j},0,dt);
%     end
        V_m{j}=@(x)exp(-dt./t_T(i)-t_T(i).*(x./(2.*lambda_T(i))).^2./dt);%t=1
        q(j)=integral(V_m{j},-Inf,Inf);
        V_m2{j}=@(t)q(j)./sqrt(t);
        p(j)=integral(V_m2{j},0,dt);
end
for j=1:t_grid
    dt=(j-1)*t_delta;
    i=14;
%     for i=1:16
%         V_m{i*j}=@(x)exp(-dt./t_T(i)-t_T(i).*(x./(2.*lambda_T(i))).^2./dt);%t=1
%         q(i*j)=integral(V_m{i*j},-Inf,Inf);
%         V_m2{i*j}=@(t)q(i*j)./sqrt(t);
%         p(i*j)=integral(V_m2{i*j},0,dt);
%     end
        V_m{j}=@(x)exp(-dt./t_T(i)-t_T(i).*(x./(2.*lambda_T(i))).^2./dt);%t=1
        q(j)=integral(V_m{j},-Inf,Inf);
        V_m2{j}=@(t)q(j)./sqrt(t);
        P(j)=integral(V_m2{j},0,dt);
end
for j=1:t_grid
    dt=(j-1)*t_delta;
    i=12;
%     for i=1:16
%         V_m{i*j}=@(x)exp(-dt./t_T(i)-t_T(i).*(x./(2.*lambda_T(i))).^2./dt);%t=1
%         q(i*j)=integral(V_m{i*j},-Inf,Inf);
%         V_m2{i*j}=@(t)q(i*j)./sqrt(t);
%         p(i*j)=integral(V_m2{i*j},0,dt);
%     end
        V_m{j}=@(x)exp(-dt./t_T(i)-t_T(i).*(x./(2.*lambda_T(i))).^2./dt);%t=1
        q(j)=integral(V_m{j},-Inf,Inf);
        V_m2{j}=@(t)q(j)./sqrt(t);
        P12(j)=integral(V_m2{j},0,dt);
end
for j=1:t_grid
    dt=(j-1)*t_delta;
    i=1;
%     for i=1:16
%         V_m{i*j}=@(x)exp(-dt./t_T(i)-t_T(i).*(x./(2.*lambda_T(i))).^2./dt);%t=1
%         q(i*j)=integral(V_m{i*j},-Inf,Inf);
%         V_m2{i*j}=@(t)q(i*j)./sqrt(t);
%         p(i*j)=integral(V_m2{i*j},0,dt);
%     end
        V_m{j}=@(x)exp(-dt./t_T(i)-t_T(i).*(x./(2.*lambda_T(i))).^2./dt);%t=1
        q(j)=integral(V_m{j},-Inf,Inf);
        V_m2{j}=@(t)q(j)./sqrt(t);
        P1(j)=integral(V_m2{j},0,dt);
end
plot(time,p)
hold on
plot(time,P)
plot(time,P12)
plot(time,P1)
leghndl = legend('0.8467','0.7933','0.6867','0.1')
%%%%%functions
        function pulse=pulseside(x)
            if x >= 0
                pulse = 1;
            else
                pulse = 0;
            end
        end