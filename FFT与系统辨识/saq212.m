% Identification using the least squares method, based on Simulink data
% dataout.data is obtained from Simulink 

% Display the data from the Simulink
Ts=0.1;                         % The sampling period
Time=[0:1:size(dataout.data(:,1))-1]*Ts;
plot(Time,dataout.data);

% Initialize the parameters
na=2;
nb=2;
n=max(na,nb);              		% Depth of historical data
N=size(dataout.data); 			% The ‘dataout” stores the input-output data
N=N(1);                         % Total number of historical data
phph=zeros(na+nb);
phy=zeros(na+nb,1);
for t=n+1:N     				% The loop part
    tmp=[];
    for j=1:na tmp=[tmp,-dataout.data(t-j,1)]; end
    for j=1:nb tmp=[tmp,dataout.data(t-j,2)]; end
    phph=phph+tmp'*tmp;    
    phy=phy+tmp'*dataout.data(t,1);
end

% Estimateθof the discrete model
theta=inv(phph)*phy			 	

% Coverte the discrete model to a continuous model
sysd=tf([theta(3),theta(4)],[1 theta(1),theta(2)],Ts);
d2c(sysd)

%% y(k) = -a1*y(k-1)-a2*y(k-2) + b1*u(k-1) + b2*u(k-2);    向前欧拉  是完全不一样的两种结果  一般向前欧拉效果会好，应为展示了系统的因果性，更接近于采用零阶保持器推状态空间的结果
%  y(k) = -a1*y(k-1)-a2*y(k-2)-a3*y(k-3) + b1*u(k)+b2*u(k-1)+b3*u(k-2);   b1=0 a3=0
%  y(k) = -a1*y(k-1)-a2*y(k-2)-a3*y(k-3) + b1*u(k)+b2*u(k-1)+b3*u(k-2)+b4*u(k-3); 
%  u(k-3)列的加入将会引发Phi矩阵的线性相关  条件数激增。why？ 

%  只用关心 要写到k减几，以及u比y滞后的阶数
%  对于有高阶微分项，必须按照中心离散  否则将产生巨大的数值误差（相当于离散后传递函数变阶）
%  y(k+1)-y(k)-(y(k)-y(k-1)) + 2(y(k)-y(k-1)) + 3*y(k) = 4*(u(k)-u(k-1)) + 5*u(k)
%  y(k+1)-y(k)-(y(k)-y(k-1)) + 2(y(k+1)-y(k)) + 3*y(k+1) = 4*(u(k+1)-u(k)) + 5*u(k+1) 不可取，误差过大了 相当于以y(k)为中心点的同时，却取了u(k+1)为中心点离散
%  y(k+1)-y(k)-(y(k)-y(k-1)) + 2(y(k+1)-y(k)) + 3*y(k) = 4*(u(k+1)-u(k)) + 5*u(k)
%  y(k)-y(k-1)-(y(k-1)-y(k-2)) + 2(y(k)-y(k-1)) + 3*y(k) = 4*(u(k)-u(k-1)) + 5*u(k)  不可取，误差过大了
%  为了避免这种问题，最优的方法是写成扩展状态空间后采用单位保持离散化状态空间（直接调用matlab）

%  拿 k = 3 举例  
Ts=0.1;                         % The sampling period
Time=[0:1:size(dataout.data(:,1))-1]*Ts;
plot(Time,dataout.data);
d = 2;                          % 用向前欧拉快速推导一下，会发现一阶系统差分方程左右时间步相差1  二阶系统时间步相差2

Phi = []; Y = [];
nn = 4;
for i = d+nn:length(Time)
    % phik = [-dataout.data(i-1,1) -dataout.data(i-2,1) dataout.data(i-1,2) dataout.data(i-2,2)];  % 向前欧拉
    % phik = [-dataout.data(i-1,1) -dataout.data(i-2,1) dataout.data(i,2) dataout.data(i-1,2)];  % 
    % phik = [-dataout.data(i-1,1) -dataout.data(i-2,1) dataout.data(i,2) dataout.data(i-1,2),dataout.data(i-2,2)];  % 向前欧拉
    % phik = [-dataout.data(i-1,1) -dataout.data(i-2,1) -dataout.data(i-3,1) dataout.data(i,2) dataout.data(i-1,2),dataout.data(i-2,2)  ];  % 向前欧拉
    % phik = [-dataout.data(i-1,1) -dataout.data(i-2,1) -dataout.data(i-3,1) -dataout.data(i-4) -dataout.data(i-5) -dataout.data(i-6) -dataout.data(i-7) dataout.data(i,2) dataout.data(i-1,2),dataout.data(i-2,2)  ];  % 向前欧拉
    % phik = [-dataout.data(i-1,1) -dataout.data(i-2,1) -dataout.data(i-3,1)  dataout.data(i,2) dataout.data(i-1,2),dataout.data(i-2,2) dataout.data(i-3,2) ]; 
    phik =[];
    for j = 1:nn
        phik = [phik -dataout.data(i-j)];
    end
    
    phik = [phik dataout.data(i,2) dataout.data(i-1,2),dataout.data(i-2,2)];
    yk = [dataout.data(i,1) ];

    Phi = [Phi;phik];
    Y = [Y;yk];
end

% phik的列越多，矩阵条件数越大，越不容易辨识，当超过一定列时，条件数激增
cond(Phi)
theta = inv(Phi'*Phi)*Phi'*Y
% sysd=tf([theta(3),theta(4)],[1 theta(1),theta(2)],Ts);
sysd=tf([theta(3),theta(4),theta(5)],[1 theta(1),theta(2)],Ts);
d2c(sysd)
%%  
% 程序功能：从传递函数生成状态空间模型，并进行离散化（修正版）

% 1. 定义连续传递函数 G(s) = num(s)/den(s)
num = [4 5];       % 分子：s + 3
den = [1 2 3];     % 分母：s^2 + 2s + 5
G_cont = tf(num, den);  % 创建连续传递函数对象

disp('连续传递函数 G(s)：');
disp(G_cont);


% 2. 将传递函数转换为状态空间模型
[A, B, C, D] = tf2ss(num, den);  % 能控标准形状态空间

disp('\n连续状态空间模型：');
disp(['A = ']); disp(A);
disp(['B = ']); disp(B);
disp(['C = ']); disp(C);
disp(['D = ']); disp(D);


% 3. 状态空间离散化（修正：先创建连续状态空间对象，再调用c2d）
Ts = 0.1;                  % 采样时间（秒）
sys_cont = ss(A, B, C, D); % 将A,B,C,D封装为连续状态空间模型对象
sys_disc = c2d(sys_cont, Ts, 'zoh');  % 离散化（ZOH方法）

% 从离散状态空间对象中提取矩阵
A_d = sys_disc.A;
B_d = sys_disc.B;
C_d = sys_disc.C;
D_d = sys_disc.D;

disp(['\n离散化状态空间（采样时间 Ts = ', num2str(Ts), ' s）：']);
disp(['A_d = ']); disp(A_d);
disp(['B_d = ']); disp(B_d);
disp(['C_d = ']); disp(C_d);
disp(['D_d = ']); disp(D_d);


% 验证：离散化后的传递函数
G_disc = tf(sys_disc);  % 从离散状态空间转换为离散传递函数
disp('\n离散化后的传递函数 G(z)：');
disp(G_disc.Numerator);
disp(G_disc.Denominator);


%%
%% 2. 向前欧拉离散化  
Ts = 0.1;
A_d_forward = eye(2) + Ts * A  % A_d = I + T*A
B_d_forward = Ts * B           % B_d = T*B
% C_d_forward = C
% D_d_forward = D
sys_cont_forward = ss(A_d_forward,B_d_forward,C_d_forward,D_d_forward,Ts);
G_disc_forward = tf(sys_disc);
G_disc_forward.Numerator
G_disc_forward.Denominator
%% 总结一个不会错的方法：
% 1. 先连续建模（最小二乘必须已知系统形态）
% 2. 写成状态空间形式
% 3. 离散化
% 4. 将离散状态空间tf转换成传递函数分子分母形式。看阶数
% 更简单的是直接看连续系统分子分母阶数之差，定时间滞后步数