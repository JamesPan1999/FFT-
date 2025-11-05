% Identification using the least squares method, based on Simulink data
% dataout.data is obtained from Simulink 

% Display the data from the Simulink
% 如果simulink的求解器采用向前欧拉ode1，且采用Ts作为采样周期，那么将得到很好的拟合效果，但是得不到准确的传递函数。这是因为用欧拉进行模型离散化太粗糙了。
% thata = [a1 b1]';
% phik = [-y(k-1) u(k-1)]';
% y(k) = -y(k-1)*a1 + u(k-1)*b1;   % 由于simulink已经帮我们计算了整个动态过程的时间序列，这个就不需要了

Ts=0.1;                         % The sampling period
Time=(0:1:size(dataout.data(:,1))-1)*Ts;
plot(Time,dataout.data);

% Initialize the parameters
na=1;
nb=1;
n=max(na,nb);              		% Depth of historical data
N=size(dataout.data); 			% The ‘dataout” stores the input-output data
N=N(1);                         % Total number of historical data
phph=zeros(na+nb);
phy=zeros(na+nb,1);
for t=n+1:N     				% The loop part
    tmp=[];
    for j=1:na tmp=[tmp,-dataout.data(t-j,1)]; end
    for j=1:nb tmp=[tmp,dataout.data(t-j,2)]; end
    phph=phph+tmp'*tmp;                % Phi'*Phi    直接吧Phi矩阵乘开了，变成phi向量相乘
    phy=phy+tmp'*dataout.data(t,1);    % Phi'*Y      按分块矩阵乘。也可以全部拼接完后统一乘  不过这种方法有递归最小二乘的味道
end

% Estimateθof the discrete model
theta=inv(phph)*phy			 	

% Coverte the discrete model to a continuous model
sysd=tf(theta(2), [1 theta(1)],Ts);
d2c(sysd)

%%
% Time
% dataout.data(:,1);
% dataout.data(:,2);

% 运行前clear
Ts=0.1;                         % The sampling period
Time=(0:1:size(dataout.data(:,1))-1)*Ts;
d = 1;

Phi = [];  Y = [];
for i = 1:length(Time)-d 
    phik = [-dataout.data(i,1) dataout.data(i,2)];
    yk = dataout.data(i+d,1);

    Phi = [Phi;phik];
    Y = [Y;yk];
end

thetaE = inv(Phi'*Phi)*Phi'*Y