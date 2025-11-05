% Batch least squares parameter estimation (LS)

clear all;
a=[1 -1.5 0.7]'; b=[1 0.5]'; d=3;       % Plant parameters
na=length(a)-1; nb=length(b)-1;         % na and nb are the order of A and B
L=300;                                  % Data length
Uk=zeros(d+nb,1);                       % Initial values of the input: Uk(i) means u(k-i)
Yk=zeros(na,1);                         % Initial values of the output 
Reg1=1; Reg2=1; Reg3=1; Reg4=0; Sq=1; 	% Initial values of the shift registers and square wave 
xi=sqrt(0.1)*randn(L,1);                  % White noise sequence
Theta=[a(2:na+1);b];                    % True values of the plant parameter

% 计算逆M序列作为输入时，输出的情况。
for k=1:L                               % Phi(k,:) is the row vector, 
    Phi(k,:)=[-Yk;Uk(d:d+nb)]';         %     which is convenient to form phi matrix
    Y(k)=Phi(k,:)*Theta+xi(k);          % Output data
    Ms=xor(Reg3,Reg4);                  % Generate M-sequence
    IMs=xor(Ms,Sq);                     % Generate inverse M-sequence
    if IMs==0  u(k)=-1;  else u(k)=1;  end
    Sq=not(Sq);                         % Generate a square wave 
    
   % Update data
    Reg4=Reg3; Reg3=Reg2; Reg2=Reg1; Reg1=Ms;
    for i=d+nb:-1:2
        Uk(i)=Uk(i-1);
    end
    Uk(1)=u(k);
    for i=na:-1:2
        Yk(i)=Yk(i-1);
    end
    Yk(1)=Y(k);
end

ThetaE=inv(Phi'* Phi)*Phi'* Y'      % Calculate the parameter estimate thetae
error = Theta - ThetaE
plot(u)
hold on
plot(Y')
%% 如果用单位阶跃  后面两个参数会辨识的很不准确, 这是应为后两个参数的随机性不够，信噪比太小. 
%  phi矩阵条件数太大。（线性映射的椭圆性太大。一种情况是多个列向量夹角很小；另一种情况是有一个列向量模长很小。）
% Batch least squares parameter estimation (LS)

clear all;

% 微分方程初值
a=[1 -1.5 0.7]'; b=[1 0.5]'; d=3;       % Plant parameters
na=length(a)-1; nb=length(b)-1;         % na and nb are the order of A and B
L=300;                                  % Data length
Uk=zeros(d+nb,1);                       % Initial values of the input: Uk(i) means u(k-i)
Yk=zeros(na,1);                         % Initial values of the output 

xi=sqrt(0.0001)*randn(L,1);                  % White noise sequence
Theta=[a(2:na+1);b];                    % True values of the plant parameter

% 构造系统阶跃输入u
u = zeros(L,1);
u(L/4:end) = 1;

for k = 1:L
    
    % one step
    Phi(k,:) = [-Yk;Uk(d:d+nb)]';    % u(k-1)  u(k-2)由于在差分方程中不存在（即系数为0，求逆将产生奇异问题。且最小二乘的原理本身构成方程组时，就没有加入这两项）
                                     % 但是在记录时，需要这两项进行迭代。
    Y(k) = Phi(k,:)*Theta+xi(k);
    
    % 更新
    for i = length(Yk):-1:2
        Yk(i) = Yk(i-1);
    end
    Yk(1) = Y(k);
    
    for i = length(Uk):-1:2
        Uk(i) = Uk(i-1);
    end
    Uk(1) = u(k);
end

ThetaE = inv(Phi'*Phi)*Phi'*Y'
error = Theta - ThetaE
plot(u)
hold on
plot(Y')
%% 使用sin函数  也不行，后两个参数的结果也非常糟糕。
%  输入信号的自相关性直接决定了频谱的情况。或者说两者相互决定。频谱分布越均匀，自相关性越低。频谱分布越不均匀（sin就会有尖峰），自相关性越高。
%  自相关性又决定了Phi矩阵的列向量的线性无关成度。

clear all;

% 微分方程初值
a=[1 -1.5 0.7]'; b=[1 0.5]'; d=3;       % Plant parameters
na=length(a)-1; nb=length(b)-1;         % na and nb are the order of A and B
L=300;                                  % Data length
Uk=zeros(d+nb,1);                       % Initial values of the input: Uk(i) means u(k-i)
Yk=zeros(na,1);                         % Initial values of the output 

xi=sqrt(0.1)*randn(L,1);                  % White noise sequence
Theta=[a(2:na+1);b];                    % True values of the plant parameter

% 构造系统正弦输入u  似乎 就是没有M序列得到的效果优秀，即使归一化或者使用SVD分解
ts = 0.01;
u = sin(0:ts:ts*L);


for k = 1:L
    
    % one step
    Phi(k,:) = [-Yk;Uk(d:d+nb)]';    % u(k-1)  u(k-2)由于在差分方程中不存在（即系数为0，求逆将产生奇异问题。且最小二乘的原理本身构成方程组时，就没有加入这两项）
                                     % 但是在记录时，需要这两项进行迭代。
    Y(k) = Phi(k,:)*Theta+xi(k);
    
    % 更新
    for i = length(Yk):-1:2
        Yk(i) = Yk(i-1);
    end
    Yk(1) = Y(k);
    
    for i = length(Uk):-1:2
        Uk(i) = Uk(i-1);
    end
    Uk(1) = u(k);
end

% 归一化
for i = 1:size(Phi,2)
    zeta(i) = (max(Phi(:,i))-min(Phi(:,i)));
end
Zeta = diag(zeta);
cond(Phi)
PPhi = [Phi(:,1)./zeta(1),Phi(:,2)./zeta(2),Phi(:,3)./zeta(3),Phi(:,4)./zeta(4)];

[U, S, V] = svd(Phi); 
n = size(Phi,2);
m = size(Phi,1);
% 用奇异值分解或者归一化Phi矩阵可以获得一个相对好的结果，但输入数据自相关性太强时，用奇异值和归一化也无济于事
% 构造奇异值矩阵的伪逆Sigma_plus
% S是m×n对角矩阵，伪逆是n×m对角矩阵（非零奇异值取倒数，零奇异值保留0）
Sigma_plus = zeros(n, m);
for i = 1:min(m, n)
    if S(i, i) > 1e-10  % 忽略接近0的奇异值（避免数值误差）
        Sigma_plus(i, i) = 1 / S(i, i);
    end
end
% 计算SVD解：x_svd = V * Sigma_plus * U' * y
x_svd = V * Sigma_plus * U' * Y'

error = Theta - x_svd
% ThetaE = inv(Phi'*Phi)*Phi'*Y'
% error = Theta - ThetaE
% plot(u)
% hold on
% plot(Y')
% cond(Phi'*Phi)


