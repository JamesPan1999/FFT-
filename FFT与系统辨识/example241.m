% Recursive generalised least squares parameter estimation (RELS)

clear all; close all;
% a=[1 -1.5 0.7]'; b=[1 0.5]'; c=[1 -1 0.2]'; d=3; 	% Plant parameters   只能辨识c(1) = 1 的情况
a=[1 1.6 0.7]'; b=[1.0 0.4]'; c=[0.9 1.2 0.3]'*1/0.9; d=1; 	% Plant parameters
na=length(a)-1; nb=length(b)-1; nc=length(c)-1; 	% na, nb and nc are the order of  A, B and C

L=10000;                         % Simulation length
t = 1:L;
M_seq13 = PRBS13(t,13,1);
M_seq8 = PRBS8(t,8,1);

uk=zeros(d+nb,1);               % Initial values of the input: uk(i) means u(k-i)
yk=zeros(na,1); 				% Initial values of the uutput 
xik=zeros(nc,1); 				% Initial values of the noise
xiek=zeros(nc,1);               % Initial value of noise estimation
u=M_seq13;                   % Input adopts white noise sequence
xi=M_seq8; 		% White noise sequence
theta=[a(2:na+1);b;c(2:nc+1)]; 	% Plant parameters
thetaE_1=zeros(na+nb+1+nc,1); 	% na+nb+1+nc is the number of 
                                %       identification parameters
P=10^10*eye(na+nb+1+nc);
for k=1:L  % k =7
    uRecoder(k) = uk(1);
    xiRecoder(k) = xik(1);

    phi=[-yk;uk(d:d+nb);xik];
    y(k)=phi'*theta+xi(k); 		% Collect output data
    phie=[-yk;uk(d:d+nb);xiek]; % Set up phie

   %Recursive augmented least squares method
    K=P*phie/(1+phie'*P*phie);
    thetaE(:,k)=thetaE_1+K*(y(k)-phie'* thetaE_1);
    P=(eye(na+nb+1+nc)-K*phie')*P;
    xie=y(k)-phie'* thetaE(:,k);  % Estimated value of white noise

    yRecoder(k) = y(k);
    %Update data
    thetaE_1=thetaE(:,k);
    for i=d+nb:-1:2  uk(i)=uk(i-1); end
    uk(1)=u(k);
    for i=na:-1:2  yk(i)=yk(i-1); end
    yk(1)=y(k);
    for i=nc:-1:2  xik(i)=xik(i-1);xiek(i)=xiek(i-1);end
    xik(1)=xi(k);
    xiek(1)=xie;
end

figure(1);
plot([1:L],thetaE(1:na,:));
xlabel('k'); 
legend('a_1','a_2'); 
axis([0 L -2 2]);

figure(2);
plot([1:L],thetaE(na+1:na+nb+1,:));
xlabel('k'); 
legend('b_0','b_1'); 
axis([0 L 0 1.5]);

figure(3);
plot([1:L],thetaE(na+nb+2:na+nb+nc+1,:));
xlabel('k'); 
legend('c_1','c_2'); 
axis([0 L -2 2]);

thetaE(:,L)



function Out = PRBS13(t,n,a)
% a: the altitude of M-sequence, 
% t为生成的时间一共多长
% n用来控制一个PRBS序列的周期长度  一个周期＝2^n-1
% Generate a PRBS signal
% 1. 自相关性尽可能为δ函数（低相关性）。  2. 尽可能均匀分布(接近白噪声）  3. 但又要可被重复实现 
% n=4;                                % the number of registers
% L=15+1;                             % the length of M-sequence to be required   越长越接近白噪声
L = length(t);

for i = 1:n  Reg(i) = 1; end        % Initialize n registers  位移寄存器全1
Out(1)=-a;

for i=2:L
        temp=Reg(1);
        % 注意反馈抽头 可用GPT帮忙写
        % 抽头为第4、5、6、8级，异或结果反馈到第1级  仅适用于n = 8的情况
        % Reg(1) = xor(xor(xor(Reg(4), Reg(5)), Reg(6)), Reg(8)); 
        %% n = 13
        Reg(1) = xor(xor(xor(Reg(1), Reg(3)), Reg(4)), Reg(13));
        %%  仅适用于n= 4  情况
        % Reg(1)=xor(Reg(n-1),Reg(n));  % Modulo 2 adder   第一个寄存器 
        j=2;
        while j<=n                  % Registers shift                
             temp1=Reg(j);
             Reg(j)=temp;
             j=j+1;
             temp=temp1;
        end        
        Out(i)=Reg(n);
        if (Reg(n)==1)  Out(i) =-a; end
        if (Reg(n)==0)  Out(i)=a; end
end   

% stairs(Out); axis([1 L -0.5 1.5]);  % Plot the PRBS signal

end


function Out = PRBS8(t,n,a)
% a: the altitude of M-sequence, 
% t为生成的时间一共多长
% n用来控制一个PRBS序列的周期长度  一个周期＝2^n-1
% Generate a PRBS signal
% 1. 自相关性尽可能为δ函数（低相关性）。  2. 尽可能均匀分布(接近白噪声）  3. 但又要可被重复实现 
% n=4;                                % the number of registers
% L=15+1;                             % the length of M-sequence to be required   越长越接近白噪声
L = length(t);

for i = 1:n  Reg(i) = 1; end        % Initialize n registers  位移寄存器全1
Out(1)=-a;

for i=2:L
        temp=Reg(1);
        % 注意反馈抽头 可用GPT帮忙写
        % 抽头为第4、5、6、8级，异或结果反馈到第1级  仅适用于n = 8的情况
        Reg(1) = xor(xor(xor(Reg(4), Reg(5)), Reg(6)), Reg(8)); 
        %% n = 13
        % Reg(1) = xor(xor(xor(Reg(1), Reg(3)), Reg(4)), Reg(13));
        %%  仅适用于n= 4  情况
        % Reg(1)=xor(Reg(n-1),Reg(n));  % Modulo 2 adder   第一个寄存器 
        j=2;
        while j<=n                  % Registers shift                
             temp1=Reg(j);
             Reg(j)=temp;
             j=j+1;
             temp=temp1;
        end        
        Out(i)=Reg(n);
        if (Reg(n)==1)  Out(i) =-a; end
        if (Reg(n)==0)  Out(i)=a; end
end   

% stairs(Out); axis([1 L -0.5 1.5]);  % Plot the PRBS signal

end