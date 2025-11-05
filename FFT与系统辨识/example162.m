% Generate a PRBS signal
% 1. 自相关性尽可能为δ函数（低相关性）。  2. 尽可能均匀分布(接近白噪声）  3. 但又要可被重复实现 
clear all;
n=4;                                % the number of registers
L=15+1;                             % the length of M-sequence to be required   越长越接近白噪声

for i = 1:n  Reg(i) = 1; end        % Initialize n registers  位移寄存器全1
Out(1)=Reg(n);

for i=2:L
        temp=Reg(1);
        Reg(1)=xor(Reg(n-1),Reg(n));  % Modulo 2 adder   第一个寄存器 
        j=2;
        while j<=n                  % Registers shift                
             temp1=Reg(j);
             Reg(j)=temp;
             j=j+1;
             temp=temp1;
        end        
        Out(i)=Reg(n);
end   

stairs(Out); axis([1 L -0.5 1.5]);  % Plot the PRBS signal

% plot(Out)