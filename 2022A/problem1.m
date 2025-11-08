clc,clear,close all;  

% import data
data3 = xlsread('D:/学习资料/数学建模/2022国赛题目/A题/附件3.xlsx', 1); % 读取数据 
data4 = xlsread('D:/学习资料/数学建模/2022国赛题目/A题/附件4.xlsx', 1);
data_1 = data3(1,:);
w = data_1(1);  % 入射波浪频率
T = 2*pi/w*40;  % 一个波浪周期
mr = data_1(3); % 纵摇附加转动惯量
F1 = data_1(6); % 垂荡激励力振幅
F2 = data_1(7); % 纵摇激励力矩振幅
mf = data4(1);  % 浮子质量
mz = data4(5);  % 振子质量
m = data_1(2);  % 垂荡附加质量
h1 = data4(3);  % 浮子圆柱部分高度
h2 = data4(4);  % 浮子圆锥部分高度
rz = data4(6);  % 振子半径
rf = data4(2);  % 浮子底半径
hz = data4(7);  % 振子高度
tho = data4(8); % 海水密度
g = data4(9);   % 重力加速度
orilen = data4(11); % 弹簧原长
k1 = data_1(4); % 垂荡兴波阻尼系数
k2 = data_1(5); % 纵摇兴波阻尼系数
k3 = data4(10); % 弹簧刚度
k4 = data4(12);  % 扭转弹簧刚度
k5 = data4(13);  % 静水恢复力矩系数
delta_t = 0.001;   % 步进0.001s

% 计算一些静态参数 以海平面作为原点 向下为正方向
Vzhui = 1/3*pi*rf^2*h2; % 圆锥壳的体积
Vpai_0 = (mz+mf)/tho; % 静态平衡时 浮子排开水的体积
Hf = (Vpai_0-Vzhui)/pi/(rf^2)+h2;  % 静态平衡时 浮子沉在海底的长度 浮子顶点深度
len0 = mz*g/k3;  % 静态时 弹簧的压缩量
len1 = h2+orilen+hz/2; % 弹簧原长时的初始距离
len2 = h2+hz/2;  % 俩者之间的最小间距
len3 = h2+hz/2+2*orilen;  % 俩者之间最大间距
Hz = Hf-hz/2-h2-(orilen-len0); % 静态时 振子的高度 以中心点为研究对象

% 调用目标函数 mode=1 即第一种情况 阻尼系数为定值
[VF_out,XF_out,AF_out,VZ_out,XZ_out,AZ_out] = Get_Endfu(delta_t,T,Hf,Hz,data_1,data4,len1,len2,len3,1);
[VF_out1,XF_out1,AF_out1,VZ_out1,XZ_out1,AZ_out1] = Get_Endfu(delta_t,T,Hf,Hz,data_1,data4,len1,len2,len3,1.5);

% % 写入结果
k = 0;   % 计次变量
n = round(T/delta_t)+1;
index = zeros(1,897); %间隔0.2s的索引值数组
Time_write = zeros(1,897);
for i = 2 : n
    t_now = (i-1)*delta_t;
    if mod(t_now,0.2) == 0
        k = k+1;
        index(k) = i;
        Time_write(k) = t_now;
    end
    if t_now == 10
        fprintf("阻尼系数为常数时:\n");
        fprintf("第%d秒时，浮子的垂荡位移为:%.4f\n",t_now,XF_out(i));
        fprintf("第%d秒时，浮子的速度为:%.4f\n",t_now,VF_out(i));
        fprintf("第%d秒时，振子的垂荡位移为:%.4f\n",t_now,XZ_out(i));
        fprintf("第%d秒时，振子的速度为:%.4f\n",t_now,VZ_out(i));
        fprintf("阻尼系数变化时:\n");
        fprintf("第%d秒时，浮子的垂荡位移为:%.4f\n",t_now,XF_out1(i));
        fprintf("第%d秒时，浮子的速度为:%.4f\n",t_now,VF_out1(i));
        fprintf("第%d秒时，振子的垂荡位移为:%.4f\n",t_now,XZ_out1(i));
        fprintf("第%d秒时，振子的速度为:%.4f\n",t_now,VZ_out1(i));
    end
     if t_now == 20
        fprintf("阻尼系数为常数时:\n");
        fprintf("第%d秒时，浮子的垂荡位移为:%.4f\n",t_now,XF_out(i));
        fprintf("第%d秒时，浮子的速度为:%.4f\n",t_now,VF_out(i));
        fprintf("第%d秒时，振子的垂荡位移为:%.4f\n",t_now,XZ_out(i));
        fprintf("第%d秒时，振子的速度为:%.4f\n",t_now,VZ_out(i));
        fprintf("阻尼系数变化时:\n");
        fprintf("第%d秒时，浮子的垂荡位移为:%.4f\n",t_now,XF_out1(i));
        fprintf("第%d秒时，浮子的速度为:%.4f\n",t_now,VF_out1(i));
        fprintf("第%d秒时，振子的垂荡位移为:%.4f\n",t_now,XZ_out1(i));
        fprintf("第%d秒时，振子的速度为:%.4f\n",t_now,VZ_out1(i));
     end
     if t_now == 40
        fprintf("阻尼系数为常数时:\n");
        fprintf("第%d秒时，浮子的垂荡位移为:%.4f\n",t_now,XF_out(i));
        fprintf("第%d秒时，浮子的速度为:%.4f\n",t_now,VF_out(i));
        fprintf("第%d秒时，振子的垂荡位移为:%.4f\n",t_now,XZ_out(i));
        fprintf("第%d秒时，振子的速度为:%.4f\n",t_now,VZ_out(i));
        fprintf("阻尼系数变化时:\n");
        fprintf("第%d秒时，浮子的垂荡位移为:%.4f\n",t_now,XF_out1(i));
        fprintf("第%d秒时，浮子的速度为:%.4f\n",t_now,VF_out1(i));
        fprintf("第%d秒时，振子的垂荡位移为:%.4f\n",t_now,XZ_out1(i));
        fprintf("第%d秒时，振子的速度为:%.4f\n",t_now,VZ_out1(i));
     end
     if t_now == 60
        fprintf("阻尼系数为常数时:\n");
        fprintf("第%d秒时，浮子的垂荡位移为:%.4f\n",t_now,XF_out(i));
        fprintf("第%d秒时，浮子的速度为:%.4f\n",t_now,VF_out(i));
        fprintf("第%d秒时，振子的垂荡位移为:%.4f\n",t_now,XZ_out(i));
        fprintf("第%d秒时，振子的速度为:%.4f\n",t_now,VZ_out(i));
        fprintf("阻尼系数变化时:\n");
        fprintf("第%d秒时，浮子的垂荡位移为:%.4f\n",t_now,XF_out1(i));
        fprintf("第%d秒时，浮子的速度为:%.4f\n",t_now,VF_out1(i));
        fprintf("第%d秒时，振子的垂荡位移为:%.4f\n",t_now,XZ_out1(i));
        fprintf("第%d秒时，振子的速度为:%.4f\n",t_now,VZ_out1(i));
     end
     if t_now == 100
        fprintf("阻尼系数为常数时:\n");
        fprintf("第%d秒时，浮子的垂荡位移为:%.4f\n",t_now,XF_out(i));
        fprintf("第%d秒时，浮子的速度为:%.4f\n",t_now,VF_out(i));
        fprintf("第%d秒时，振子的垂荡位移为:%.4f\n",t_now,XZ_out(i));
        fprintf("第%d秒时，振子的速度为:%.4f\n",t_now,VZ_out(i));
        fprintf("阻尼系数变化时:\n");
        fprintf("第%d秒时，浮子的垂荡位移为:%.4f\n",t_now,XF_out1(i));
        fprintf("第%d秒时，浮子的速度为:%.4f\n",t_now,VF_out1(i));
        fprintf("第%d秒时，振子的垂荡位移为:%.4f\n",t_now,XZ_out1(i));
        fprintf("第%d秒时，振子的速度为:%.4f\n",t_now,VZ_out1(i));
     end  
end
VF_write = VF_out(index);
XF_write = XF_out(index);
VZ_write = VZ_out(index);
XZ_write = XZ_out(index);
VF_write1 = VF_out1(index);
XF_write1 = XF_out1(index);
VZ_write1 = VZ_out1(index);
XZ_write1 = XZ_out1(index);

% xlswrite('D:/学习资料/数学建模/2022国赛题目/A题/result1-1.xlsx',Time_write', ...
%   1, 'A3');
% xlswrite('D:/学习资料/数学建模/2022国赛题目/A题/result1-1.xlsx',XF_write', ...
%   1, 'B3');
% xlswrite('D:/学习资料/数学建模/2022国赛题目/A题/result1-1.xlsx',VF_write', ...
%   1, 'C3');
% xlswrite('D:/学习资料/数学建模/2022国赛题目/A题/result1-1.xlsx',XZ_write', ...
%   1, 'D3');
% xlswrite('D:/学习资料/数学建模/2022国赛题目/A题/result1-1.xlsx',VZ_write', ...
%   1, 'E3');
% 
% xlswrite('D:/学习资料/数学建模/2022国赛题目/A题/result1-2.xlsx',Time_write', ...
%   1, 'A3');
% xlswrite('D:/学习资料/数学建模/2022国赛题目/A题/result1-2.xlsx',XF_write1', ...
%   1, 'B3');
% xlswrite('D:/学习资料/数学建模/2022国赛题目/A题/result1-2.xlsx',VF_write1', ...
%   1, 'C3');
% xlswrite('D:/学习资料/数学建模/2022国赛题目/A题/result1-2.xlsx',XZ_write1', ...
%   1, 'D3');
%xlswrite('D:/学习资料/数学建模/2022国赛题目/A题/result1-2.xlsx',VZ_write1', ...
%   1, 'E3');


% 画图 检测结果 情况1
t = linspace(0,T,n);
figure(1)
plot(t, VF_out);
xlabel("时间轴");
ylabel("浮子的速度");
title("浮子速度和时间的关系");
figure(2)
plot(t, XF_out);
xlabel("时间轴");
ylabel("浮子的位移");
title("浮子位移和时间的关系");
figure(3)
plot(t, VZ_out);
xlabel("时间轴");
ylabel("振子的速度");
title("振子速度和时间的关系");
figure(4)
plot(t, XZ_out);
xlabel("时间轴");
ylabel("振子的位移");
title("振子位移和时间的关系");

% 画图 检测结果 情况2
t = linspace(0,T,n);
figure(5)
plot(t, VF_out1);
xlabel("时间轴");
ylabel("浮子的速度");
title("浮子速度和时间的关系");
figure(6)
plot(t, XF_out1);
xlabel("时间轴");
ylabel("浮子的位移");
title("浮子位移和时间的关系");
figure(7)
plot(t, VZ_out1);
xlabel("时间轴");
ylabel("振子的速度");
title("振子速度和时间的关系");
figure(8)
plot(t, XZ_out1);
xlabel("时间轴");
ylabel("振子的位移");
title("振子位移和时间的关系");
find(VZ_out==max(VZ_out))
find(VF_out==max(VF_out))
find(VZ_out==min(VZ_out))
find(VF_out==min(VF_out))



%% 功能函数和目标函数
function Vpai = Get_Vpai(xf)  % xf为浮子顶点的坐标 即现在在水里的深度
   h0 = 0.8; % 浮子圆锥部分的高度
   r0 = 1;
   if xf>=0 && xf<=h0
       Vpai = (xf/(0.8))^2*pi/3*xf;
   elseif xf>h0 && xf<= 3.8
       Vpai = (h0/(0.8))^2*pi/3*h0 + (xf-h0)*pi*r0^2;
   else
       Vpai = (h0/(0.8))^2*pi/3*h0 + (3.8-h0)*pi*r0^2;    
   end
end


function [Endvf,Endxf,Endaf,Endvz,Endxz,Endaz] = Get_Endfu(delta_t,T,Hf,Hz,data_1,data4,len1,len2,len3,mode)  
  % 计算浮子沿着步进的加速度 速度 和位移  Hf为静态参数
  % 输出为 振子的速度 位移 浮子的速度位移
  % 先读入题干参数
   F1 = data_1(6); % 垂荡激励力振幅
   mf = data4(1);  % 浮子质量
   m = data_1(2);  % 垂荡附加质量
   n = round(T/delta_t)+1;   % 计算数组的大小
   afu = zeros(1,n);
   Vfu = zeros(1,n);
   Xfu = zeros(1,n);
   azhen = zeros(1,n);
   Vzhen = zeros(1,n);
   Xzhen = zeros(1,n);
   % 先进行参数初始化
   afu(1) = F1*cos(0)/(mf+m);
   Vfu(1) = 0;
   Xfu(1) = Hf;
   azhen(1) = 0;
   Vzhen(1) = 0;
   Xzhen(1) = Hz;
   for i = 2 : n
       time_now = (i-1)*delta_t;  % 时间从0.2s开始
       Vfu(i) = afu(i-1)*delta_t+Vfu(i-1);  % v=a*t
       Xfu(i) = (Vfu(i)+Vfu(i-1))/2*delta_t+Xfu(i-1); % 平均速度算位移
       if i == 2
           Vzhen(2) = 0;
           Xzhen(2) = Xzhen(1);
       else
           Vzhen(i) = azhen(i-1)*delta_t+Vzhen(i-1);
           Xzhen(i) = (Vzhen(i)+Vzhen(i-1))/2*delta_t+Xzhen(i-1);
       end
       azhen(i) = Get_azhen(Xfu(i),Xzhen(i),Vfu(i),Vzhen(i),len1,len2,len3,mode);
       afu(i) =  Get_afu(Xfu(i),Xzhen(i),Vfu(i),Vzhen(i),len1,time_now,len2,len3,mode);
   end
   
   Endvf = Vfu-Vfu(1);
   Endxf = Xfu-Xfu(1);
   Endaf = afu-afu(1);
   Endvz = Vzhen-Vzhen(1);
   Endxz = Xzhen-Xzhen(1);
   Endaz = azhen-azhen(1);
end

function Ft = Get_Ft(Xfu,Xzhen,len1,len2,len3)  % 求弹力 以浮子为主要研究对象 输入为某一个时刻的Xfu 和 Xzhen
   length = abs(Xfu-Xzhen);
   if length <= len2
       length = len2;
   elseif length >= len3
       length = len3;
   else
       length = length;
   end
   sub = abs(length-len1);
   if length < len1
       Ft = sub*80000;
   else
       Ft = -sub*80000;
   end
end

function [flag_fu,flag_zhen] = Judge(Vfu,Vzhen)  %通过速度来判断阻尼器给的力的方向
    if Vfu*Vzhen < 0   %俩者方向不同 均为阻力
        if Vfu < 0
            flag_fu = 1;
            flag_zhen = -1;
        else
            flag_fu = -1;
            flag_zhen = 1;
        end
    else   % 俩者方向相同
        if abs(Vfu)>abs(Vzhen)   % 浮子速度大于振子 浮子为阻力 振子为动力
            if Vfu<0
                flag_fu = 1;
                flag_zhen = -1;
            else 
                flag_fu = -1;
                flag_zhen = 1;
            end
        elseif abs(Vfu)<=abs(Vzhen)
            if Vfu<0
                flag_fu = -1;
                flag_zhen = 1;
            elseif Vfu>=0
                flag_fu = 1;
                flag_zhen = -1;
            end 
        end
    end       
end

function afu = Get_afu(Xfu,Xzhen,Vfu,Vzhen,len1,time_now,len2,len3,mode) % 输入为某一个时刻俩者速度 以及Xfu
   m = 1335.535;  % 垂荡附加质量
   mf = 4866;  % 浮子质量
   g = 9.8;   % 重力加速度
   tho = 1025; % 海水密度
   k1 = 656.3616; % 垂荡兴波阻尼系数
   F1 = 6250; % 垂荡激励力振幅
   w = 1.4005;  % 入射波浪频率
   [flag,~]=Judge(Vfu,Vzhen);  % 接收标志位
   Fzn = 10000*flag*(abs(Vfu-Vzhen)^mode);
   Ffu = mf*g-tho*g*Get_Vpai(Xfu)-Vfu*k1+Fzn...
   +Get_Ft(Xfu,Xzhen,len1,len2,len3)+F1*cos(w*time_now);
   afu = Ffu/(m+mf);
end

function azhen = Get_azhen(Xfu,Xzhen,Vfu,Vzhen,len1,len2,len3,mode) % 输入为某一个时刻俩者速度 以及Xfu
   mz = 2433;  % 振子质量
   g = 9.8;   % 重力加速度
   [~,flag]=Judge(Vfu,Vzhen);  % 接收标志位
   Fzn_0 = 10000*(abs(Vfu-Vzhen)^mode);
   Fzn = 10000*(abs(Vfu-Vzhen)^mode)*flag;
   if Vzhen == 0
       Ffu = mz*g-Get_Ft(Xfu,Xzhen,len1,len2,len3)+Fzn_0;
   else
       Ffu = mz*g-Get_Ft(Xfu,Xzhen,len1,len2,len3)+Fzn;
   end
   azhen = Ffu/mz;
end


