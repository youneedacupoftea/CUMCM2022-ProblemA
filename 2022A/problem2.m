clc,clear,close all;  

% import data
data3 = xlsread('D:/学习资料/数学建模/2022国赛题目/A题/附件3.xlsx', 1); % 读取数据 
data4 = xlsread('D:/学习资料/数学建模/2022国赛题目/A题/附件4.xlsx', 1);
data_1 = data3(2,:);
w = data_1(1);  % 入射波浪频率
T = 2*pi/w*40;     % 一个波浪周期
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
% [VF_out,XF_out,AF_out,VZ_out,XZ_out,AZ_out,P_out] = Get_Endfu(delta_t,T,Hf,Hz,data_1,data4,len1,len2,len3,1);
% [VF_out1,XF_out1,AF_out1,VZ_out1,XZ_out1,AZ_out1,P_out1] = Get_Endfu(delta_t,T,Hf,Hz,data_1,data4,len1,len2,len3,0.5);
% 调用步进可变的遍历搜索
% temp = -inf;
% Y_Show = zeros(1,10000);
% mm=1;
% x_show = linspace(1,100000,100000);
% for i = 1:100000
%    [~,~,~,~,~,~,P_out] = Get_Endfu(delta_t,T,Hf,Hz,data_1,data4,len1,len2,len3,1,x_show(i));
%    Y_Show(mm) = P_out;
%    mm = mm + 1;
%    if P_out > temp
%        temp = P_out;
%        k = i;
%    end
%    fprintf("i=%d\n",i);
% end
%  [~,~,~,~,~,~,P_out] = Get_Endfu(delta_t,T,Hf,Hz,data_1,data4,len1,len2,len3,1,43062);
%  fprintf("找到的最大功率为:%.2f\n",P_out);
%  fprintf("找到的最大功率为:%.2f\n",temp);
%  fprintf("最合适的阻尼系数为:%.2f\n",k);
% plot(x_show,Y_Show);
% 找到的最大功率为:3544533.65
% 最合适的阻尼系数为:43217.00
 
% 调用随机搜索函数进行求解
%    s=0;i=1;N=1000;
%    while(s<N)
%        kzn = 2*rand+43063;
%        [~,~,~,~,~,~,P_out] = Get_Endfu(delta_t,T,Hf,Hz,data_1,data4,len1,len2,len3,1,kzn);
%        e(i)=P_out;
%        maxe1(i) = max(e);
%        if i == 1
%           emax=maxe1(1); 
%           kzn0 = kzn;
%        else
%            if(maxe1(i)>maxe1(i-1))
%               emax=maxe1(i);  
%               kzn0 = kzn;
%               s = 0;
%            else
%                s=s+1;
%             end 
%        end
%        i = i + 1;
% %        fprintf("迭代次数为:%d \n",i);
%    end
%   fprintf("找到的最大功率为:%.4f\n",emax);
%   fprintf("最合适的阻尼系数为:%.4f\n",kzn0);
%   figure(1)
%   plot(e);  
%   xlabel("迭代次数");
%   ylabel("目标函数值(W)");
%   title("目标函数关于迭代次数的变化");

%第二问的第二小问，步进搜素求解
% maxP = -inf;
% mode=1:0.1:2;
% kzn =0:10:100000;
% k=1;
%  [~,~,~,~,~,~,P_out]= Get_Endfu(delta_t,T,Hf,Hz,data_1,data4,len1,len2,len3,1,43064)
%  [~,~,~,~,~,~,P_out]= Get_Endfu(delta_t,T,Hf,Hz,data_1,data4,len1,len2,len3,1,40956)
% for j=1:10^4
%    for i = 1:11
%          [~,~,~,~,~,~,P_out]= Get_Endfu(delta_t,T,Hf,Hz,data_1,data4,len1,len2,len3,mode(i),kzn(j));
%          if maxP<P_out
%              maxP=P_out;
%              best_mode=mode(i)-1;
%              best_kzn=kzn(j);
%          end   
%          x(k)=mode(i)-1;
%          y(k)=kzn(j);
%          z(k)=P_out;
%          k=k+1; 
%     end
%     fprintf("j=%d\n",j);
% end
%  fprintf("找到的最大功率为:%.4f\n",maxP);
%  fprintf("最合适的阻尼系数为:%.4f\n",best_kzn);
%  fprintf("最合适的幂系数为:.%4f\n",best_mode);
%  plot3(x,y,z);
%  xlabel("幂系数");
%  ylabel("阻尼系数");
%  zlabel("功率");
 
%第二问的第二小问，随机搜索算法求解
   s=0;i=1;N=1000;
   while(s<N)
       mode = rand*0.2+0.3;
       kzn = rand*20+94990;
       [~,~,~,~,~,~,P_out]= Get_Endfu(delta_t,T,Hf,Hz,data_1,data4,len1,len2,len3,mode+1,kzn);
       e1(i)=P_out;
       maxe1(i) = max(e1);
       if i == 1
          emax1=maxe1(1); 
          kzn1 = kzn;
          mode0 = mode;
       else
           if(maxe1(i)>maxe1(i-1))
              emax1=maxe1(i);  
              kzn1 = kzn;
              mode0 = mode;
              s = 0;
           else
               s=s+1;
            end 
       end
       i = i + 1;
       fprintf("迭代次数为：%d\n",i);
   end
  fprintf("找到的最大功率为:%.4f\n",emax1);
  fprintf("最合适的阻尼系数为:%.4f\n",kzn1);
  fprintf("最合适的幂系数为:.%4f\n",mode0);
  %figure(2)
  for i =1: length(e1)
      fprintf("P=%4f",e1(i));
  end
  plot(e1);
  xlabel("迭代次数");
  ylabel("输出功率/W");
  %title("目标函数关于迭代次数的变化");

   
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

function [Endvf,Endxf,Endaf,Endvz,Endxz,Endaz,P_out] = Get_Endfu(delta_t,T,Hf,Hz,data_1,data4,len1,len2,len3,mode,kzn)  
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
   index0 = 0;  
   index1 = 0;  
   for i = 2 : n
       time_now = (i-1)*delta_t;  % 时间从0.2s开始
       if time_now == 80
           index0 = i;
       elseif time_now == 100
           index1 = i;
       end
       
       Vfu(i) = afu(i-1)*delta_t+Vfu(i-1);  % v=a*t
       Xfu(i) = (Vfu(i)+Vfu(i-1))/2*delta_t+Xfu(i-1); % 平均速度算位移
       if i == 2
           Vzhen(2) = 0;
           Xzhen(2) = Xzhen(1);
       else
           Vzhen(i) = azhen(i-1)*delta_t+Vzhen(i-1);
           Xzhen(i) = (Vzhen(i)+Vzhen(i-1))/2*delta_t+Xzhen(i-1);
       end
       azhen(i) = Get_azhen(Xfu(i),Xzhen(i),Vfu(i),Vzhen(i),len1,len2,len3,mode,kzn);
       afu(i) =  Get_afu(Xfu(i),Xzhen(i),Vfu(i),Vzhen(i),len1,time_now,len2,len3,mode,kzn);
   end
   
   num = index1-index0;
   W_out = zeros(1,num);
   k = 1;
   for j = index0  : index1-1  %计算总功
       W_out(k) = kzn*(abs(Vfu(j)-Vzhen(j))^mode)*abs((Vfu(j)+Vfu(j+1))/2-(Vzhen(j)+Vzhen(j+1))/2)*delta_t;
       k = k + 1;
   end
   P_out = sum(W_out)/20;
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
    flag_fu = 0;
    flag_zhen = 0;
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
            else
                flag_fu = 1;
                flag_zhen = -1;
            end 
        end
    end  
    
end

function afu = Get_afu(Xfu,Xzhen,Vfu,Vzhen,len1,time_now,len2,len3,mode,kzn) % 输入为某一个时刻俩者速度 以及Xfu
   m = 1165.992;  % 垂荡附加质量
   mf = 4866;  % 浮子质量
   g = 9.8;   % 重力加速度
   tho = 1025; % 海水密度
   k1 = 167.8395; % 垂荡兴波阻尼系数
   F1 = 4890; % 垂荡激励力振幅
   w = 2.2143;  % 入射波浪频率
   [flag,~]=Judge(Vfu,Vzhen);  % 接收标志位
   Fzn = kzn*flag*(abs(Vfu-Vzhen)^mode);
   Ffu = mf*g-tho*g*Get_Vpai(Xfu)-Vfu*k1+Fzn...
   +Get_Ft(Xfu,Xzhen,len1,len2,len3)+F1*cos(w*time_now);
   afu = Ffu/(m+mf);
end

function azhen = Get_azhen(Xfu,Xzhen,Vfu,Vzhen,len1,len2,len3,mode,kzn) % 输入为某一个时刻俩者速度 以及Xfu
   mz = 2433;  % 振子质量
   g = 9.8;   % 重力加速度
   [~,flag]=Judge(Vfu,Vzhen);  % 接收标志位
   Fzn_0 = kzn*(abs(Vfu-Vzhen)^mode);
   Fzn = kzn*(abs(Vfu-Vzhen)^mode)*flag;
   if Vzhen == 0
       Ffu = mz*g-Get_Ft(Xfu,Xzhen,len1,len2,len3)+Fzn_0;
   else
       Ffu = mz*g-Get_Ft(Xfu,Xzhen,len1,len2,len3)+Fzn;
   end
   azhen = Ffu/mz;
end