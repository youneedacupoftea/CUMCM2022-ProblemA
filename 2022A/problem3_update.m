clc,clear,close all;  

% import data
data3 = xlsread('D:/学习资料/数学建模/2022国赛题目/A题/附件3.xlsx', 1); % 读取数据 
data4 = xlsread('D:/学习资料/数学建模/2022国赛题目/A题/附件4.xlsx', 1);
data_1 = data3(3,:);
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
delta_t = 0.0001;   % 步进0.001s

% 计算一些静态参数 以海平面作为原点 向下为正方向
Vzhui = 1/3*pi*rf^2*h2; % 圆锥壳的体积
Vpai_0 = (mz+mf)/tho; % 静态平衡时 浮子排开水的体积
Hf = (Vpai_0-Vzhui)/pi/(rf^2);  % 静态平衡时 浮子沉在海底的长度 转轴连接点的坐标
len0 = mz*g/k3;  % 静态时 弹簧的压缩量
len1 = orilen+hz/2; % 弹簧原长时的初始距离
len2 = hz/2;  % 俩者之间的最小间距
len3 = hz/2+2*orilen;  % 俩者之间最大间距
Hz = Hf-hz/2-(orilen-len0); % 静态时 振子的高度 以中心点为研究对象 未发生变化 跟前面情况一致
tho_0 = abs(Hz-Hf);    %xtho的初始长度

% 调用函数 
[xfc,xfz,vfc,vfz,W_show,Y_show,V_show,Xita_show] = Get_Endfu(delta_t,T,Hf,Hz,len1,10000,1000);

% % 写入结果
k = 0;   % 计次变量
m = round(T/delta_t)+1;
index = zeros(1,m); %间隔0.2s的索引值数组
Time_write = zeros(1,m);
for i = 2 : m
    t_now = (i-1)*delta_t;
    if mod(t_now,0.2) == 0
        k = k+1;
        index(k) = i;
        Time_write(k) = t_now;
    end
    if t_now == 10
        fprintf("对于振子来说:\n");
        %fprintf("第%d秒时，垂荡位移为:%.4f\n",t_now,Y_show(i));
        fprintf("第%d秒时，速度为:%.4f\n",t_now,V_show(i));
        fprintf("第%d秒时，角度位移为:%.4f\n",t_now,Xita_show(i));
        fprintf("第%d秒时，角速度为:%.4f\n",t_now,W_show(i));
        fprintf("对于浮子来说:\n");
        fprintf("第%d秒时, 垂荡位移为:%.4f\n",t_now,xfc(i));
        fprintf("第%d秒时，速度为:%.4f\n",t_now,vfc(i));
        fprintf("第%d秒时，角度位移为:%.4f\n",t_now,xfz(i));
        fprintf("第%d秒时，角速度为:%.4f\n",t_now,vfz(i));
    end
     if t_now == 20
         fprintf("对于振子来说:\n");
        fprintf("第%d秒时，垂荡位移为:%.4f\n",t_now,Y_show(i));
        fprintf("第%d秒时，速度为:%.4f\n",t_now,V_show(i));
        fprintf("第%d秒时，角度位移为:%.4f\n",t_now,Xita_show(i));
        fprintf("第%d秒时，角速度为:%.4f\n",t_now,W_show(i));
        fprintf("对于浮子来说:\n");
        fprintf("第%d秒时, 垂荡位移为:%.4f\n",t_now,xfc(i));
        fprintf("第%d秒时，速度为:%.4f\n",t_now,vfc(i));
        fprintf("第%d秒时，角度位移为:%.4f\n",t_now,xfz(i));
        fprintf("第%d秒时，角速度为:%.4f\n",t_now,vfz(i));
     end
     if t_now == 40
        fprintf("对于振子来说:\n");
        fprintf("第%d秒时，垂荡位移为:%.4f\n",t_now,Y_show(i));
        fprintf("第%d秒时，速度为:%.4f\n",t_now,V_show(i));
        fprintf("第%d秒时，角度位移为:%.4f\n",t_now,Xita_show(i));
        fprintf("第%d秒时，角速度为:%.4f\n",t_now,W_show(i));
        fprintf("对于浮子来说:\n");
        fprintf("第%d秒时, 垂荡位移为:%.4f\n",t_now,xfc(i));
        fprintf("第%d秒时，速度为:%.4f\n",t_now,vfc(i));
        fprintf("第%d秒时，角度位移为:%.4f\n",t_now,xfz(i));
        fprintf("第%d秒时，角速度为:%.4f\n",t_now,vfz(i));
     end
     if t_now == 60
        fprintf("对于振子来说:\n");
        fprintf("第%d秒时，垂荡位移为:%.4f\n",t_now,Y_show(i));
        fprintf("第%d秒时，速度为:%.4f\n",t_now,V_show(i));
        fprintf("第%d秒时，角度位移为:%.4f\n",t_now,Xita_show(i));
        fprintf("第%d秒时，角速度为:%.4f\n",t_now,W_show(i));
        fprintf("对于浮子来说:\n");
        fprintf("第%d秒时, 垂荡位移为:%.4f\n",t_now,xfc(i));
        fprintf("第%d秒时，速度为:%.4f\n",t_now,vfc(i));
        fprintf("第%d秒时，角度位移为:%.4f\n",t_now,xfz(i));
        fprintf("第%d秒时，角速度为:%.4f\n",t_now,vfz(i));
     end
     if t_now == 100
        fprintf("对于振子来说:\n");
        fprintf("第%d秒时，垂荡位移为:%.4f\n",t_now,Y_show(i));
        fprintf("第%d秒时，速度为:%.4f\n",t_now,V_show(i));
        fprintf("第%d秒时，角度位移为:%.4f\n",t_now,Xita_show(i));
        fprintf("第%d秒时，角速度为:%.4f\n",t_now,W_show(i));
        fprintf("对于浮子来说:\n");
        fprintf("第%d秒时, 垂荡位移为:%.4f\n",t_now,xfc(i));
        fprintf("第%d秒时，速度为:%.4f\n",t_now,vfc(i));
        fprintf("第%d秒时，角度位移为:%.4f\n",t_now,xfz(i));
        fprintf("第%d秒时，角速度为:%.4f\n",t_now,vfz(i));
     end  
end

% 画图 检测结果 情况1
t = linspace(0,T,m);
figure(1)
plot(t, xfc);
xlabel("时间轴");
ylabel("浮子的垂荡位移");

figure(2)
plot(t, xfz);
xlabel("时间轴");
ylabel("浮子的偏移角");

figure(3)
plot(t, vfc);
xlabel("时间轴");
ylabel("浮子的垂荡速度");

figure(4)
plot(t, vfz);
xlabel("时间轴");
ylabel("浮子的纵摇速度");

figure(5)
plot(t, Y_show);
xlabel("时间轴");
ylabel("浮子的垂荡位移");

figure(6)
plot(t, Xita_show);
xlabel("时间轴");
ylabel("振子的偏移角");

figure(7)
plot(t, V_show);
xlabel("时间轴");
ylabel("浮子的垂荡速度");

figure(8)
plot(t, W_show);
xlabel("时间轴");
ylabel("浮子的角速度");




%% 功能函数和目标函数
function Vpai = Get_Vpai(Xfu,xitaF)  % Xfu 为连接点的坐标 xitaZ为浮子偏移角度
   xitaF = pi/2-xitaF;
   X_head = Xfu+0.8*cos(abs(xitaF));  % 浮子顶点坐标 纵坐标
   len_temp = sqrt(1+0.64);           % 圆锥母线的长度
   angle_head = atan(1.25);          %顶角的一半
   alpha = pi/2-abs(xitaF)-angle_head; 
   h_temp = len_temp*abs(sin(alpha));
   h_temp1 = X_head-h_temp;
   long = h_temp1/abs(cos(xitaF));  % 梯形的长边
   short = long-abs(2*tan(xitaF));  % 梯形的短边
   Vpai =1/3*pi*1+pi*1*long-(long-short)*pi*1/2;
end

function [xfc,xfz,vfc,vfz,W_out,Y_out,V_out,Xita_out] = Get_Endfu(delta_t,T,Hf,Hz,len1,kn,kn1)  
  % 计算浮子沿着步进的加速度 速度 和位移  Hf为静态参数
  % 输出为 振子的速度 位移 浮子的速度位移
  % 先读入题干参数
   n = round(T/delta_t)+1;   % 计算数组的大小
   tho_0 = abs(Hz-Hf);    %xtho的初始长度
   
   afu = zeros(1,n);
   Vfu = zeros(1,n);
   Xfu = zeros(1,n);
   
   xitaF = zeros(1,n);
   bfu = zeros(1,n);
   Wfu = zeros(1,n);

   xitaZ = zeros(1,n);
   bzhen = zeros(1,n);
   Wzhen = zeros(1,n);
  
   xtho = zeros(1,n);
   Vtho = zeros(1,n);
   atho = zeros(1,n);
   
   loca = zeros(2,n);      %储存极坐标原点信息 俩维度信息 1表示横坐标 2表示纵坐标
   loca_zhen = zeros(2,n);  %无需初值
   delta_y = zeros(1,n);    %储存特定时刻对应的振子的垂荡位移 要初值
   delta_v = zeros(1,n);    %储存特定时刻的垂荡速度 无初值
   my_w = zeros(1,n);       %纵摇方向角速度 无初值
   my_angle = zeros(1,n);   %纵摇方向角位移 要初值
   
   my_angle(1) = pi/2;
   my_w(1) = 0;
   delta_v(1) = 0;
   delta_y(1) = 0;
   loca(1,1)= 0;  % 初始研究点横坐标为零
   loca(2,1)= Xfu(1); % 初始研究点纵坐标为Hf 相当于loca(:,1)即初始时刻坐标已经知道
 
   
  % 先进行参数初始化
  afu(1) = 1690*cos(0)/(4866+1028.876);
  I = 4866*2^2; %转动惯量 mr^2
  Iadd = 7001.914;  %纵摇附加转动惯量
  Vfu(1) = 0;
  Xfu(1) = round(Hf,3); %初始浮子研究点的坐标
  xitaZ(1) = pi/2;
  xitaF(1) = pi/2;
  bfu(1) = 3640*cos(0)/(I+Iadd);
  bzhen(1) = 0;
  Wzhen(1) = 0;
  Wfu(1) = 0;
  xtho(1) = round(tho_0,3);  %保留三位小数 初始径向长度
  Vtho(1) = 0;  % 初始径向速度
  atho(1) = 0;
 
  
   for i = 2 : n
       time_now = (i-1)*delta_t;    %现在的实际时间
       Vfu(i)=afu(i-1)*delta_t+Vfu(i-1);
       Xfu(i)=(Vfu(i)+Vfu(i-1))/2*delta_t+Xfu(i-1);
       Wfu(i)=Wfu(i-1)+bfu(i-1)*delta_t;
       xitaF(i)=xitaF(i-1)+(Wfu(i)+Wfu(i-1))/2*delta_t;
       if i == 2
           Vtho(2)=0;
           Wzhen(2)=0;
           xtho(2)=xtho(1);
           xitaZ(2)=xitaZ(1);
       else
           Wzhen(i)=Wzhen(i-1)*bzhen(i-1)*delta_t;
           xitaZ(i)=xitaZ(i-1)+(Wzhen(i)+Wzhen(i-1))/2*delta_t;
           Vtho(i)=Vtho(i-1)+atho(i-1)*delta_t;
           xtho(i)=xtho(i-1)+(Vtho(i-1)+Vtho(i))/2*delta_t;
       end
       afu(i)=Get_afu(Xfu(i),Vfu(i),len1,time_now,kn,xitaF(i),xtho(i),xitaZ(i),Vtho(i),Wzhen(i));
       bfu(i)=Get_bfu(Hf,Wfu(i),xitaF(i),len1,xtho(i),xitaZ(i),kn1,Wzhen(i),kn,time_now,Vtho(i));
       bzhen(i)=Get_bzhen(xtho(i),Wzhen(i),Wfu(i),xitaZ(i),xitaF(i),kn1);
       atho(i)= Get_atho(len1,kn,Vtho(i),xitaZ(i),Wzhen(i),xtho(i),Wfu(i),Hf,xitaF(i));
       
       loca(1,i)=loca(1,i-1)+2*(cos(xitaF(i))-cos(xitaF(i-1)));                              % 横坐标
       loca(2,i)=loca(2,i-1)+Xfu(i)-Xfu(i-1)+2*(abs(sin(xitaF(i)))-abs(sin(xitaF(i-1))));    % 纵坐标
       loca_zhen(1,i-1)=loca(1,i-1)+ abs(xtho(i))*Judge_5(xitaZ(i))*cos(xitaZ(i));         % 求横坐标
       loca_zhen(2,i-1)=loca(2,i-1)+ abs(xtho(i))*abs(sin(xitaZ(i)));                      %  求纵坐标
       delta_y(i) = loca_zhen(2,i)-loca_zhen(2,i-1)+delta_y(i-1);            %每时刻              %垂荡方向位移
       % 以下为垂荡方向速度 每时刻
        if xitaZ(i) >pi/2   
            delta_v(i) = Vtho(i)*abs(sin(xitaZ(i)))+Wzhen(i)*abs(cos(xitaZ(i)))*abs(xtho(i));   %矢量式 研究垂荡方向的阻尼力
        else
            delta_v(i) = Vtho(i)*abs(sin(xitaZ(i)))-Wzhen(i)*abs(cos(xitaZ(i)))*abs(xtho(i));   %矢量式 研究垂荡方向的阻尼力
        end
        
        my_w(i) = Wzhen(i);   %角速度转移 每时刻
        [temp,~] = cart2pol(loca_zhen(1,i-1),loca_zhen(2,i-1));    %转到极坐标方便比较
        [temp1,~] = cart2pol(loca_zhen(1,i),loca_zhen(2,i));
         my_angle(i)=my_angle(i-1)+temp1-temp;   %角度位移 结果是每一个时刻的值
   end    
   
   
   % 计算振子的相对位移以及速度等参数 记得减去初值
   xfz=xitaF-xitaF(1);
   xfc=Xfu-Xfu(1);
   vfc=Vfu-Vfu(1);
   vfz=Wfu-Wfu(1);
   W_out=my_w;
   Y_out=delta_y;
   V_out =delta_v;
   Xita_out =my_angle;
end



%% 功能函数和目标函数
function flag_judge = Judge_5(xita)  % 判断相对坐标
    if xita > pi/2
        flag_judge = -1;
    else
        flag_judge = 1;
    end       
end

function Ft = Get_Ft(len1,xtho) % 振子采用极坐标 故直接用xtho来进行计算 理解成中心到接触点的距离
       Ft = 0;
       if abs(xtho)< len1  % 正常压缩的时候 看情况进行加减 研究浮子和振子时
           Ft = abs(abs(xtho)-len1)*80000;  % 缩短的时候 为正
       elseif abs(xtho)>=len1
           Ft = -abs(abs(xtho)-len1)*80000;
       end
end

function [flag_fu,flag_zhen] = Judge_1(Vfu,Vtho,Wzhen,xitaZ,xtho)  %普通阻尼器 垂荡方向
   if xitaZ >pi/2  
       Vzhen = Vtho*abs(sin(xitaZ))+Wzhen*abs(cos(xitaZ))*abs(xtho);   %矢量式 研究垂荡方向的阻尼力
   else
       Vzhen = Vtho*abs(sin(xitaZ))-Wzhen*abs(cos(xitaZ))*abs(xtho);   %矢量式 研究垂荡方向的阻尼力 
   end
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
        else
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

function [flag_fu,flag_zhen] = Judge_2(Vtho,Wzhen,xtho,xitaZ,Wfu,Hf,xitaF) % 普通阻尼器 纵摇研究Vh
   flag_fu = 0;
   flag_zhen = 0;
   if xitaZ >pi/2
       Vzhen = -Vtho*abs(cos(xitaZ))+Wzhen*abs(sin(xitaZ))*abs(xtho);
   else
       Vzhen = Vtho*abs(cos(xitaZ))+Wzhen*abs(sin(xitaZ))*abs(xtho);
   end
   Vfu = Wfu*Hf*abs(sin(xitaF));
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

function [flag_fu,flag_zhen] = Judge3(Wzhen,Wfu)  % 判断旋转阻尼力矩的方向
   flag_fu=0;
   flag_zhen=0;
    Vfu = Wfu;
    Vzhen = Wzhen;
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

function [flag_F,flag_Z]=Judge_4(xitaZ,xitaF)  % 判断弹簧扭矩的方向
   flag_F=0;
   flag_Z=0;
   if xitaZ*xitaF < 0
        if xitaZ < 0
            flag_F = -1;
            flag_Z =1;
        else
            flag_F = 1;
            flag_Z = -1;
        end 
   else  % 同向
       if abs(xitaZ)<abs(xitaF)  %对XZ为动力 对XZ为阻力
           if xitaZ < 0
               flag_Z = -1;
               flag_F = 1;
           else
               flag_Z = 1;
               flag_F =-1;
           end
       elseif abs(xitaF)<= abs(xitaZ)
             if xitaZ < 0
               flag_Z = 1;
               flag_F = -1;
           else
               flag_Z = -1;
               flag_F =1;
             end    
       end
   end
end


function afu = Get_afu(Xfu,Vfu,len1,time_now,kzn,xitaF,xtho,xitaZ,Vtho,Wzhen) % 输入为某一个时刻俩者速度 以及Xfu
   m = 1028.876;  % 垂荡附加质量
   mf = 4866;  % 浮子质量
   g = 9.8;   % 重力加速度
   tho = 1025; % 海水密度
   k1 = 683.4558; % 垂荡兴波阻尼系数
   F1 = 3640; % 垂荡激励力振幅
   w = 1.7152;  % 入射波浪频率
   [flag,~]=Judge_1(Vfu,Vtho,Wzhen,xitaZ,xtho); % 接收标志位
   if xitaZ >pi/2  
       Vzhen = Vtho*abs(sin(xitaZ))+Wzhen*abs(cos(xitaZ))*abs(xtho);  %矢量式 研究垂荡方向的阻尼力
   else
       Vzhen = Vtho*abs(sin(xitaZ))-Wzhen*abs(cos(xitaZ))*abs(xtho);  %矢量式 研究垂荡方向的阻尼力 
   end
   Fzn = kzn*flag*(abs(Vfu-Vzhen));
   Ffu = mf*g-tho*g*Get_Vpai(Xfu,xitaF)-Vfu*k1+Fzn...
   +Get_Ft(len1,xtho)*sin(xitaZ)+F1*cos(w*time_now);
   afu = Ffu/(m+mf);
end


function atho = Get_atho(len1,kzn,Vtho,xitaZ,Wzhen,xtho,Wfu,Hf,xitaF) % 计算沿着杆的加速度
   g = 9.8;   % 重力加速度
   mz = 2433;  % 振子质量
   if xitaZ >pi/2  % 求纵摇方向 bfu 里
       Vzhen = -Vtho*abs(cos(xitaZ))+Wzhen*abs(sin(xitaZ))*abs(xtho);
   else
       Vzhen = Vtho*abs(cos(xitaZ))+Wzhen*abs(sin(xitaZ))*abs(xtho);
   end
   Vfu = Wfu*Hf*abs(sin(xitaF));
   [~,flag_fu] = Judge_2(Vtho,Wzhen,xtho,xitaZ,Wfu,Hf,xitaF);
   Fzong = kzn*flag_fu*(abs(Vfu-Vzhen));  
   
   [~,flag]=Judge_1(Vfu,Vtho,Wzhen,xitaZ,xtho); % 接收标志位
   if xitaZ >pi/2  
       Vzhen = Vtho*abs(sin(xitaZ))-Wzhen*abs(cos(xitaZ))*abs(xtho);   %矢量式 研究垂荡方向的阻尼力
   else
       Vzhen = Vtho*abs(sin(xitaZ))+Wzhen*abs(cos(xitaZ))*abs(xtho);   %矢量式 研究垂荡方向的阻尼力 
   end
   Fheng = kzn*flag*(abs(Vfu-Vzhen));  
   Ftotal = sqrt(Fheng^2+Fzong^2);  %总合力大小
   if flag_fu==-1 && flag ==-1
        Ftho = mz*g*abs(sin(xitaZ))-Ftotal-Get_Ft(len1,xtho);
   else
        Ftho = mz*g*abs(sin(xitaZ))+Ftotal-Get_Ft(len1,xtho);
   end
   atho = Ftho/mz;
end


% function azhen = Get_azhen(Vfu,len1,len2,len3,kzn,Vtho,Wzhen,xitaZ,xtho) % 输入为某一个时刻俩者速度 以及Xfu
%    mz = 2433;  % 振子质量
%    g = 9.8;   % 重力加速度
%    [~,flag]=Judge_1(Vfu,Vtho,Wzhen,xitaZ,tho);  % 接收标志位
%    if xitaZ >pi/2  
%        Vzhen = Vtho*abs(sin(xitaZ))+Wzhen*abs(cos(xitaZ))*abs(xtho);   %矢量式 研究垂荡方向的阻尼力
%    else
%        Vzhen = Vtho*abs(sin(xitaZ))-Wzhen*abs(cos(xitaZ))*abs(xtho);   %矢量式 研究垂荡方向的阻尼力 
%    end
%    Fzn_0 = kzn*(abs(Vfu-Vzhen));
%    Fzn = kzn*(abs(Vfu-Vzhen))*flag;
%    if Vzhen == 0
%        Ffu = mz*g-Get_Ft(len1,len2,len3,xtho)*sin(xitaZ)+Fzn_0;
%    else
%        Ffu = mz*g-Get_Ft(len1,len2,len3,xtho)*sin(xitaZ)+Fzn;
%    end
%    azhen = Ffu/mz;
% end

function bfu = Get_bfu(Hf,Wfu,xitaF,len1,xtho,xitaZ,kn1,Wzhen,kn,time_now,Vtho)   % 得到浮子的角加速度 纵摇运动
   mf = 4866;  % 浮子质量
%  I = mf*round(Hf,3)^2; %转动惯量 mr^
   I = 24499.16;  %转动惯量
   Iadd = 7001.914;  %纵摇附加转动惯量
   k1 = 654.3383; %纵摇兴波阻尼系数
   k2 = 8890.7;   % 静水恢复力系数
   F1 = 1690;    %纵摇激励力振幅
   k3 = 250000;  % 扭矩弹簧刚度
   w = 1.7152;   % 波浪频率
   [flag,~]=Judge_2(Vtho,Wzhen,xtho,xitaZ,Wfu,Hf,xitaF);
   [flag1,~] = Judge3(Wzhen,Wfu);
   [flag2,~] = Judge_4(xitaZ,xitaF);
   if xitaZ >pi/2
       Vzhen = -Vtho*abs(cos(xitaZ))+Wzhen*abs(sin(xitaZ))*abs(xtho);
   else
       Vzhen = Vtho*abs(cos(xitaZ))+Wzhen*abs(sin(xitaZ))*abs(xtho);
   end
    Vfu = Wfu*Hf*abs(sin(xitaF));
    Mzn = kn1*abs(Wzhen)*flag1;
    Mt = k3*abs(xitaZ-xitaF)*flag2;  
       if xitaZ > pi/2
           if xitaF == 0
               Mfu = -Wfu*k1-Get_Ft(len1,xtho)*abs(xtho*cos(xitaZ))+flag...
                   *abs(Vzhen-Vfu)*kn*abs(xtho*cos(xitaZ))+Mzn*2+2*Mt+F1*cos(w*time_now);
           else
               Mfu = -Wfu*k1-xitaF*k2-Get_Ft(len1,xtho)*abs(xtho*cos(xitaZ))+flag...
                   *abs(Vzhen-Vfu)*kn*abs(xtho*cos(xitaZ))+Mzn*2+2*Mt+F1*cos(w*time_now);
           end
       else
           if xitaF == 0
               Mfu = -Wfu*k1-abs(Get_Ft(len1,xtho))*abs(xtho*cos(xitaZ))+flag...
                    *abs(Vzhen-Vfu)*kn*abs(xtho*cos(xitaZ))+Mzn*2+2*Mt+F1*cos(w*time_now);
           else
               Mfu = -Wfu*k1-xitaF*k2-abs(Get_Ft(len1,xtho))*abs(xtho*cos(xitaZ))+flag...
                   *abs(Vzhen-Vfu)*kn*abs(xtho*cos(xitaZ))+Mzn*2+2*Mt+F1*cos(w*time_now);   
           end
       end
   bfu = Mfu/(I+Iadd);   % 最后的总角加速度
end


function bzhen = Get_bzhen(xtho,Wzhen,Wfu,xitaZ,xitaF,kn1)
   mz = 2433;  % 振子质量
 %  I = mz*(xtho)^2;
   I = 13000;
   g = 9.8;
   k3 = 250000;  % 扭矩弹簧刚度
   [~,flag1] = Judge3(Wzhen,Wfu);
   [~,flag2] = Judge_4(xitaZ,xitaF);
   Mzn = kn1*abs(Wzhen)*flag1;
   Mt = k3*abs(xitaZ-xitaF)*flag2;  
   if xitaZ > pi/2
       Mzhen = 2*Mzn+2*Mt+mz*g*xtho*abs(cos(xitaZ));
   else
       Mzhen = 2*Mzn+2*Mt-mz*g*xtho*abs(cos(xitaZ));
   end
   bzhen = Mzhen/I;
end