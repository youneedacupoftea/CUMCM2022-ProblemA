%蒙特卡洛求排水体积
%角度单位均为弧度制

function V=MCvolume(h,h0,angle)
   a1=atan(1.25);
   a2=pi/2-a1-angle;
   x=sin(a2)*(1+0.64)^(1/2);
   x1=h+h0-x;
   a=x1/cos(angle);%半圆柱一个边
   x2=cos(atan(1.25)-angle)*(1+0.64)^(1/2);
   x3=h+h0-x2;
   b=x3/cos(angle);%半圆柱一个边
   V1=(1/3)*0.8*pi*1;%圆锥体积
   %切面方程
   
end