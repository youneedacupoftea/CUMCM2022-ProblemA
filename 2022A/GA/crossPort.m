% 
function [y1, y2] = crossPort(x1, x2)  %交叉函数 此处采用单点交叉 输入为俩个个体的x 把俩个个体x即基因交叉
n = numel(x1);
s = randi([1,n-1]);  %n-1是因为最后一个点交换没意义
y1 = [x1(1 : s) x2(s+1 : end)];  %直接用矩阵把俩个部分拼接起来 从s+1 开始交叉
y2 = [x2(1 : s) x1(s+1 : end)];
end

% %% 交叉函数 基于方向的启发式算子 基于排序选择展开
% function [y1, y2] = crossplus(input, nPop, nVat)   %输入参数为 排序选择后的结果  以及种群规模n 
% a = 1;  %表示步进 一般取1
% R1 = rand(1, nVat);   %表示生成 跟标量维数一致的均匀分布的系数
% R2 = rand(1, nVat);   %表示生成 跟标量维数一致的均匀分布的系数
% for i = 1 : nPop/2
%     D = input(i+nPop/2) - input(i);
%     y1 = input(i) + a*R1.*D;
%     y2 = input(i+nPop/2) + a*R2.*D;
% end
% end


