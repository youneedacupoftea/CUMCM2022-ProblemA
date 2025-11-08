%% 交叉函数 基于方向的启发式算子 基于排序选择展开
function [Offspring] = crossplus(Parent, nPop, nVat)   %输入参数为 排序选择后的结果  以及种群规模nPop
a = 1;  %表示步进 一般取1
R1 = rand(1, nVat);   %表示生成 跟标量维数一致的均匀分布的系数
R2 = rand(1, nVat);   %表示生成 跟标量维数一致的均匀分布的系数
for i = 1 : nPop/2
    D =  Parent(i).x - Parent(i+nPop/2).x;     %求最大值 就是优向量减去不优的向量
    Offspring(i, 1).x = Parent(i).x + a*R1.*D;
    Offspring(i, 2).x = Parent(i+nPop/2).x + a*R2.*D;
end
end

