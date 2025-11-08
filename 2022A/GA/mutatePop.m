function p = mutatePop(x, mu)  %变异函数 返回一个单个个体 此处采用单点变异 mu为变异概率
if(rand <= mu)
n = numel(x);
s = randi(n);
if x(s) == 0
    x(s) = 1;
elseif x(s) == 1
    x(s) = 0;
end
end
p = x;  %结果输出
end

% %% 变异函数 选用三个变异算子 分阶段进行变异 能够同时兼顾多个方面
% function p = mutatePop(x, it, x_best, nVat, a, b)  %变异函数  mu为变异概率 it为迭代次数 nVat为变量个数
% if(mod(it, 3) == 0)
%     sd = sqrt(abs(x_best-x)/6);    %sd为正态分布的标准差 
%     p = x + normrnd(x, sd, 1, nVat); 
% elseif(mod(it, 3) == 1)
%     p = x + trnd(1, 1, nVat);
% elseif(mod(it, 3) == 2)
%     Rand = -1 + 2*rand(1, nVat);
%     p = x + Rand.*(b-a)/it;
% end
% end