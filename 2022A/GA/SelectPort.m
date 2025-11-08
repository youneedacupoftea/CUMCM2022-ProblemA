
function p = SelectPort(Parent)  %选择函数 输入为父代种群 这里采用锦标赛选择法 返回一个个体 优胜个体
n = numel(Parent);       %先计算父代种群的大小
index = randperm(n);     %将三十个下标打乱成随机序列 每次运行结果不同 可以认为前俩个即为选中的个体
p1 = Parent(index(1));
p2 = Parent(index(2));
if p1.y <= p2.y   %选取个小的P出来
    p = p1;
else
    p = p2;
end
end

% function p = SelectPort(Parent) %输入种群个体数 和父代信息 然后进行排序选择 返回为排序的组合
% [~, p] = sort([Parent.y], 'descend');   %对输入父代进行降序排列 即从大到小 P储存索引信息
% end

