% %% 选择函数 改进版 放弃轮盘赌和竞标赛的方法 选用排序分组的方法
% function p = selectplus(Parent, nC) %输入种群个体数 和父代信息 然后进行排序选择 返回为排序的组合
% for i = 1 : nC
%     [Parent(i).z, Parent(i).y] = goal(Parent(i).x(1), Parent(i).x(2));    
% end
% [~, p] = sort([Parent.y], 'descend');   %对输入父代进行降序排列 即从大到小 P储存索引信息 记得paret.y加中括号
% end