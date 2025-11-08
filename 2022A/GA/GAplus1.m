clc,clear;
%% 初始化种群
nVat = 2; %自变量个数 也可以看成自变量染色体上的基因数 一个就是一个基因
nPop = 160; %种群规模大小 不能小
maxIt = 4000; %最大进化次数
s = 80;   %精英个体的数目 不能太小
nPc = 1.0;  %相当于繁育下一代的比例 交叉率 每次俩个个体
nMu = 0.5; %变异概率 不建议太大 最好不要过0.5
nC = round(nPop*nPc/2)*2;  %技术处理 确保交叉生成的子代规模为偶数
a = -10*ones(1,2);   %自变量下界  %% 需要更改！！！！！！！！！！！
b = 10*ones(1,2);  %自变量上界    %% 需要更改！！！！！！！！！！
M1 = 1e+7;
M2 = 1e+6;

template.x = [];   %定义空的结构体 绑定x和y的关系
template.y = [];   %y表示temp 
template.z = [];   %z表示实际的所求

Parent = repmat(template, nPop, 1);   %将结构体体拷贝 生成30*1的父代群体 每一个个体有x y俩个变量
Best = repmat(template, 1, 1);  %种群中的最优个体

%% 遗传算法
%先构建初始父代种群
for i = 1 : nPop
    % 需要更改  ！！！！！！！！！
    Parent(i).x = -10 + 20*rand(1, 2);  %每一个个体的自变量是一个 1*2 的行向量 每一个xi 在0到6之间
    [Parent(i).z, Parent(i).y] = goal(Parent(i).x(1), Parent(i).x(2), M2);     %求解每一个y
end


%% 开始遗传选择部分
for It = 1 : maxIt

    Offspring = repmat(template, nC/2, 2);   %每次It迭代 刷新方便重新生成 写成俩列方便查看交叉结果
    Best_pa = repmat(template, s, 1);       %生成父代精英个体
    Best_son = repmat(template, s, 1);       %生成子代精英个体
    
    
    [~, so0] = sort([Parent.y], 'descend');  %首先进行排序 方便选择和交叉
    Parent = Parent(so0);      %根据排序重新赋值
    Best_pa = Parent(1 : s);  %选出五十个精英个体
    
    % 种群交叉
    [Offspring] = crossplus(Parent, nPop, nVat);  %交叉操作
%    
%     % disp(['交叉结束：', num2str(It), ', 最优解:', num2str(Parent(1).x)]); 
    
    %种群排序 方便进行替换
    Offspring = Offspring(:);    %结构体结构变换
    for i = 1 : nC
    [Offspring(i).z, Offspring(i).y] = goal(Offspring(i).x(1), Offspring(i).x(2), M2);     %求解每一个y
    end
    [~, so1] = sort([Offspring .y], 'descend');
    Offspring =  Offspring(so1);       %根据排序结果重组
    Best_son = Offspring(1 : s);      %选出五十个精英个体

    
    %种群替代 增加多样性
    for j = 1 : nC-1           %根据重组的结果 简单判断并且进行替代操作
        i = j + 1;
        if(Offspring(j).x == Offspring(i).x)
            % 需要更改 ！！！！！！！！！！！
            Offspring(j).x = -10 + 20*rand(1, 2);
            [Offspring(j).z, Offspring(j).y] = goal(Offspring(j).x(1), Offspring(j).x(2), M2);
        end
    end
%     
%     
    %种群变异
    [~, so2] = sort([Offspring .y], 'descend');
    Offspring =  Offspring(so2);       %根据排序结果重组
    Best = Offspring(1);         %选出最优个体
    for j = 1 : nC
        if rand < nMu
        Offspring(j).x = mutateplus(Offspring(j).x, It, Best.x, nVat, a, b);
        [Offspring(j).z, Offspring(j).y] = goal(Offspring(j).x(1), Offspring(j).x(2), M2);
        end
    end
%     
%     
    %种群精英淘汰
    [~, so3] = sort([Offspring .y], 'descend');
    Offspring =  Offspring(so3);       %根据排序结果重组
    Offspring(nC - s + 1 : end) =  Best_son;  %替换成精英种群
    Parent = Offspring;
%     if(abs((Parent(1).z)-13.59)<0.001)
%         break;
%     end
end
 disp(['迭代次数：', num2str(It), ', 最优解:', num2str(Parent(1).y)]); 
 
 
%% 目标函数
function[out, temp] = goal(input1, input2, M2)
temp = 0.5 - (sin(sqrt(input1^2+input1^2))^2-0.5)/(1+0.001*(input1^2+input2^2))^2;
out = 0;
end
% function [out,temp] = goal(input1, input2, M2)
% h1 = (input1-0.05)^2 + (input2-2.5)^2-4.84;
% h2 = 4.84 - input1^2 -(input2-2.5)^2;
% h3 = input1;
% h4 = input2;
% h5 = 6-input1;
% h6=  6-input2;
% h1 = -h1;
% h2 = -h2;
% % out = (input1*input1+input2-11)^2 + (input1+input2^2-7)^2;
% out = (input1*input1+input2-11)^2 + (input1+input2^2-7)^2 + ...
% M2*((min(h1,0))^2+(min(h2,0))^2+(min(h3,0))^2+(min(h4,0))^2+...
% (min(h5,0))^2+(min(h6,0))^2); %实际求out最小值 转换成求temp最大值
% % fprintf("1 \n");
% % min(h1,0)
% % min(h2,0) 
% % min(h3,0) 
% % min(h4,0)
% % min(h5,0) 
% % min(h6,0)
% % fprintf("6 \n");
% temp = - out;
% end