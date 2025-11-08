clc,clear;
%% 初始化种群 每个个体相当于一组解
nVat = 2; %自变量个数 也可以看成自变量染色体上的基因数 一个就是一个基因
nPop = 30; %种群规模大小
maxIt = 2000; %最大进化次数
nPc = 0.8;  %相当于繁育下一代的比例 交叉率 每次俩个个体
nMu = 0.01;
% nC = nPop*nPc; %父代乘以繁育率得到下一代的规模大小 
nC = round(nPop*nPc/2)*2;  %技术处理 确保交叉生成的子代规模为偶数

template.x = [];   %定义空的结构体 绑定x和y的关系
template.y = [];

Parent = repmat(template, nPop, 1);   %将结构体体拷贝 生成30*1的父代群体 每一个个体有x y俩个变量
% Offspring = repmat(template, nC, 1);   %将结构体体拷贝 生成30*1的父代群体 每一个个体有x y俩个变量

%% 遗传算法
%先构建初始父代种群
for i = 1 : nPop
%     Parent(i).x = randi([0 1],1,nVat);   %每一个个体的自变量是一个 100*1的行向量 每一个xi 均是0或者1
    Parent(i).x = 6*rand(1, nVat);
    Parent(i).y = goal(Parent(i).x);     %求解每一个y
end
% 开始遗传选择部分
for It = 1 : maxIt
    
    Offspring = repmat(template, nC/2, 2);   %每次It迭代 刷新方便重新生成 写成俩列方便查看交叉结果
    
    for j = 1 : nC/2                %交叉操作 父代交叉生成子代
        p1 =  SelectPort(Parent);   %随机选出俩个 P1 P2优胜个体
        p2 =  SelectPort(Parent);
        %以下先对俩个优胜个体x进行交叉 然后得到了俩个子代 依次迭代 得到所有的子代
        [Offspring(j, 1).x, Offspring(j, 2).x] = crossPort(p1.x, p2.x);   
    end
    
    Offspring = Offspring(:);  %相当于按顺序列出所有元素 依次赋给 然后变成24*1的结构体
    
    for k = 1 : nC           %变异操作 生成的子代个体依次发生变异
        Offspring(k).x =  mutatePop(Offspring(k).x, nMu);
        Offspring(k).y =  goal( Offspring(k).x);
    end
    
    %淘汰并且生成新的父代
    newPop = [Parent; Offspring];
    [~,so] = sort([newPop.y], 'ascend');   %so 相当于按照升序排列储存的索引下标数组 好用
    newPop = newPop(so);   %赋值后 相当于newPOP中的元素 按照指定大小顺序排列
    Parent = newPop(1 : nPop);  %只取最小的三十 其余的去掉 相当于进行了筛选 生成新的父代
    disp(['迭代次数：', num2str(It), ', 最优解:', num2str(Parent(1).y)]); 
end

% 
%% 目标函数
function out = goal(input)

out = (input(1)*input(1)+input(2)-11)^2 + (input(1)+input(2)^2-7)^2;
end