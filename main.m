clc,clear,close all
warning off
format longG
% BFO优化算法 参数
sizepop = 20;              % 种群数量
Nc = 50;                   % 趋化次数
Ns = 4;                    % 游动次数
C(:,1) = 0.001*ones(sizepop,1);            % 翻转选定方向后，单个细菌前进的步长
Nre = 4;                   % 复制次数
Ned = 2;                   % 驱散(迁移)次数
Sr = ceil( sizepop/2 );    % 复制（分裂）次数
Ped = 0.25;                             % 细菌驱散(迁移)概率
d_attract = 0.05;          % 吸引剂的数量
w_attract = 0.05;          % 吸引剂的释放速度
h_repellant = 0.05;        % 排斥剂的数量
w_repellant = 0.05;        % 排斥剂的释放速度
ww = 0.5;                  % 适应度增量因子

nvar = 2;   % 2个未知量
popmin1 = -1;  popmax1 = 1; % x1
popmin2 = -1;  popmax2 = 1; % x2
Cmin = -1;  % 最小步长
Cmax = 1;   % 最大步长
%% 初始化种群
for i=1:sizepop
    x1 = popmin1 + (popmax1-popmin1)*rand;
    x2 = popmin2 + (popmax2-popmin2)*rand;
    pop(i,1) = x1;
    pop(i,2) = x2;                     % 初始化个体
    fitness(i) = fun([x1,x2]);         % 适应度值
    C(i,1) = Cmin + (Cmax-Cmin)*rand;  % 步长
end
clear x1 x2
%% 记录一组最优值
[bestfitness,bestindex]=min(fitness);
zbest=pop(bestindex,:);   % 全局最佳
fitnesszbest=bestfitness; % 全局最佳适应度值
NcSizepop = 0;            % 记录最优适应度值（函数值）
%% 迭代寻优
for i = 1:Ned                   % 驱散(迁移)次数
    for k = 1:Nre               % 复制次数
        
        for m = 1:Nc            % 趋化次数
            for j=1:sizepop     % 种群
                % Jcc计算
                Jcc = sum( -d_attract*exp( w_attract*sum(pop(j,1)-pop(:,1).^2) ) + ...
                    h_repellant*exp( w_repellant*sum(pop(j,2)-pop(:,2).^2) ));
               
                poplast = pop(j,:);  % 当前的种群个体
                fitness(j) = fitness(j) + ww*Jcc;
                Jlast = fitness(j);  % 当前适应度值
      
                % 翻转
                delta = 2*rand(1,nvar)-0.5;
                pop(j,:) = pop(j,:) + C(j,:).*delta./(sqrt( delta*delta' ));
               
                % x1
                if pop(j,1)>popmax1
                    pop(j,1)=popmax1;
                end
                if pop(j,1)<popmin1
                    pop(j,1)=popmin1;
                end
                % x2
                if pop(j,2)>popmax2
                    pop(j,2)=popmax2;
                end
                if pop(j,2)<popmin2
                    pop(j,2)=popmin2;
                end
               
                % 更新当前适应度值？
                fitness(j) = fun(pop(j,:));
                % 游动
                n=0;
                while(n<Ns)     % 游动次数
                    if fitness(j)<Jlast
                        Jlast = fitness(j);
                        poplast = pop(j,:);
                    else       % 不更新适应度值
                        n=Ns;
                    end
                end

                % 适应度更新
                % 比较 个体间比较
                if Jlast<bestfitness
                    bestfitness = Jlast;
%                     zbest =  pop(j,:);
                    zbest = poplast;
                end
            end   % sizepop  种群数量
            
            % 记录最优适应度值
            NcSizepop = NcSizepop+1;
            fitness_iter(NcSizepop) = bestfitness;
        end       % Nc       趋化次数
        
        % 复制操作
        [maxF,index] = sort(fitness,'descend');  % 降序排列
        for Nre2 = 1:Sr   % 将最大适应度值的Sr个种群，进行更新
            pop(index(Nre2),1) = popmin1 + (popmax1-popmin1)*rand;
            pop(index(Nre2),2) = popmin2 + (popmax2-popmin2)*rand;
            fitness(index(Nre2)) = fun(pop(index(Nre2),:));
            C(index(Nre2),1) = Cmin + (Cmax-Cmin)*rand;  % 步长
            % 比较 个体间比较
            if fitness(index(Nre2))<bestfitness
                bestfitness = fitness(index(Nre2));
                zbest =  pop(index(Nre2),:);
            end
        end
    end   % Nre  复制操作
   
    for j=1:sizepop     % 种群
        if Ped>rand
            pop(j,1) = popmin1 + (popmax1-popmin1)*rand;
            pop(j,2) = popmin2 + (popmax2-popmin2)*rand;
            fitness(j) = fun(pop(j,:));
            % 比较 个体间比较
            if fitness(j)<bestfitness
                bestfitness = fitness(j);
                zbest =  pop(j,:);
            end
        end
    end
   
end       % Ned   驱散(迁移)次数

disp('最优解')
disp(zbest)
fprintf('\n')

figure('color',[1,1,1])
plot(fitness_iter,'ro-','linewidth',2)
% loglog(fitness_iter,'ro-','linewidth',2)
axis tight
grid on