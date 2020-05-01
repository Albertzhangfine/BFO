clc,clear,close all
warning off
format longG
% BFO�Ż��㷨 ����
sizepop = 20;              % ��Ⱥ����
Nc = 50;                   % ��������
Ns = 4;                    % �ζ�����
C(:,1) = 0.001*ones(sizepop,1);            % ��תѡ������󣬵���ϸ��ǰ���Ĳ���
Nre = 4;                   % ���ƴ���
Ned = 2;                   % ��ɢ(Ǩ��)����
Sr = ceil( sizepop/2 );    % ���ƣ����ѣ�����
Ped = 0.25;                             % ϸ����ɢ(Ǩ��)����
d_attract = 0.05;          % ������������
w_attract = 0.05;          % ���������ͷ��ٶ�
h_repellant = 0.05;        % �ų��������
w_repellant = 0.05;        % �ų�����ͷ��ٶ�
ww = 0.5;                  % ��Ӧ����������

nvar = 2;   % 2��δ֪��
popmin1 = -1;  popmax1 = 1; % x1
popmin2 = -1;  popmax2 = 1; % x2
Cmin = -1;  % ��С����
Cmax = 1;   % ��󲽳�
%% ��ʼ����Ⱥ
for i=1:sizepop
    x1 = popmin1 + (popmax1-popmin1)*rand;
    x2 = popmin2 + (popmax2-popmin2)*rand;
    pop(i,1) = x1;
    pop(i,2) = x2;                     % ��ʼ������
    fitness(i) = fun([x1,x2]);         % ��Ӧ��ֵ
    C(i,1) = Cmin + (Cmax-Cmin)*rand;  % ����
end
clear x1 x2
%% ��¼һ������ֵ
[bestfitness,bestindex]=min(fitness);
zbest=pop(bestindex,:);   % ȫ�����
fitnesszbest=bestfitness; % ȫ�������Ӧ��ֵ
NcSizepop = 0;            % ��¼������Ӧ��ֵ������ֵ��
%% ����Ѱ��
for i = 1:Ned                   % ��ɢ(Ǩ��)����
    for k = 1:Nre               % ���ƴ���
        
        for m = 1:Nc            % ��������
            for j=1:sizepop     % ��Ⱥ
                % Jcc����
                Jcc = sum( -d_attract*exp( w_attract*sum(pop(j,1)-pop(:,1).^2) ) + ...
                    h_repellant*exp( w_repellant*sum(pop(j,2)-pop(:,2).^2) ));
               
                poplast = pop(j,:);  % ��ǰ����Ⱥ����
                fitness(j) = fitness(j) + ww*Jcc;
                Jlast = fitness(j);  % ��ǰ��Ӧ��ֵ
      
                % ��ת
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
               
                % ���µ�ǰ��Ӧ��ֵ��
                fitness(j) = fun(pop(j,:));
                % �ζ�
                n=0;
                while(n<Ns)     % �ζ�����
                    if fitness(j)<Jlast
                        Jlast = fitness(j);
                        poplast = pop(j,:);
                    else       % ��������Ӧ��ֵ
                        n=Ns;
                    end
                end

                % ��Ӧ�ȸ���
                % �Ƚ� �����Ƚ�
                if Jlast<bestfitness
                    bestfitness = Jlast;
%                     zbest =  pop(j,:);
                    zbest = poplast;
                end
            end   % sizepop  ��Ⱥ����
            
            % ��¼������Ӧ��ֵ
            NcSizepop = NcSizepop+1;
            fitness_iter(NcSizepop) = bestfitness;
        end       % Nc       ��������
        
        % ���Ʋ���
        [maxF,index] = sort(fitness,'descend');  % ��������
        for Nre2 = 1:Sr   % �������Ӧ��ֵ��Sr����Ⱥ�����и���
            pop(index(Nre2),1) = popmin1 + (popmax1-popmin1)*rand;
            pop(index(Nre2),2) = popmin2 + (popmax2-popmin2)*rand;
            fitness(index(Nre2)) = fun(pop(index(Nre2),:));
            C(index(Nre2),1) = Cmin + (Cmax-Cmin)*rand;  % ����
            % �Ƚ� �����Ƚ�
            if fitness(index(Nre2))<bestfitness
                bestfitness = fitness(index(Nre2));
                zbest =  pop(index(Nre2),:);
            end
        end
    end   % Nre  ���Ʋ���
   
    for j=1:sizepop     % ��Ⱥ
        if Ped>rand
            pop(j,1) = popmin1 + (popmax1-popmin1)*rand;
            pop(j,2) = popmin2 + (popmax2-popmin2)*rand;
            fitness(j) = fun(pop(j,:));
            % �Ƚ� �����Ƚ�
            if fitness(j)<bestfitness
                bestfitness = fitness(j);
                zbest =  pop(j,:);
            end
        end
    end
   
end       % Ned   ��ɢ(Ǩ��)����

disp('���Ž�')
disp(zbest)
fprintf('\n')

figure('color',[1,1,1])
plot(fitness_iter,'ro-','linewidth',2)
% loglog(fitness_iter,'ro-','linewidth',2)
axis tight
grid on