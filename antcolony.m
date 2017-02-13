clear all;
close all;
clc
tic
tw=[10 30 70 95%%%%时间窗
    15 35 50 80
    12 32 60 96
    5  25 55 70
    17 36 40 62
    8  50 55 70
    30 55 80 93
    22 43 66 88
    20 50 80 99
    18 38 45 78];%%%%%每一行代表一个卫星的时间窗口，两两表示一个时间窗
sumtw=20;%%%%%%%%%%时间窗的总数量
twhang=2;%%%%%每个卫星具有的时间窗
sat=10;%%%%%卫星数量

taskt=[11 9 6%%%%%%不同的任务在矩阵中位置为任务编号，数值为任务持续时间
       13 8 7
       10 4 5
       7 6 10
       11 6 4
       8 5  8
       10 8 5
       9 9 6
       7 12 9
       10 5 8];
taskp=[4 6 6%%%%%任务权值，与任务时间矩阵一一对应
       13 5 8
       2 14 8
       5 9 10
       11 6 4
       7 7 7
       9 8 5
       6 15 9
       5 12 7
       10 3 8]; 
sumtask=30;%%%%%总的任务数量
taskhang=3;%%%%%任务表每行的任务数

N=50;%%%%%%%%%迭代次数
anssnum=100;%%%%%%%%蚂蚁的数量
a=1.5;
b=1.5;
c=1;
p1=0.1;%%%%%%%%局部信息素更新，所有蚂蚁中最优解加p1
p2=0.5;%%%%%%%%全局信息素更新，所有迭代中最优解加p2
pmax=8;%Q/(1-p)*Popt/100;%%%信息素最大值
pmin=1;%TTmax/800;%%%%%信息素最小值
p=0.95;%%%%%%%%信息素挥发系数
T0=2;%%Q/500/(1-p);%%%%%%%%%初始信息素量
T=T0*ones(sat,taskhang);%%%%%%%%%初始化信息素矩阵
Popt_show=zeros(1,N);%%%%%%记录N次循环里每次循环得到的最优的权值
taskopt=zeros(1,2*sumtask);%%%%%记录最终最优解的调度任务顺序
twtask=zeros(sat,2*sumtask,1);%%%%%%记录最终最优解的调度任务时间
Popt=0;%%%%%最优的权值

for all=1:1:N%%%%%%%%循环N次(50)
    probb=(T.^a).*(taskp.^b)./(taskt.^c);%%%%%%%转移概率公式
    sump=zeros(1,anssnum);   %%%%%%%用于记录每个蚂蚁的总任务权值
    taskji=zeros(anssnum,2*sumtask);%%%%%%记录每只蚂蚁行走路线
    tasktime=zeros(sat,2*sumtask,anssnum);%%%%%记录每只蚂蚁任务时间表
    for anss=1:1:anssnum
        tabu=ones(sat,taskhang);%%%%%禁忌表:第几个元素为0，就说明第几个走过了，1代表可行
        time=5;   %%%%%%当前时刻，出发点为最早的时间窗开始时间
      
        for taskk=1:1:sumtask%%%%对一只蚂蚁进行任务调度循环
           while(time<100) 
               allow=zeros(sat,taskhang);%%%%%记录可行解，1代表可行，用于运算的中间变量
               for i=1:1:sat      %%%%用于找到可行解范围
                    for j=1:1:twhang
                        if time>=tw(i,2*j-1)&&time<tw(i,2*j)%%%%任务开始时间要在时间窗口内
                            for jj=1:1:taskhang
                                if taskt(i,jj)+time<=tw(i,2*j)%%%%任务结束时间要在时间窗口内
                                    allow(i,jj)=tabu(i,jj);%%%%可行解更新用于运算转移概率，中间变量
                                end
                            end
                            break;
                        end
                    end
                end
                if(sum(sum(allow))==0)%%%%allow中所有元素的和为0，该时刻没有可行的任务
                    time=time+1;
                else 
                    break;
                end
            end
            if(time==100)
                break;
            end 
            probbb=probb.*allow;%%%%%%%转移概率矩阵更新，可行的转移路线
            prob=probbb./(sum(sum(probbb))+eps);%%%%%%%计算转移概率
            q0=rand;%%%%%%%%%%%%根据城市转移概率得到下一个任务（伪随机规则）
            if(q0<=0.5)%%%%%%按照最大的转移概率路线前进
                [Q,W]=max(prob);%%%%%%%%%%%%得到转移概率矩阵各列最大值及对应在所在行
                [TT,R]=max(max(prob));%%%%%%%%%%%%得到转移概率矩阵最大值及对应所在列
                taskji(anss,2*taskk-1)=W(R);%%%%%%%%%%%%得到转移概率矩阵最大值及对应所在行并赋给taskji
                taskji(anss,2*taskk)=R;
                tabu(W(R),R)=0;%%%%%%%%%%%%更新禁忌表
                tasktime(taskji(anss,2*taskk-1),2*taskk-1,anss)=time;
                tasktime(taskji(anss,2*taskk-1),2*taskk,anss)=time+taskt(W(R),R);
           else%%%%%%不按照最大的转移概率路线前进，重新决定前进路线
                r=rand;
                P1=0;P2=0;
                for m=1:1:sat
                    for n=1:1:taskhang
                        P1=P2;
                        P2=P2+prob(m,n);
                        if r>=P1&&r<P2%%%%%%%%利用随机数r比较每一个任务的转移概率，直到满足条件，选择该任务为下一个路径
                            taskji(anss,2*taskk-1)=m;
                            taskji(anss,2*taskk)=n;
                            tabu(m,n)=0;%%%%%%%%%%%%更新禁忌表
                            tasktime(taskji(anss,2*taskk-1),2*taskk-1,anss)=time;
                            tasktime(taskji(anss,2*taskk-1),2*taskk,anss)=time+taskt(m,n);
                        end
                    end
                end%%%%%%%%%%%%根据城市转移概率得到下一个任务
           end
           if taskji(anss,2*taskk-1)~=0%%%%%该蚂蚁经过了此条路径即（2*taskk-1，2*taskk）
                time=time+taskt(taskji(anss,2*taskk-1),taskji(anss,2*taskk));%%%%%更新下一个任务开始时间（当前时刻）
                sump(anss)=sump(anss)+taskp(taskji(anss,2*taskk-1),taskji(anss,2*taskk));   %%%%%更新该蚂蚁目前的任务权值
           end
        end%%%%一只蚂蚁寻找路径的循环完毕
    end%%%%%%%%%%%%Num只蚂蚁循环结束 
    [C,I]=max(sump);
    T=T.*p;%%%%%%%%%%%%%%%%%%%%%%信息素量更新
    for nn=1:1:sumtask
       if taskji(I,2*nn-1)~=0%%%%%该蚂蚁经过了此条路径即（2*nn-1，2*nn）
           T(taskji(I,2*nn-1),taskji(I,2*nn))=T(taskji(I,2*nn-1),taskji(I,2*nn))+C*p1;%%%更新该点的信息素（局部），?sump(I)/50?;
       end
    end
    if C>Popt%%%%%%%%%%%%%%更新最优路径
        Popt=C;
        taskopt(1:2*sumtask)=taskji(I,:);%taskopt中从1到2*sumtask的所有元素赋值taskji第I行的所有元素
        twtask(:,:,1)=tasktime(:,:,I);%%%%%将最优路径蚂蚁的路线记录在最终最优路径时间表
    end%%%%%%%%%%%%%%更新最优路径
     for nn=1:1:sumtask
          if taskopt(2*nn-1)~=0
               T(taskopt(2*nn-1),taskopt(2*nn))=T(taskopt(2*nn-1),taskopt(2*nn))+Popt*p2;%记录该点为最优路径更新该点的信息素（全局）,sump(I)/50;
          end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%信息素范围控制
   for iii=1:sat
       for jjj=1:taskhang
           if T(iii,jjj)>pmax
              T(iii,jjj)=pmax;
           elseif T(iii,jjj)<pmin
               T(iii,jjj)=pmin;
           end
       end
   end
  %%%%%%%%%%%%%%%%%%%%%%%%%%信息素范围控制
    Popt_show(all)=C;%记录本次迭代的最优权值
end%%%%%%%%%%%%%%%%%N此循环结束

figure(1);
for i=1:1:sat;
    for j=1:1:twhang;
plot([tw(i,2*j-1),tw(i,2*j)],[i,i],'r-','LineWidth',6);
hold on;
    end
    for e=1:1:sumtask
    if(twtask(i,2*e-1,1)~=0)
plot([twtask(i,2*e-1,1),twtask(i,2*e,1)],[i,i],'LineWidth',2);
hold on;
    end
    end
end
xlabel('Time/minute ');
ylabel('Identification Number of Satelite');
axis([0 100 0 12]);

figure(2);
i=1:1:N;
plot(i,Popt_show(i),'k.');
fprintf('最佳总优先级为%f,',Popt);
fprintf('收益比为%d\n',Popt/sum(sum(taskp)));
for e=1:1:sumtask
    if(taskopt(2*e-1)~=0)
        fprintf('第%d调度卫星上的',taskopt(2*e-1));
        fprintf(' 第%d个任务',taskopt(2*e));
        fprintf(' 优先级为%d\n',taskp(taskopt(2*e-1),taskopt(2*e)));
    end
end
xlabel('The number of Iterations');
ylabel('Total Priority');
axis([0 50 0 125]);
toc

    
