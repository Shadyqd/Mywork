clear all;
close all;
clc
tic
tw=[10 30 70 95%%%%ʱ�䴰
    15 35 50 80
    12 32 60 96
    5  25 55 70
    17 36 40 62
    8  50 55 70
    30 55 80 93
    22 43 66 88
    20 50 80 99
    18 38 45 78];%%%%%ÿһ�д���һ�����ǵ�ʱ�䴰�ڣ�������ʾһ��ʱ�䴰
sumtw=20;%%%%%%%%%%ʱ�䴰��������
twhang=2;%%%%%ÿ�����Ǿ��е�ʱ�䴰
sat=10;%%%%%��������

taskt=[11 9 6%%%%%%��ͬ�������ھ�����λ��Ϊ�����ţ���ֵΪ�������ʱ��
       13 8 7
       10 4 5
       7 6 10
       11 6 4
       8 5  8
       10 8 5
       9 9 6
       7 12 9
       10 5 8];
taskp=[4 6 6%%%%%����Ȩֵ��������ʱ�����һһ��Ӧ
       13 5 8
       2 14 8
       5 9 10
       11 6 4
       7 7 7
       9 8 5
       6 15 9
       5 12 7
       10 3 8]; 
sumtask=30;%%%%%�ܵ���������
taskhang=3;%%%%%�����ÿ�е�������

N=50;%%%%%%%%%��������
anssnum=100;%%%%%%%%���ϵ�����
a=1.5;
b=1.5;
c=1;
p1=0.1;%%%%%%%%�ֲ���Ϣ�ظ��£��������������Ž��p1
p2=0.5;%%%%%%%%ȫ����Ϣ�ظ��£����е��������Ž��p2
pmax=8;%Q/(1-p)*Popt/100;%%%��Ϣ�����ֵ
pmin=1;%TTmax/800;%%%%%��Ϣ����Сֵ
p=0.95;%%%%%%%%��Ϣ�ػӷ�ϵ��
T0=2;%%Q/500/(1-p);%%%%%%%%%��ʼ��Ϣ����
T=T0*ones(sat,taskhang);%%%%%%%%%��ʼ����Ϣ�ؾ���
Popt_show=zeros(1,N);%%%%%%��¼N��ѭ����ÿ��ѭ���õ������ŵ�Ȩֵ
taskopt=zeros(1,2*sumtask);%%%%%��¼�������Ž�ĵ�������˳��
twtask=zeros(sat,2*sumtask,1);%%%%%%��¼�������Ž�ĵ�������ʱ��
Popt=0;%%%%%���ŵ�Ȩֵ

for all=1:1:N%%%%%%%%ѭ��N��(50)
    probb=(T.^a).*(taskp.^b)./(taskt.^c);%%%%%%%ת�Ƹ��ʹ�ʽ
    sump=zeros(1,anssnum);   %%%%%%%���ڼ�¼ÿ�����ϵ�������Ȩֵ
    taskji=zeros(anssnum,2*sumtask);%%%%%%��¼ÿֻ��������·��
    tasktime=zeros(sat,2*sumtask,anssnum);%%%%%��¼ÿֻ��������ʱ���
    for anss=1:1:anssnum
        tabu=ones(sat,taskhang);%%%%%���ɱ�:�ڼ���Ԫ��Ϊ0����˵���ڼ����߹��ˣ�1�������
        time=5;   %%%%%%��ǰʱ�̣�������Ϊ�����ʱ�䴰��ʼʱ��
      
        for taskk=1:1:sumtask%%%%��һֻ���Ͻ����������ѭ��
           while(time<100) 
               allow=zeros(sat,taskhang);%%%%%��¼���н⣬1������У�����������м����
               for i=1:1:sat      %%%%�����ҵ����нⷶΧ
                    for j=1:1:twhang
                        if time>=tw(i,2*j-1)&&time<tw(i,2*j)%%%%����ʼʱ��Ҫ��ʱ�䴰����
                            for jj=1:1:taskhang
                                if taskt(i,jj)+time<=tw(i,2*j)%%%%�������ʱ��Ҫ��ʱ�䴰����
                                    allow(i,jj)=tabu(i,jj);%%%%���н������������ת�Ƹ��ʣ��м����
                                end
                            end
                            break;
                        end
                    end
                end
                if(sum(sum(allow))==0)%%%%allow������Ԫ�صĺ�Ϊ0����ʱ��û�п��е�����
                    time=time+1;
                else 
                    break;
                end
            end
            if(time==100)
                break;
            end 
            probbb=probb.*allow;%%%%%%%ת�Ƹ��ʾ�����£����е�ת��·��
            prob=probbb./(sum(sum(probbb))+eps);%%%%%%%����ת�Ƹ���
            q0=rand;%%%%%%%%%%%%���ݳ���ת�Ƹ��ʵõ���һ������α�������
            if(q0<=0.5)%%%%%%��������ת�Ƹ���·��ǰ��
                [Q,W]=max(prob);%%%%%%%%%%%%�õ�ת�Ƹ��ʾ���������ֵ����Ӧ��������
                [TT,R]=max(max(prob));%%%%%%%%%%%%�õ�ת�Ƹ��ʾ������ֵ����Ӧ������
                taskji(anss,2*taskk-1)=W(R);%%%%%%%%%%%%�õ�ת�Ƹ��ʾ������ֵ����Ӧ�����в�����taskji
                taskji(anss,2*taskk)=R;
                tabu(W(R),R)=0;%%%%%%%%%%%%���½��ɱ�
                tasktime(taskji(anss,2*taskk-1),2*taskk-1,anss)=time;
                tasktime(taskji(anss,2*taskk-1),2*taskk,anss)=time+taskt(W(R),R);
           else%%%%%%����������ת�Ƹ���·��ǰ�������¾���ǰ��·��
                r=rand;
                P1=0;P2=0;
                for m=1:1:sat
                    for n=1:1:taskhang
                        P1=P2;
                        P2=P2+prob(m,n);
                        if r>=P1&&r<P2%%%%%%%%���������r�Ƚ�ÿһ�������ת�Ƹ��ʣ�ֱ������������ѡ�������Ϊ��һ��·��
                            taskji(anss,2*taskk-1)=m;
                            taskji(anss,2*taskk)=n;
                            tabu(m,n)=0;%%%%%%%%%%%%���½��ɱ�
                            tasktime(taskji(anss,2*taskk-1),2*taskk-1,anss)=time;
                            tasktime(taskji(anss,2*taskk-1),2*taskk,anss)=time+taskt(m,n);
                        end
                    end
                end%%%%%%%%%%%%���ݳ���ת�Ƹ��ʵõ���һ������
           end
           if taskji(anss,2*taskk-1)~=0%%%%%�����Ͼ����˴���·������2*taskk-1��2*taskk��
                time=time+taskt(taskji(anss,2*taskk-1),taskji(anss,2*taskk));%%%%%������һ������ʼʱ�䣨��ǰʱ�̣�
                sump(anss)=sump(anss)+taskp(taskji(anss,2*taskk-1),taskji(anss,2*taskk));   %%%%%���¸�����Ŀǰ������Ȩֵ
           end
        end%%%%һֻ����Ѱ��·����ѭ�����
    end%%%%%%%%%%%%Numֻ����ѭ������ 
    [C,I]=max(sump);
    T=T.*p;%%%%%%%%%%%%%%%%%%%%%%��Ϣ��������
    for nn=1:1:sumtask
       if taskji(I,2*nn-1)~=0%%%%%�����Ͼ����˴���·������2*nn-1��2*nn��
           T(taskji(I,2*nn-1),taskji(I,2*nn))=T(taskji(I,2*nn-1),taskji(I,2*nn))+C*p1;%%%���¸õ����Ϣ�أ��ֲ�����?sump(I)/50?;
       end
    end
    if C>Popt%%%%%%%%%%%%%%��������·��
        Popt=C;
        taskopt(1:2*sumtask)=taskji(I,:);%taskopt�д�1��2*sumtask������Ԫ�ظ�ֵtaskji��I�е�����Ԫ��
        twtask(:,:,1)=tasktime(:,:,I);%%%%%������·�����ϵ�·�߼�¼����������·��ʱ���
    end%%%%%%%%%%%%%%��������·��
     for nn=1:1:sumtask
          if taskopt(2*nn-1)~=0
               T(taskopt(2*nn-1),taskopt(2*nn))=T(taskopt(2*nn-1),taskopt(2*nn))+Popt*p2;%��¼�õ�Ϊ����·�����¸õ����Ϣ�أ�ȫ�֣�,sump(I)/50;
          end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%��Ϣ�ط�Χ����
   for iii=1:sat
       for jjj=1:taskhang
           if T(iii,jjj)>pmax
              T(iii,jjj)=pmax;
           elseif T(iii,jjj)<pmin
               T(iii,jjj)=pmin;
           end
       end
   end
  %%%%%%%%%%%%%%%%%%%%%%%%%%��Ϣ�ط�Χ����
    Popt_show(all)=C;%��¼���ε���������Ȩֵ
end%%%%%%%%%%%%%%%%%N��ѭ������

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
fprintf('��������ȼ�Ϊ%f,',Popt);
fprintf('�����Ϊ%d\n',Popt/sum(sum(taskp)));
for e=1:1:sumtask
    if(taskopt(2*e-1)~=0)
        fprintf('��%d���������ϵ�',taskopt(2*e-1));
        fprintf(' ��%d������',taskopt(2*e));
        fprintf(' ���ȼ�Ϊ%d\n',taskp(taskopt(2*e-1),taskopt(2*e)));
    end
end
xlabel('The number of Iterations');
ylabel('Total Priority');
axis([0 50 0 125]);
toc

    
