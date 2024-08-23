%_______________________________________________________________________________________%
%  Sea-Horse optimizer (SHO)
%  Developed in MATLAB R2018a
%
%  programmer: Shijie Zhao and Tianran Zhang   
%  E-mail: zhaoshijie@lntu.edu.cn
%          ztr20010118@126.com
%  The code is based on the following papers.
%  Shijie Zhao, Tianran Zhang, Shilin Ma, Mengchen Wang 
%  Sea-horse optimizer: a novel nature-inspired meta-heuristic for global
%  optimization problems.
%  Applied Intelligence
%_______________________________________________________________________________________%
clear all 
clc

pop=20; % Number of search agents
Max_iter=50; % Maximum iteration
%F_name='F1'; % Name of the test function that can be from F1 to F23 (Table 2,3,4 in the paper)

%[LB,UB,Dim,fobj]=BenchmarkFunctions(F_name);% Load details of the selected benchmark function
LB=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
LB=LB';
        UB=[5,5,5,5,5,5,5,5,1,1,1,1,1,1,20,20];
        UB=UB';
        Dim=16;

tic
%[ObjectiveFitness,ObjectivePosition,Convergence_curve,Trajectories,fitness_history, population_history]=SHO(popsize,Max_iter,LB,UB,Dim);
Sea_horses=initialization_FOPID(pop,Dim,UB,LB);
Sea_horsesFitness = zeros(1,pop);
fitness_history=zeros(pop,Max_iter);
population_history=zeros(pop,Dim,Max_iter);
Convergence_curve=zeros(1,Max_iter);
Trajectories=zeros(pop,Max_iter);
%%Calculate the fitness value ofinitial sea horses
for i=1:pop
    Kp1=Sea_horses(i,1);
Kd1=Sea_horses(i,2);
Ki1=Sea_horses(i,3);
Kt1=Sea_horses(i,4);
Kp2=Sea_horses(i,5);
Kd2=Sea_horses(i,6);
Ki2=Sea_horses(i,7);
Kt2=Sea_horses(i,8);
l1=Sea_horses(i,9);
l2=Sea_horses(i,10);
m1=Sea_horses(i,11);
m2=Sea_horses(i,12);
n1=Sea_horses(i,13);
n2=Sea_horses(i,14);
N1=Sea_horses(i,15);
N2=Sea_horses(i,16);
[tout] = sim('btech_ev_1');
o=tout.e;
    Sea_horsesFitness(1,i)=o;
    fitness_history(i,1)=Sea_horsesFitness(1,i);
    population_history(i,:,1)=Sea_horses(i,:);
    Trajectories(:,1)=Sea_horses(:,1);
   
end
Sea_horsesFitness;
%%Find the best seahorses in the initial population
[~,sorted_indexes]=sort(Sea_horsesFitness);

TargetPosition=Sea_horses(sorted_indexes(1),:);
TargetFitness = Sea_horsesFitness(sorted_indexes(1));
Convergence_curve(1)=TargetFitness;
t=1;
u=0.05;
v=0.05;
l=0.05;
while t<Max_iter+1
    beta=randn(pop,Dim);
    Elite=repmat(TargetPosition,pop,1);
    %The motor behavior of sea horses
    r1=randn(1,pop);
    Step_length=levy_FOPID(pop,Dim,1.5);
    for i=1:pop
        for j=1:Dim
            if r1(i)>0  
                r=rand();
                theta=r*2*pi;
                row=u*exp(theta*v);
                x=row*cos(theta);
                y=row*sin(theta);
                z=row*theta;
                Sea_horses_new1(i,j)=Sea_horses(i,j)+Step_length(i,j)*((Elite(i,j)-Sea_horses(i,j)).*x.*y.*z+Elite(i,j));%Eq.(4)
            else
                Sea_horses_new1(i,j)=Sea_horses(i,j)+rand()*l*beta(i,j)*(Sea_horses(i,j)-beta(i,j)* Elite(i,j));%Eq.(7)
            end
        end
    end
     for i=1:pop
        %Bound the variable
        Tp=Sea_horses_new1(i,:)>UB';
        Tm=Sea_horses_new1(i,:)<LB';
        Sea_horses_new1(i,:)=(Sea_horses_new1(i,:).*(~(Tp+Tm)))+UB'.*Tp+LB'.*Tm;
      
    end
    %The predation behavior of sea horses
    for i=1:pop
        for j=1:Dim
            r1=rand();
            r2(i)=rand();
            alpha=(1-t/Max_iter)^(2*t/Max_iter);
            if r2(i)>=0.1
                Sea_horses_new2(i,j)=alpha*(Elite(i,j)-rand()*Sea_horses_new1(i,j))+(1-alpha)*Elite(i,j);  %Eq.(10)
            else
                Sea_horses_new2(i,j)=(1-alpha)*(Sea_horses_new1(i,j)-rand()*Elite(i,j))+alpha*Sea_horses_new1(i,j);  %Eq.(10)
            end
        end
    end
    for i=1:pop
        %Bound the variable
        Tp=Sea_horses_new2(i,:)>UB';
        Tm=Sea_horses_new2(i,:)<LB';
        Sea_horses_new2(i,:)=(Sea_horses_new2(i,:).*(~(Tp+Tm)))+UB'.*Tp+LB'.*Tm;
        %Calculated all sea horses' fitness values
         Kp1=Sea_horses(i,1);
Kd1=Sea_horses(i,2);
Ki1=Sea_horses(i,3);
Kt1=Sea_horses(i,4);
Kp2=Sea_horses(i,5);
Kd2=Sea_horses(i,6);
Ki2=Sea_horses(i,7);
Kt2=Sea_horses(i,8);
l1=Sea_horses(i,9);
l2=Sea_horses(i,10);
m1=Sea_horses(i,11);
m2=Sea_horses(i,12);
n1=Sea_horses(i,13);
n2=Sea_horses(i,14);
N1=Sea_horses(i,15);
N2=Sea_horses(i,16);
[tout] = sim('btech_ev_1');
         o1=tout.e;
        Sea_horsesFitness1(1,i)=o1;
    end
    [~,index]=sort(Sea_horsesFitness1);
    
    %The reproductive behavior of sea horses
    Sea_horses_father=Sea_horses_new2(index(1:pop/2),:);    %Eq.(12)
    Sea_horses_mother=Sea_horses_new2(index(pop/2+1:pop),:);  %Eq.(12)
    for k=1:pop/2
        r3=rand();
        Si(k,:)=r3*Sea_horses_father(k,:)+(1-r3)*Sea_horses_mother(k,:);  %Eq.(13)
    end
    Sea_horses_offspring=Si;
    for i=1:pop*0.5
        %Bound the variable
        Tp=Sea_horses_offspring(i,:)>UB';
        Tm=Sea_horses_offspring(i,:)<LB';
        Sea_horses_offspring(i,:)=(Sea_horses_offspring(i,:).*(~(Tp+Tm)))+UB'.*Tp+LB'.*Tm;
        %Calculated all offspring' fitness values
         Kp1=Sea_horses(i,1);
Kd1=Sea_horses(i,2);
Ki1=Sea_horses(i,3);
Kt1=Sea_horses(i,4);
Kp2=Sea_horses(i,5);
Kd2=Sea_horses(i,6);
Ki2=Sea_horses(i,7);
Kt2=Sea_horses(i,8);
l1=Sea_horses(i,9);
l2=Sea_horses(i,10);
m1=Sea_horses(i,11);
m2=Sea_horses(i,12);
n1=Sea_horses(i,13);
n2=Sea_horses(i,14);
N1=Sea_horses(i,15);
N2=Sea_horses(i,16);
[tout] = sim('btech_ev_1');
         o2=tout.e;
        Sea_horsesFitness2(1,i)=o2;
        
    end
    %Sea horses selection
    Sea_horsesFitness=[Sea_horsesFitness1,Sea_horsesFitness2];
    Sea_horses_new=[Sea_horses_new2;Sea_horses_offspring];
    
    [~,sorted_indexes]=sort(Sea_horsesFitness);
    
    %Select the best pop and sort it
    Sea_horses=Sea_horses_new(sorted_indexes(1:pop),:);
    
    SortfitbestN = Sea_horsesFitness(sorted_indexes(1:pop));
    fitness_history(:,t)=SortfitbestN';
    population_history(:,:,t)=Sea_horses;
    Trajectories(:,t)=Sea_horses(:,1);
    %Update the optimal solution
    if SortfitbestN(1)<TargetFitness
        TargetPosition=Sea_horses(1,:);
        TargetFitness=SortfitbestN(1);
    end
    
    Convergence_curve(t)=TargetFitness;
    t = t + 1;
    display(['The Current Iteration is:',num2str(t)]);
    
    
end


time=toc;

semilogy(1:Max_iter,Convergence_curve,'color','r','linewidth',2.5);
title('Convergence curve');
xlabel('Iteration');
ylabel('Best score obtained so far')
display(['The running time is:', num2str(time)]);
display(['The best solution obtained by SHO is : ', num2str(TargetFitness)]);
display(['The best optimal sea horse of the objective funciton found by SHO is : ', num2str(TargetPosition)]);