clear all;
data_exists = exist('data.mat'); %sprawdza czy istnieje cos takiego jak data.mat

if data_exists==2   %2 oznacza ze istnieje taki plik z rozszerzeniem .mat
   load('data')
   i=size(allData, 1) + 1;
elseif data_exists==0
   i=1;
end
% NOTE 3:
% Perez and Behdinan (2007) demonstrated that the particle swarm is only
% stable if the following conditions are satisfied:
% Given that C0 = particle inertia %default 0.9
%            C1 = options.SocialAttraction
%            C2 = options.CognitiveAttraction
%    1) 0 < (C1 + C2) < 4
%    2) (C1 + C2)/2 - 1 < C0 < 1


%parametry optymalizacji - przypisanie tych zmiennych do problemu jest nizej
C1 = 1.6; %default 0.5
C2 = 2.2; %default 1.25
popSize = 100; %default 40

%parametry obiektu
global alfa1 alfa2 beta1 beta2;
alfa1 = -1.413505;
alfa2 = 0.462120;
beta1 = 0.016447;
beta2 = 0.012722;

%w sumie to teraz pso_scale juz nie jest potrzebne - dlatego =1
global pso_scale;
pso_scale=1;

%wybor algorytmu
global algorytm;
algorytm='NPL'; %GPC NPL PID oraz NO, ktory nie dziala

if strcmp (algorytm, 'NPL')
    global w10 w1 w20 w2;
    model8              %wczytanie modelu
    lb=[5/pso_scale; 5/pso_scale; 0.2];
    %ub=[100; 100; 20];
    popinit = [[5;20],[5;20], [5;20]]; %[15/pso_scale; 5/pso_scale; 3];
elseif strcmp(algorytm, 'PID')
    lb=[];
    %popinit =[[-10;5],[0;20], [0;5]];% [-0.6; 10; 0.5];
end
global u y; %do podgl¹du ostatecznych wartoœci
global kk; kk = 800;
global yzad; %Yzad;
yzad(1:round(kk/10)) = 0;
yzad(round(kk/6):round(2*kk/6)) = 1.3;
yzad(round(2*kk/6):round(3*kk/6)) = 2.2;
yzad(round(3*kk/6):round(4*kk/6)) = 0.7;
yzad(round(4*kk/6):round(5*kk/6)) = -0.1;
yzad(round(5*kk/6):round(6*kk/6)) =-0.3;


%fitnessfcn = @(x)(x(1)^2+x(2)-11)^2+(x(1)+x(2)^2-7)^2 ;
nvars = 3 ;

options.CognitiveAttraction = C1; %default 0.5
options.SocialAttraction = C2; %default 1.25
options.PopulationSize = popSize; %default 40
options.PopInitRange=popinit;%ones(nvars,1);%[-5; 5; 5];
% problem.options.KnownMin=[0,0];
problem.options = options;
problem.Aineq = []; problem.bineq =[];
problem.Aeq=[]; problem.beq =[];
problem.LB=lb; problem.UB=[];
problem.fitnessfcn.nonIcon = [] ;
problem.nonlcon = [] ;
problem.fitnessfcn = @PID_NPL_fun;%fitnessfcn;
problem.nvars = nvars; 
problem.options.DemoMode='pretty';
problem.options.PlotFcns={@psoplotbestf};%,@psoplotswarmsurf};
% problem.options.PlotFcns={@psoplotbestf};
% problem.options.Display = 'iter';



problem.options.Generations = 100;

[x_return,fmin] = pso(problem)

x = -5:0.1:5;
y_pso = -5:0.1:5;
[x1,x2] = meshgrid(x,y_pso);
y_pso = (x1.^2+x2-11).^2+(x1+x2.^2-7).^2;
figure;
contour(x1,x2,y_pso,30);
title(sprintf('x_1 = %f x_2 = %f y = %e',x_return(1),x_return(2),fmin))
xlabel('x_1','FontSize',18);
ylabel('x_2','FontSize',18);
colormap(jet)
hold on; 
plot(x_return(1),x_return(2),'x');



E = fmin;
str = sprintf('Algorytm %s    E=%e', algorytm, E);
if strcmp(algorytm, 'PID')
    K=x_return(1);
    Ti=x_return(2);
    Td=x_return(3);
    str2 = sprintf('K=%.2f  Ti=%d  Td=%.2f', K, Ti, Td);
    filename = sprintf('%s__K=%.2f__Ti=%d__Td=%.2f', algorytm, K, Ti, Td);
    filename = strrep(filename,'.',',');
else
    N = round(x_return(1)*pso_scale);
    Nu = round(x_return(2)*pso_scale);
    lambda=x_return(3);
    str2 = sprintf('N=%d  Nu=%d  lambda=%d', N, Nu, lambda);
    filename = sprintf('%s__N=%d__Nu=%d__l=%d', algorytm, N, Nu, lambda);
end

%%%%%%%%prezentacja wyników symulacji%%%%%%%%
figure; 
subplot(2,1,1);
stairs(y); hold on;
stairs(yzad,':'); xlim([0, kk]); ylim([-0.5, 2.4]);
title({'Wyjœcie regulatora',str, str2}); xlabel('k'); ylabel('wartoœæ sygna³u');
legend('wyjœcie y(k)','wartoœæ zadana');

subplot(2,1,2); 
stairs(u);hold on
title({'Sygna³ sterowania'}); xlabel('k'); ylabel('wartoœæ sygna³u'); xlim([0, kk]);

%%%%zapis wyników
if strcmp(algorytm, 'PID')
    dataVec={i, algorytm, E, K, Ti, Td, C1, C2, popSize};
elseif strcmp(algorytm, 'NPL')
    dataVec={i, algorytm, E, N, Nu, lambda, C1, C2, popSize};
end
allData(i,:)=dataVec;

save('data.mat', 'allData')