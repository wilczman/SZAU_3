clear all;

fitnessfcn = @(x)(x(1)^2+x(2)-11)^2+(x(1)+x(2)^2-7)^2 ;
nvars = 2 ;

options.PopInitRange=[-5; 5];
% problem.options.KnownMin=[0,0];
problem.options = options;
problem.Aineq = []; problem.bineq =[];
problem.Aeq=[]; problem.beq =[];
problem.LB=[]; problem.UB=[];
problem.fitnessfcn.nonIcon = [] ;
problem.nonlcon = [] ;
problem.fitnessfcn = fitnessfcn;
problem.nvars =2; 
problem.options.DemoMode='pretty';
problem.options.PlotFcns={@psoplotbestf,@psoplotswarmsurf};
% problem.options.PlotFcns={@psoplotbestf};
% problem.options.Display = 'iter';
problem.options.Generations = 100;

[x_return,fmin] = pso(problem)

x = -5:0.1:5;
y = -5:0.1:5;
[x1,x2] = meshgrid(x,y);
y = (x1.^2+x2-11).^2+(x1+x2.^2-7).^2;
figure;
contour(x1,x2,y,30);
title(sprintf('x_1 = %f x_2 = %f y = %e',x_return(1),x_return(2),fmin))
xlabel('x_1','FontSize',18);
ylabel('x_2','FontSize',18);
colormap(jet)
hold on; 
plot(x_return(1),x_return(2),'x');