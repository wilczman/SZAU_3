x = -5:0.1:5;
y = -5:0.1:5;
[x1,x2] = meshgrid(x,y);
y = (x1.^2+x2-11).^2+(x1+x2.^2-7).^2;
figure;
contour(x1,x2,y,30);
xlabel('x_1','FontSize',18);
ylabel('x_2','FontSize',18);
colormap(jet)
figure;
meshc(x1,x2,y);
xlabel('x_1','FontSize',18);
ylabel('x_2','FontSize',18);
zlabel('y','FontSize',17)
colormap(jet)
 