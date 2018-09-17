clear all;
A = load('derdata.mat');
x1 = A.X;
y1 = A.Y;

dy1 = diff(y1(:))./diff(x1(:));  % since Der yields dy1 and dx1 of different sizes
dx1 = x1(1:length(x1)-1);
dy2 = diff(dy1(:))./diff(dx1(:));
dx2 = dx1(1:length(dx1)-1);


figure;
subplot(1,2,1)
plot(x1,y1,dx1,dy1,dx2,dy2)
legend('Y','Y''(X)','Y''''(X)')

subplot(1,2,2)
plot(dx1,dy1,dx2,dy2)
legend('Y''(X)','Y''''(X)')




