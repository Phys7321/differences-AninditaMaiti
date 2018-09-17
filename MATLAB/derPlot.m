clear all;

x = linspace(-2,4,400);
y = @(x) sin(1 ./(x .*(2-x))).^2;
[dyc,dxc] = Der(y,x,'dc');
[dyf,dxf] = Der(y,x,'df');
[dyb,dxb] = Der(y,x,'db');

[dyc2,dxc2] = Der(y,x,'2dc');
[dyc3,dxc3] = Der(y,x,'3dc');
[dyc4,dxc4] = Der(y,x,'4dc');
[dyc5,dxc5] = Der(y,x,'5dc');


[dyf2,dxf2] = Der(y,x,'2df');
[dyb2,dxb2] = Der(y,x,'2db');

figure
subplot(2,2,1)
plot(x,y(x))
legend('y(x)','location','southoutside')

subplot(2,2,2)
plot(dxc,dyc,dxf,dyf,dxb,dyb)
legend('central','forward','backward','location','southoutside')
title('Comparison of 1st derivatives')

subplot(2,2,3)
plot(dxc2,dyc2,dxf2,dyf2,dxb2,dyb2)
legend('central','forward','backward','location','southoutside')
title('Comparison of 2nd derivatives')

subplot(2,2,4)
plot(dxc,dyc,dxc2,dyc2,dxc3,dyc3,dxc4,dyc4,dxc5,dyc5)
legend('y''(x)','y''''(x)','y''''''(x)','y''''''''(x)','y''''''''''(x)','location','southoutside')
title('Upto 5th order derivatives by central differentiation')

