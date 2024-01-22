function genModel = parameterizedGenModel(a, b, c, q, p, N,tMax)

y0 = zeros(2*N,1);
y0(1:N) = linspace(-1,1,N);
y0(N+1:2*N) = rand(N,1)-0.5;

tspan = [0 tMax];
if c>1 && mod(c,2) == 0
	sol = ode45(@(t,y) parameterizedSystemGrad(t,y,a,b,c,q,p,N), tspan, y0);
elseif c==1
	sol = ode45(@(t,y) paramL1Grad(t,y,a,b,c,q,p,N), tspan, y0);
else
	sol = ode45(@(t,y) paramGeneralLcGrad(t,y,a,b,c,q,p,N), tspan, y0);
end

xpoints = sol.y(1:N,:);
ypoints = sol.y(N+1:2*N,:);
m = size(xpoints,2);

xLow = min(min(xpoints));
xHigh = max(max(xpoints));
yLow = min(min(ypoints));
yHigh = max(max(ypoints));

stringName1 = sprintf('p=%d,q=%d,c=%d,a=%d,b=%d,N=%d,tMax=%d,iteration=%d',p,q,c,a,b,N,tMax,1);

figure(1);
plot(xpoints(:,1),ypoints(:,1),'k.','MarkerSize',10);
axis([-1 1 -1 1])
axis equal
axis on
print('-dpng','-r150',stringName1)
close(1)

for i=0:m
    xcol = xpoints(:,m-i);
    ycol = ypoints(:,m-i);
    if ~(any(isnan(xcol)) || any(isnan(ycol)))        
        figure(1);
        plot(xpoints(:,m-i),ypoints(:,m-i),'k.','MarkerSize',10);
        axis([xLow xHigh yLow yHigh])
        axis equal
        axis on
	stringName2 = sprintf('p=%d,q=%d,c=%d,a=%d,b=%d,N=%d,tMax=%d,iteration=%d',p,q,c,a,b,N,tMax,m-i);
        print('-dpng','-r150',stringName2)
        close(1)
        break
    end
end