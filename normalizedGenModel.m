tic

global p q a b c N

% Interaction kernel parameters
% q > p > -2
p = 1;
q = 6;

% Norm
c = 2;

% Elliptical parameters
a = 1/3;
b = 1;

% Number of particles
N = 1000;

y0 = zeros(2*N,1);
y0(1:N) = linspace(-1,1,N);
y0(N+1:2*N) = rand(N,1)-0.5;

% Time length
tMax = 100;

tspan = [0 tMax];
sol = ode45(@normalizedODEsolverSystemGrad, tspan, y0);

xpoints = sol.y(1:N,:);
ypoints = sol.y(N+1:2*N,:);
m = size(xpoints,2);

xLow = min(min(xpoints));
xHigh = max(max(xpoints));
yLow = min(min(ypoints));
yHigh = max(max(ypoints));

figure(1);
plot(xpoints(:,1),ypoints(:,1),'k.','MarkerSize',10);
axis([-1 1 -1 1])
axis equal
axis on
print('-dpng','-r150',num2str(1))
close(1)

figure(1);
plot(xpoints(:,m),ypoints(:,m),'k.','MarkerSize',10);
axis([xLow xHigh yLow yHigh])
axis equal
axis on
print('-dpng','-r150',num2str(m))
close(1)

E = zeros(N,m);

for i=1:N
    for j=1:N
        if j~=i
            x=xpoints(i,:)-xpoints(j,:);
            y=ypoints(i,:)-ypoints(j,:);           
            % J is norm
            J = ( (a*x).^c + (b*y).^c ).^(1/c);
            K = 1/q*(J.^q);
            if p~=0
                E(i,:) = E(i,:) + (-1/p*(J.^p)+K);
            else
                E(i,:) = E(i,:) + (-log(J)+K); 
            end
        end
    end
end

Etotal = (1/(a*a*b*b))*sum(E,1);

figure(2); clf; 
plot(sol.x,Etotal,'b');
print('energy','-dpng','-r150')

toc