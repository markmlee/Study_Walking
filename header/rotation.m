p1 = [0.615; 0.553];
p2 = [0.714; 0.553];
p3 = [0.714; 0.503];
p4 = [0.614; 0.503];

th = 60*pi/180;

R = [cos(th) -sin(th); sin(th) cos(th)];

pNew1 = R*p1;
pNew2 = R*p2;
pNew3 = R*p3;
pNew4 = R*p4;

x = [p1(1) p2(1) p3(1) p4(1)];
y = [p1(2) p2(2) p3(2) p4(2)];

plot(x,y,'o')

hold on 

xN = [pNew1(1) pNew2(1) pNew3(1) pNew4(1)];
yN = [pNew1(2) pNew2(2) pNew3(2) pNew4(2)];

scatter(xN,yN,'o', 'k')
