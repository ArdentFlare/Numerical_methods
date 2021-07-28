proc iml;
/* #############################: EUELERS METHOD */
/* create module for equation */
start f(y,t);
return y - t##2 +1;
finish;


/* create Euler's method module  */
start Euler(t0,t1,n);
/* declare given values */
a = 0;
b = 2;
y0 = 0.5;
/* stepsize equation */
h = (b - a)/n;
wo = y0;
/* create table to show Euler's algorithm solutions */
T = j(n+1,5);

do i=1 to n+1;
T[i,1]=i-1;
end;
/* fill in  rows using Euler's Method */
T[1,2]=t0;
T[1,3]=wo;
T[1,4]=h#f(wo,t0);
T[1,5]=wo+h#f(wo,t0);

do j = 2 to n+1;
T[j,2] = t1;
T[j,3] = T[j-1,5];
T[j,4] = h#f(T[j-1,5],T[j,2]);
T[j,5] = T[j-1,5] + h#f(T[j-1,5],T[j,2]);
t1 = t1+h;
end;

l = {"i" "ti" "wi" "hf(ti,wi)" "wi+1"};
name = {"Eulers Method"};
print T[colname =l label = name];
return T;
finish;

E = Euler(0,0.5,4);
print E;


/* ################: NEWTONS METHOD FOR SYSTEMS */

start f(x,y);
y1 =x#y - sin(x##2+y##2)-log(22/7);
return y1;
finish;
pie = constant("pi");

start g(x,y);
y2 = exp(x#y)+x##2+y##2-9;
return y2;
finish;

start fx(x,y);
y3 = y - cos(x##2+y##2)#2#x;
return y3;
finish;

start fy(x,y);
y4 =  x - cos(x##2+y##2)#2#y;
return y4;
finish;

start gx(x,y);
y5 = exp(x#y)#(y)+2#x ;
return y5;
finish;

start gy(x,y);
y6 =  exp(x#y)#(x)+2#y;
return  y6;
finish;

a = 2.5;
b = 0.5;
n = 12;

start NRaph(a,b,n,e);

t = j(n,2);
t[1,1] = a;
t[1,2] = b;

do i = 2 to n;

sol1 = t[i-1,1]//t[i-1,2];
print sol1;
f = f(t[i-1,1],t[i-1,2])//g(t[i-1,1],t[i-1,2]);
print f;
J = (fx(t[i-1,1],t[i-1,2])||fy(t[i-1,1],t[i-1,2]))//(gx(t[i-1,1],t[i-1,2])||gy(t[i-1,1],t[i-1,2]));
print j;
Jinverse = inv(j);
print Jinverse;
sol2 = sol1 -inv(J)*f;
print sol2;
a2 = sol2[1,1];
b2 = sol2[2,1];

t[i,1] = a2;
t[i,2] = b2;

end;
return t;
finish;

tay = NRaph(a,b,n,0.00001);
print tay;

/* ###############: SIMPSONS RULE */

start f(x);
return sin(x##2)+sin(x##2)+((2/(3#constant('pi')))##0.5)#x;
finish;

start f4(x);
return 16#x##4#sin(x##2)-48#x##2#cos(x##2)-12#sin(x##2);
finish;

x = do(0,(((3#constant('pi'))/2)##0.5),0.05);
y = f4(x);
call Series(x,y) grid = {X Y};
/* if the max is not at "a" or "b"(difficult to tell where) then find next derivative and use froot function */
a=0;
b=((3#constant('pi'))/2)##0.5;
e4 = abs(f4(b)#(b-a)##5/2880);
print e4;
/* no. of iterations #Uses error module above*/

start ni(e4,e);
return sqrt(e4/e);
finish;

n = ni(e4,0.00001);
print n;

/* SIMPSONS RULE */

start simpson(a,b);
return (b-a)#(f(a)+4#f((a+b)/2)+f(b))/6;
finish;

s = simpson(0,(3#constant('pi')/2)##0.5);
print s;

start cSimpson(a,b,n);
p = {"i" "xi(h)" "f(xi)" "simpson"};
Tay = j(n+1,4);
h = (b-a)/n;

do i=1 to n+1;
Tay[i,1] = i-1;
Tay[i,2] = a+(i-1)#h;
Tay[i,3] = f(a+(i-1)#h);
end;

Tay[1,4] = Tay[1,3];
Tay[n+1,4] = Tay[n+1,3];

do i = 2 to n by 2;
Tay[i,4] = 4#Tay[i,3];
end;
do i=3 to n-1 by 2;
Tay[i,4] = 2#Tay[i,3];
end;


sumsim = sum(Tay[,4])#(h/3);
print sumsim;
print Tay[colname = p];
finish;

run cSimpson(0,(3#constant('pi')/2)##0.5,ceil(n));


/* #################: SECANT METHOD */
 
start g(x);
y = cos(x) - x;
return y;
finish;
pie = constant("pi");

start p(pn,pnminus);
pnplus = pn - (g(pn)*(pn - pnminus))/(g(pn) - g(pnminus));
return pnplus;
finish;

pa = 0.5;
pb = pie/4;
n = 7;

T = j(n,3);
colhead = {"i" "pn" "g(pn)"};
TableHead = {"SECANT METHOD"};
T[1,2] = pa;
T[1,3] = g(pa);
T[2,2] = pb;
T[2,3] = g(pb);
 
do i = 1 to n;
T[i,1] = i;
end;

do j = 3 to n;
pnplusone = p(pb,pa);
T[j,2] = pnplusone;
T[j,3] = g(pnplusone);
pa = pb;
pb = pnplusone;
end;

print T[colname = colhead LABEL = TableHead];
  

/* ###############: RUNGE KUTTA MIDPOINT */

/* create module for equation */
start f(y,t);
return y-t##2+1;
finish;

/* module for order 2 */
start f2(y,t);
return (y-t##2+1)-2#t;
finish;


/* create Euler's method module  */
start RKMid(t0,t1,n,y0);

/* stepsize equation */
h = (t1 - t0)/n;
print h;
wo = y0;

/* create table to show Euler's algorithm solutions */
T = j(n+1,8);

do i=1 to n+1;
T[i,1]=i-1;
end;
/* fill in  rows using Euler's Method */
do j = 1 to n+1;
T[j,2]=t0;
T[j,3]=t0+(h/2);
T[j,4]=wo;
T[j,5]=f(wo,t0);
T[j,6]=wo+(h/2)#f(wo,t0);
T[j,7]=f(wo+(h/2)#f(wo,t0),t0+h/2);
T[j,8]=wo+(h/2)#f(wo+(h/2),t0+(h/2)#f(wo,t0));
t0=t0+h;
wo=T[j,8];

end;

l = {"i" "ti" "ti+h/2" "wi" "f(ti,wi)" "wo+(h/2)f(ti,wi)" "f(ti+h/2,wi+(h/2)f(ti,wi)" "wi+1"};
name = {"Runge Kutta Midpoint"};
print T[colname =l label = name];

finish;
run RKMid(0,2,10,0.5);









 




 

