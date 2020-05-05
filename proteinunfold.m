function proteinunfold(k,n)
A=  [3,2,2,3,3,2,2,3,4,4,4,4,5,5,5,5]; %native protein structure
B = [4,4,3,3,2,2,1,1,1,2,3,4,4,3,2,1];


global E;
E= k; %defining global variable for energy of interaction
h= zeros(n,1);% initializing array for energy matrix
v= zeros(n,1);% initializing array for iterations

for t=1:n
C=A; % lattice point matrix for previous iteration
D=B;
global i; % defining global variable for randomly picked lattice point
i= randi(16);
if (i==16) || (i==1) % if corner most amino acid undergo cornermove
 [A,B]= cornermove(A,B,i);
else
   [A,B]=  crankshaft(A,B,i);
   
end
 [A,B]=accept(A,B,C,D,E); % calling function to accept or reject the transformation
 disp(A);
 disp(B);
 energy(A,B,E)
plot(A,B,'-o'); %ploting modified structure after every iteration
axis([-5 10 -5 10]);
e= energy(A,B,E);
txt = ['Energy:' num2str(e) ' units   moves:' num2str(t)];
text(-4,10.5,txt);
drawnow limitrate
disp(t); % displaying no of iteration in command window

h(t)= energy(A,B,E); % updating matrix for energy at each iteration
v(t)= t; % updating iteration matrix
end
plot(v,h); %plotting energy vs monte carlo steps
axis([1 1000000 -9 0]);
%plot(A,B,'-o'); %ploting modified structure after every iteration
%axis([-5 10 -5 10]);
%txt = ['Energy:' num2str(energy(A,B,E)) ' units   moves:' num2str(t)];

end




function [A,B]= crankshaft(A,B,i)
x= A(i-1)+A(i+1)-A(i); % equivalent to finding fourth coordiante of square
y= B(i-1)+B(i+1)-B(i);
c= zeros(16,1);
for p= 1:16 % checking for overlap
    if x==A(p) && y==B(p)
        c(p)=1;
    else
        c(p)=0;
    end
end
if c == zeros(16,1)
    A(i)= x;
    B(i)= y;
end 
end


function [A,B] = cornermove(A,B,i)
v=[pi/2,-pi/2];
theta= v(randi(2)); %randomly choosing clockwise or anticlockwise rotation
if i==16
    a=15;
else
    a=2;
    
end
x= round(cos(theta)*(A(i)-A(a))- sin(theta)*(B(i)-B(a))+A(a)); % rotation of lattice point(corner point) wrt adjacent point
y= round(sin(theta)*(A(i)-A(a))- cos(theta)*(B(i)-B(a))+B(a));

c= zeros(16,1);
for p= 1:16 %checking for overlap
    if x==A(p) && y==B(p)
        c(p)=1;
    else
        c(p)=0;
    end
end
if c == zeros(16,1)    
    A(i)= x;
   
    B(i)= y;
     
end
end




function [A,B]= accept(A,B,C,D,E)
w= exp(-(energy(A,B,E)- energy(C,D,E))); % defining w on the basis of energy of interaction
if w<1 % accepting move if w>1 if w<1 choose r randomly
    r= rand(1);
    if r>w % if r>w accept the moove
        A=C;
        B=D;
    end
end
end

function z= energy(A,B,E)
% energy function to calculate total energy of interaction respective to native state
u=0;
 if distl(A(1),B(1),A(4),B(4)) == 1  %corresponding interactions in native state 
     u=u+1; % if specific interaction present in modified matrix increment the counter by 1
 end
 if distl(A(1),B(1),A(12),B(12)) == 1
     u=u+1;
 end
 
 if  distl(A(4),B(4),A(11),B(11))==1
   u=u+1;
 end
     
 if  distl(A(5),B(5),A(8),B(8)) ==1
     u=u+1;
 end
 
 if  distl(A(5),B(5),A(10),B(10)) == 1
     u=u+1;
 end
 if  distl(A(9),B(9),A(16),B(16)) == 1
     u=u+1;
 end
 if distl(A(10),B(10),A(15),B(15)) == 1
     u=u+1;
 end
 if   distl(A(11),B(11),A(14),B(14)) == 1
     u=u+1;
 end
 if  distl(A(3),B(3),A(6),B(6)) == 1
     u=u+1;
 end
 z= u*E;
end

function y= distl(x1,y1,x2,y2) % function to calculate distance between two lattice points
y= round(sqrt((x1-x2)^2 + (y1-y2)^2));
end



