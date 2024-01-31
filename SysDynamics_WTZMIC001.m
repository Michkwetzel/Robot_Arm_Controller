%% Control of a 3-Link arm
%%
% Setting up symbolic variables
syms th0 th1 th2 th3 dth0 dth1 dth2 dth3 ddth0 ddth1 ddth2 ddth3 L g m t1 t2 t3 I I0 b
q = [th1; th2; th3];
dq = [dth1; dth2; dth3];
ddq = [ddth1; ddth2; ddth3];
% *Define positions, velocities, energies*

%% Calculate Position Vectors
% Assume the COM is the origin
p1 = [L/2*cos(th1); L/2*sin(th1);0]
p2 = [L*cos(th1); L*sin(th1);0] + [(L/2)*cos(th1 + th2); (L/2)*sin(th1+th2);0]
p3 = [L*cos(th1); L*sin(th1);0] + [L*cos(th1 + th2); L*sin(th1+th2);0] + [(L/2)*cos(th1+th2+th3); (L/2)*sin(th1+th2+th3);0]

%% Velocity
v1 = simplify(jacobian(p1,q)*dq);
v2 = simplify(jacobian(p2,q)*dq);
v3 = simplify(jacobian(p3,q)*dq);


%% Body Energy
% Kinetic
T1 = 0.5*m*transpose(v1)*v1 + 0.5*I*(dth1)^2;
T2 = 0.5*m*transpose(v2)*v2 + 0.5*I*(dth1 + dth2)^2;
T3 = 0.5*m*transpose(v3)*v3 + 0.5*I*(dth1 + dth2 + dth3)^2;

Ttot = T1 + T2 + T3;
Ttot = simplify(Ttot);

%% Potential energy - dependent on g so if you're in space...
V = 0;
Vtot = 0; %We are in actual space thus there is no potential energy

%% Mass Matrix 
% Verbose way but could also use M = hessian(Ttot, dq)
for i = 1 : length(q)
    for j = 1 : length(q)
        M(i,j) = diff(diff(Ttot,dq(i)),dq(j));
    end
end

M = simplify(M)

%% Derivative of Mass Matrix  
dM = sym(zeros(length(M),length(M)));
for i=1:length(M)
    for j=1:length(M)
        dM(i,j) = jacobian(M(i,j),q)*dq;
    end
end
dM = simplify(dM)

%% C Matrix  
% Contains the centrifugal and coriolis accelerations
C = dM*dq - transpose(jacobian(Ttot,q));
C = simplify(C)

%% B input matrix 
B = eye(3);
% t1, t2 and t3 are the applied torques at each link
u = [t1; t2; t3];

%M*ddq + C = B*u + Q

b = 4.3179
Q = [-b*dth1 ; -b*dth2 ; -b*dth3]
ddq_equation = simplify(M \ (-C - G + B*u + Q))

% Feedback Linearization

%% Output variable - we want to control the angles of the links

Y = [th1;th2;th3]

J = jacobian(Y,q)
J = simplify(J)

%dJ = sym('dJ', [3 3])
dJ(:,1) = jacobian(J(:,1),q)*dq
dJ(:,2) = jacobian(J(:,2),q)*dq
dJ(:,3) = jacobian(J(:,3),q)*dq

% The following code is to calculate the tau equatiosn symbolically
syms Kp Kd e1 e2 e3 de1 de2 de3 v1 v2 v3
b = 4.3179  %coefficient of friction
Q = [-b*dth1 ; -b*dth2 ; -b*dth3]   %Frictino force
e = [e1;e2;e3]      %error terms
de = [de1;de2;de3]    %dirrivative of error

    
if abs(det(J)) < 0.01       %this is to avoid infinite errors
    tau = M*(J*(Kp*e + Kd*de - dJ*dq)) - Q + C
else
    tau = M*(J\(Kp*e + Kd*de - dJ*dq)) - Q + C
end

%%
% Finding the Value of b

%use for loop to get all possible b values. then average that.
%the 1st value of b is left out because it throws off the calculation.

b_array = zeros(1,201)
Q = M*ddq + C - B*u
Q = subs(Q,{t1,t2,t3},{0.1,0.1,0.1})
Q = subs(Q,{L,g,m,I},{0.5,9.81,1,0.1})

for t = 2:1:length(out.dth1)
    Q_array = Q
    Q_array = subs(Q_array,{th1,th2,th3},{out.th1(t),out.th2(t),out.th3(t)});
    Q_array = subs(Q_array,{dth1,dth2,dth3},{out.th1(t),out.dth2(t),out.dth3(t)});
    Q_array = subs(Q_array,{ddth1,ddth2,ddth3},{out.ddth1(t),out.ddth2(t),out.ddth3(t)});

    %Q_eval = simplify(Q_eval)

    b_array(t) = Q_array(1)/(-out.dth1(t)) + Q_array(2)/(-out.dth2(t)) + Q_array(3)/(-out.dth3(t))
    b_array(t) = b_array(t)/3
    
end

b_final = mean(b_array) = 4.3179