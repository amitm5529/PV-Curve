%190588_Parikshit Ghosh
Tc = 460.4;   %Tc in K                                                                                                Vc = 0.25;    
Pc = 5270;    %Pc in KPa
Vc = 0.20;    %Vc in L/mol
R = 8.314;    %R in L.KPa/(K.mol)
w = 0.190;    %acentric factor

a = 0.42747*R^2*Tc^2/Pc;
b = 0.08664*R*Tc/Pc;

plot(Vc, Pc, '.', 'color', 'red');
hold on

V = [];
P = [];

for p = 3000:10:8000     %Plotting Pv curve at Tc
    Tr = Tc/Tc; 
    alp = (1+ (0.48 + 1.574*w - 0.176*w^2)*( 1 - Tr^0.5))^2;
    v = [1 -R*Tc/p a*alp/p-b*R*Tc/p-b^2 -a*alp*b/p];
    v = roots(v);
    v = v(imag(v) == 0);
    V = [V v];
    P = [P p];
end
plot(V ,P,'color','blue');


for T = Tc+10: 15 : Tc+85   %Plotting Pv curve above Tc
    V = [];
    P = [];
    
    for p = 3000:10:8000 
        Tr = T/Tc;
        alp = (1+ (0.48 + 1.574*w - 0.176*w^2)*( 1 - Tr^0.5))^2;
        v = [1 -R*T/p a*alp/p-b*R*T/p-b^2 -a*alp*b/p];
        v = roots(v);
        v = v(imag(v) == 0);
        V = [V v];
        P = [P p];
    end
    plot(V, P, 'color', 'red');
end


Psat = [];

for T = Tc-3 : -3 : Tc-36       %Finding Psat by convergence. Reference: Asg3 Q4
   pg = Pc-100;
   Eold = 100;
   e = 0.01;
   
   while 1
       Tr = T/Tc;
       alp = (1+ (0.48 + 1.574*w - 0.176*w^2)*( 1 - Tr^0.5))^2;
       v = [1 -R*T/pg a*alp/pg-b*R*T/pg-b^2 -a*alp*b/pg];
       v = roots(v);
       l = 1/v(3);   %l represents density of liquid
       g = 1/v(1);   %g represents desnity of vapor
       
       
       LHS = -log(1 - l*b) - (a*alp)*log(1 + l*b)/(b*R*T) + pg/(l*R*T) - 1 - log(pg/(l*R*T));
       RHS = -log(1 - g*b) - (a*alp)*log(1 + g*b)/(b*R*T) + pg/(g*R*T) - 1 - log(pg/(g*R*T));
       
       Enew = LHS - RHS;
       
       if abs(Enew) < e
           Psat = [Psat pg];
            break;
       elseif abs(Enew) < abs(Eold)
           pg = pg - 10;
           Eold = Enew;
       end
   end
end

n=0;
i=1;
V1 = [Vc];
V2 = [Vc];
Pt = [Pc];
for T = Tc-3: -3: Tc-36    %Constructing the dome
    Tr = T/Tc;
    alp = (1+ (0.48 + 1.574*w - 0.176*w^2)*( 1 - Tr^0.5))^2;
    v = [1 -R*T/Psat(i) a*alp/Psat(i)-b*R*T/Psat(i)-b^2 -a*alp*b/Psat(i)];
    v = roots(v);
    V1 = [V1 v(3)];
    V2 = [V2 v(1)];
    Pt = [Pt Psat(i)];
    n = n+1;
    i = i+1;
end

plot(real(V1), Pt, 'color', 'green', 'LineWidth', 4);
plot(real(V2), Pt, 'color', 'green', 'LineWidth', 4);
     

for m = 2:n        %inside the dome
    plot( real([real(V1(m)) real(V2(m))]), [Pt(m),Pt(m)], 'b');
end

     j=1;
for T = Tc-3 : -3 : Tc-36     %liquid region
    V = [];
    P = [];
    
    for p = Psat(j) + 3 : 100 : 8000
        Tr = T/Tc;
        alp = (1+ (0.48 + 1.574*w - 0.176*w^2)*( 1 - Tr^0.5))^2;
        v = [1 -R*T/p a*alp/p-b*R*T/p-b^2 -a*alp*b/p];
        v = roots(v);
        V = [V v(3)];
        P = [P p];
    end
    j = j+1;
    plot(V, P, 'color', 'blue');
end


k = 1;
for T = Tc-3 : -3 : Tc-36       %vapor region
    V = [];
    P = [];
    
    for p = Psat(k) - 3 : -100 : 3000
        Tr = T/Tc;
        alp = (1+ (0.48 + 1.574*w - 0.176*w^2)*( 1 - Tr^0.5))^2;
        v = [1 -R*T/p a*alp/p-b*R*T/p-b^2 -a*alp*b/p];
        v = roots(v);
        V = [V v(1)];
        P = [P p];
    end
    k = k+1;
    plot(real(V), P, 'color', 'blue');
end

title('Pv curve for Soava Redlich Kwong (SRK) EOS');
xlabel('v(L/mol)');
ylabel('P(KPa)');
hold off

    
    
    
    
    

        
        
    

    
    
    


        
    
    
    
    
    
    
    