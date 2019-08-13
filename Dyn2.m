
b_initial= [75.0; 5.0; 80.0; 60.0; 70.0; 80.0; 20.0; 90.0]
rho_ij= [0.5; 0.5; 0.5; 0.5; 0.5; 0.5]
rho= 10
A= [0 1 0 1; 1 0 1 0; 0 1 0 1; 1 0 1 0]

%delta_initial= euclidean_distances(b_initial)
[a_vector, mu_initial]= inizialize_mu(A, rho_ij, rho)
%ode_solver(a_vector, b_initial, mu_initial, rho, rho_ij)
[t,y] = ode45(@vdp1,[0 1],[b_initial; mu_initial]);
plot(y(:,1),y(:,2),'-r', y(:,3),y(:,4), '-g',y(:,5),y(:,6), '-b', y(:,7),y(:,8), '-y')
y_size= size(y);
t=[1:y_size(1)];
plot(t,y(:,9),'-r',t, y(:,10), '-g',t,y(:,11), '-b',t, y(:,12), '-y',t,y(:,13), '-c',t,y(:,14), '-m')

function delta= euclidean_distances(b)
	x_matrix = [1 0 -1 0 0 0 0 0; 1 0 0 0 -1 0 0 0; 1 0 0 0 0 0 -1 0; 0 0 1 0 -1 0 0 0; 0 0 1 0 0 0 -1 0; 0 0 0 0 1 0 -1 0];
	y_matrix = [0 1 0 -1 0 0 0 0 ; 0 1 0 0 0 -1 0 0; 0 1 0 0 0 0 0 -1; 0 0 0 1 0 -1 0 0; 0 0 0 1 0 0 0 -1; 0 0 0 0 0 1 0 -1];

	delta= sqrt((x_matrix*b).^2+(y_matrix*b).^2);
end

function [a_vector, mu_initial]= inizialize_mu(A, rho_ij, rho)
	p=[];
	for r_i= 1: size(A, 1)
		for r_j= r_i+1: size(A, 2)
			p= [p; A(r_i, r_j)];
		end
	end
	a_vector= p;
	mu_initial= a_vector.*(rho- rho_ij)+ (1-a_vector).*(rho+ rho_ij);
  end

%_______________POTENTIAL FUNCTION ___________________

function beta= beta_function(delta, rho)
	beta_ij = delta.^2-rho.^2;
	beta = prod(beta_ij);
end

%beta= beta_function(delta, rho_ij)


function gamma= gamma_function(a_vector, delta, mu, rho)
	t_square = rho^2 - mu.^2; 
	v_square = rho^2 + mu.^2; 
	gamma_ij= a_vector.*((delta.^2- t_square).^2)+ (1-a_vector).*((delta.^2- v_square).^2);
	gamma= sum(gamma_ij);
end

%gamma= gamma_function(a_vector, delta, mu_initial, rho)


function phi= artificial_function(a_vector, b, mu, rho, rho_ij)
	delta= euclidean_distances(b);
	beta= beta_function(delta, rho_ij);
	gamma= gamma_function(a_vector, delta, mu, rho);
	K=6;
	phi_cappello= gamma^K /beta;
	phi= phi_cappello^(1/K); %(phi_cappello/(1+phi_cappello))^(1/K);
end

%phi= artificial_function(a_vector, delta, mu_initial, rho, rho_ij)


function ode_solver(a_vector, b_initial, mu_initial, rho, rho_ij)
	syms b(t) [8 1];
	syms m(t) [6 1];
        syms b2 [8 1];
        syms m2 [6 1];
        
        H=  artificial_function(a_vector, b2, m2, rho, rho_ij);
        Db2 =  gradient(H, b2);
        Dm2 = gradient(H, m2);
        Db= subs(Db2, b2, b);
        Dm= subs(Dm2, m2, m);
        disp("diff")
 
        ode1 = diff(b,t) == - Db
	ode2 = diff(m,t) == - Dm
	odes = [ode1; ode2]

	cond1 = b(0) == b_initial
	cond2 = m(0) == mu_initial

	conds = [cond1 cond2]

	S = dsolve(odes, conds)

end


function dydt = vdp1(t,y)
    t
    b_initial= [75.0; 5.0; 80.0; 60.0; 70.0; 80.0; 20.0; 90.0];
    rho_ij= [0.5; 0.5; 0.5; 0.5; 0.5; 0.5];
    rho= 10;
    A= [0 1 0 1; 1 0 1 0; 0 1 0 1; 1 0 1 0];
    [a_vector, mu_initial]= inizialize_mu(A, rho_ij, rho);
    
    
    b= y(1:8,1);
    m= y(9:14,1);
    syms sb [8 1];
    syms sm [6 1];
    H = artificial_function(a_vector, sb, sm, rho, rho_ij);
    r = [-gradient(H, sb); -gradient(H, sm)];
    r1= subs(r, sb, b);
    r2= subs(r1, sm, m);
    dydt= double(r2)
   
end



