% compute and save the date to file
L = 2*pi;
N_theta = 64;
a = 0.025;
b = 10*a;
M = 401;

warning('off');
theta = (L/N_theta)*[0:N_theta-1]';

%% N=8;
N=8;
Eps = 0.01*a;

f = exp(cos(theta));
name = 'expcos100padeN8';
FE_applicationPade(M,f,N_theta,a,b,N,Eps,L,theta,name);
name = 'expcos100N8';
FE_application(M,f,N_theta,a,b,N,Eps,L,theta,name);

f = cos(2*theta);
name = 'cos2100padeN8';
FE_applicationPade(M,f,N_theta,a,b,N,Eps,L,theta,name);
name = 'cos2100N8';
FE_application(M,f,N_theta,a,b,N,Eps,L,theta,name);

f = cos(4*theta);
name = 'cos4100padeN8';
FE_applicationPade(M,f,N_theta,a,b,N,Eps,L,theta,name);
name = 'cos4100N8';
FE_application(M,f,N_theta,a,b,N,Eps,L,theta,name);

f = cos(8*theta);
name = 'cos8100padeN8';
FE_applicationPade(M,f,N_theta,a,b,N,Eps,L,theta,name);
name = 'cos8100N8';
FE_application(M,f,N_theta,a,b,N,Eps,L,theta,name);


Eps = 0.1*a;
f = exp(cos(theta));
name = 'expcos10padeN8';
FE_applicationPade(M,f,N_theta,a,b,N,Eps,L,theta,name);
name = 'expcos10N8';
FE_application(M,f,N_theta,a,b,N,Eps,L,theta,name);

f = cos(2*theta);
name = 'cos210padeN8';
FE_applicationPade(M,f,N_theta,a,b,N,Eps,L,theta,name);
name = 'cos210N8';
FE_application(M,f,N_theta,a,b,N,Eps,L,theta,name);

f = cos(4*theta);
name = 'cos410padeN8';
FE_applicationPade(M,f,N_theta,a,b,N,Eps,L,theta,name);
name = 'cos410N8';
FE_application(M,f,N_theta,a,b,N,Eps,L,theta,name);

f = cos(8*theta);
name = 'cos810padeN8';
FE_applicationPade(M,f,N_theta,a,b,N,Eps,L,theta,name);
name = 'cos810N8';
FE_application(M,f,N_theta,a,b,N,Eps,L,theta,name);


%% N=12;
N=12;
Eps = 0.01*a;

f = exp(cos(theta));
name = 'expcos100padeN12';
FE_applicationPade(M,f,N_theta,a,b,N,Eps,L,theta,name);
name = 'expcos100N12';
FE_application(M,f,N_theta,a,b,N,Eps,L,theta,name);

f = cos(2*theta);
name = 'cos2100padeN12';
FE_applicationPade(M,f,N_theta,a,b,N,Eps,L,theta,name);
name = 'cos2100N12';
FE_application(M,f,N_theta,a,b,N,Eps,L,theta,name);

f = cos(4*theta);
name = 'cos4100padeN12';
FE_applicationPade(M,f,N_theta,a,b,N,Eps,L,theta,name);
name = 'cos4100N12';
FE_application(M,f,N_theta,a,b,N,Eps,L,theta,name);

f = cos(8*theta);
name = 'cos8100padeN12';
FE_applicationPade(M,f,N_theta,a,b,N,Eps,L,theta,name);
name = 'cos8100N12';
FE_application(M,f,N_theta,a,b,N,Eps,L,theta,name);


Eps = 0.1*a;
f = exp(cos(theta));
name = 'expcos10padeN12';
FE_applicationPade(M,f,N_theta,a,b,N,Eps,L,theta,name);
name = 'expcos10N12';
FE_application(M,f,N_theta,a,b,N,Eps,L,theta,name);

f = cos(2*theta);
name = 'cos210padeN12';
FE_applicationPade(M,f,N_theta,a,b,N,Eps,L,theta,name);
name = 'cos210N12';
FE_application(M,f,N_theta,a,b,N,Eps,L,theta,name);

f = cos(4*theta);
name = 'cos410padeN12';
FE_applicationPade(M,f,N_theta,a,b,N,Eps,L,theta,name);
name = 'cos410N12';
FE_application(M,f,N_theta,a,b,N,Eps,L,theta,name);

f = cos(8*theta);
name = 'cos810padeN12';
FE_applicationPade(M,f,N_theta,a,b,N,Eps,L,theta,name);
name = 'cos810N12';
FE_application(M,f,N_theta,a,b,N,Eps,L,theta,name);

