clear
close all

iter = 1e5; %Iteraciones
dt = .01; %Paso
%Opciones de ejecucion del programa 
random_initial_conditions = false; 
freeze = true; %Evita calcular trayectorias de particulas fuera de la botella
freeze_threshold = 4; %Alejamiento maximo proporcional a "a" antes del congelamiento
min_particles = 1; %Minimo de particulas en la botella para continuar la simulacion
testAxis = false; %Test para revisar el comportamiento en cada eje
testE = false; %Test para revisar la interaccion electrica

botella = true; %Campo magnetico externo
particle_E = true; %Interaccion electrica entre las particulas
particle_B = true; %Interaccion magnetica entre las particulas

%Opciones de animacion
animation = true; 
speed = 50; %Velocidad de animacion
rotSpeed = 5; %Velocidad de rotacion de la camara
tail = 5000; %Puntos ploteados de cada particula en un dado instante
lag = speed*50; lag_count = 0; %NOT WORKING Retraso para detener el programa si todas las particulas escapan

k = 1; %1/4pi*eps0 
k2 = 1; %mu0/4pi
q = 1; m = 1; %Propiedades de las particulas
c = 7; a = 5; %Propiedades del campo magnetico externo
if random_initial_conditions
    rng(100) %Semilla
    N = 7; %Numero de particulas
    L = .9*a; %Maxima distancia inicial de las particulas
    pos = 2*L*rand(N,3)-L; 
    vel = .2*L*rand(N,3)-.1*L; 
    x0 = pos(:,1); y0 = pos(:,2); z0 = pos(:,3); 
    dx0 = vel(:,1); dy0 = vel(:,2); dz0 = vel(:,3);
else
    %Condiciones iniciales manuales
    x0 = [1,-1]; y0 = [1,-1]; z0 = [1,-1];
    dx0 = [.1,-.1]; dy0 = [.1,-.1]; dz0 = [.1,-.1];
    N = length(x0);
end
save ParteC_1

%Se definen los vectores a utilizar
[x, y, z, dx, dy, dz, ddx, ddy, ddz] = deal(zeros(N,iter+1));
%Se asignan condiciones inciales
x(:,1) = x0; y(:,1) = y0; z(:,1) = z0;
dx(:,1) = dx0; dy(:,1) = dy0; dz(:,1) = dz0;

%Se define la interaccion electrica entre particulas
if particle_E 
    Ex = @(x1,y1,z1,x2,y2,z2) ...
        -k*(x2-x1) / (((x2-x1)^2+(y2-y1)^2+(z2-z1)^2)^(3/2));
    Ey = @(x1,y1,z1,x2,y2,z2) ...
        -k*(y2-y1) / (((x2-x1)^2+(y2-y1)^2+(z2-z1)^2)^(3/2));
    Ez = @(x1,y1,z1,x2,y2,z2) ...
        -k*(z2-z1) / (((x2-x1)^2+(y2-y1)^2+(z2-z1)^2)^(3/2));
else
    Ex = @(x1,y1,z1,x2,y2,z2) 0; Ey = @(x1,y1,z1,x2,y2,z2) 0; Ez = @(x1,y1,z1,x2,y2,z2) 0; 
end

%Se define el campo magnetico externo
if botella 
    Bx = @(x,y,z) ...
        (c*a^(3/2))*(((x+a)/(((x+a)^2+y^2+z^2)^(2))) - ((x-a)/(((x-a)^2+y^2+z^2)^(2))));
    By = @(x,y,z) ...
        (c*a^(3/2))*((y/(((x+a)^2+y^2+z^2)^(2))) - (y/(((x-a)^2+y^2+z^2)^(2))));
    Bz = @(x,y,z) ...
        (c*a^(3/2))*((z/(((x+a)^2+y^2+z^2)^(2))) - (z/(((x-a)^2+y^2+z^2)^(2))));
elseif testE
    Bx = @(x,y,z) 0; By = @(x,y,z) 0; Bz = @(x,y,z) 0;
else
    Bx = @(x,y,z) 1; By = @(x,y,z) 0; Bz = @(x,y,z) 0; %Campo magnetico uniforme para pruebas
end

%Se define la interaccion magnetica entre particulas
if particle_B
    Bpx = @(x1,y1,z1,x2,y2,z2,dx2,dy2,dz2) ...
        k2*((dy2*(z2-z1)-dz2*(y2-y1))/((x2-x1)^2+(y2-y1)^2+(z2-z1)^2)^(3/2));
    Bpy = @(x1,y1,z1,x2,y2,z2,dx2,dy2,dz2) ...
        k2*((dz2*(x2-x1)-dx2*(z2-z1))/((x2-x1)^2+(y2-y1)^2+(z2-z1)^2)^(3/2));
    Bpz = @(x1,y1,z1,x2,y2,z2,dx2,dy2,dz2) ...
        k2*((dx2*(y2-y1)-dy2*(x2-x1))/((x2-x1)^2+(y2-y1)^2+(z2-z1)^2)^(3/2));
else
    Bpx = @(x1,y1,z1,x2,y2,z2,dx2,dy2,dz2) 0; Bpy = @(x1,y1,z1,x2,y2,z2,dx2,dy2,dz2) 0; Bpz = @(x1,y1,z1,x2,y2,z2,dx2,dy2,dz2) 0;
end

NN = 1:N; %Lista de particulas en la botella
for i = 1:iter+1
    %Si alguna particula sale del limite permitido, es eliminada de la lista
    if freeze
        if any(abs(x(NN,i))>freeze_threshold*a)
            NN = setdiff(NN,NN(abs(x(NN,i))>freeze_threshold*a)); 
        elseif any(abs(y(NN,i))>freeze_threshold*a)
            NN = setdiff(NN,NN(abs(y(NN,i))>freeze_threshold*a));
        elseif any(abs(z(NN,i))>freeze_threshold*a)
            NN = setdiff(NN,NN(abs(z(NN,i))>freeze_threshold*a));
        end
        %Si quedan menos de "min_particles", se detienen los calculos
        if length(NN) < min_particles
            %Despues de un retraso
            lag_count = lag_count + 1;
            if lag_count == lag
                last_i = i %Ultima iteracion calculada
                break
            end
        end
    end
    
    for  n = 1:N
        %Si una particula no esta en la lista, su posicion se congela
        if ~any(n==N)
            x(n,i+1) = x(n,i);
            y(n,i+1) = y(n,i);
            z(n,i+1) = z(n,i);
        end
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ COORDENADA X ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
        E = [0;0;0]; %Campo electrico total en la posicion de la particula n
        Bp = [0;0;0]; %Campo magnetico debido a otras particulas en la posicion de la particula n
        
        %Se suma la interaccion de las demas particulas sobre n
        for p = NN(NN~=n) %Para los calculos en x (temporal)
            E(1) = E(1) + Ex(x(n,i),y(n,i),z(n,i),x(p,i),y(p,i),z(p,i)); 
            Bp(2) = Bp(2) + Bpy(x(n,i),y(n,i),z(n,i),x(p,i),y(p,i),z(p,i),dx(p,i),dy(p,i),dz(p,i));
            Bp(3) = Bp(3) + Bpz(x(n,i),y(n,i),z(n,i),x(p,i),y(p,i),z(p,i),dx(p,i),dy(p,i),dz(p,i));
        end
        %Campo magnetico total en la posicion de la particula n
        B = [0;
             By(x(n,i),y(n,i),z(n,i)) + Bp(2); 
             Bz(x(n,i),y(n,i),z(n,i)) + Bp(3)]; 
        
        ddx(n,i+1) = (q/m)*(q*E(1) ... %Nueva aceleracion en x (temporal)
                + dy(n,i)*B(3) ...
                - dz(n,i)*B(2) );
        dx(n,i+1) = dx(n,i) + dt*ddx(n,i+1); %Nueva aelocidad en x (temporal)
        x(n,i+1) = x(n,i) + dt*dx(n,i+1); %Nueva posicion en x (temporal)
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ COORDENADA Y ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
        E = [0;0;0]; %Campo electrico total en la posicion de la particula n
        Bp = [0;0;0]; %Campo magnetico debido a otras particulas en la posicion de la particula n
        
        %Se suma la interaccion de las demas particulas sobre n
        for p = NN(NN~=n) %Para los calculos en y, utilizando la nueva x temporal
            E(2) = E(2) + Ey(x(n,i+1),y(n,i),z(n,i),x(p,i+1),y(p,i),z(p,i));
            Bp(1) = Bp(1) + Bpx(x(n,i+1),y(n,i),z(n,i),x(p,i+1),y(p,i),z(p,i),dx(p,i+1),dy(p,i),dz(p,i));
            Bp(3) = Bp(3) + Bpz(x(n,i+1),y(n,i),z(n,i),x(p,i+1),y(p,i),z(p,i),dx(p,i+1),dy(p,i),dz(p,i));
        end
        %Campo magnetico total en la posicion de la particula n
        B = [Bx(x(n,i+1),y(n,i),z(n,i)) + Bp(1);
             0; 
             Bz(x(n,i+1),y(n,i),z(n,i)) + Bp(3)];
         
        ddy(n,i+1) = (q/m)*(q*E(2) ... %Nueva aceleracion en y
                    + dz(n,i)*B(1)... 
                    - dx(n,i+1)*B(3) );
        dy(n,i+1) = dy(n,i) + dt*ddy(n,i+1); %Nueva velocidad en y
        y(n,i+1) = y(n,i) + dt*dy(n,i+1); %Nueva posicion en y
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ COORDENADA Z ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
        E = [0;0;0]; %Campo electrico total en la posicion de la particula n
        Bp = [0;0;0]; %Campo magnetico debido a otras particulas en la posicion de la particula n
        
        %Se suma la interaccion de las demas particulas sobre n
        for p = NN(NN~=n) %Para los calculos en z, utilizando la nueva y
            E(3) = E(3) + Ez(x(n,i),y(n,i+1),z(n,i),x(p,i),y(p,i+1),z(p,i));
            Bp(1) = Bp(1) + Bpx(x(n,i),y(n,i+1),z(n,i),x(p,i),y(p,i+1),z(p,i),dx(p,i),dy(p,i+1),dz(p,i));
            Bp(2) = Bp(2) + Bpy(x(n,i),y(n,i+1),z(n,i),x(p,i),y(p,i+1),z(p,i),dx(p,i),dy(p,i+1),dz(p,i));
        end
        %Campo magnetico total en la posicion de la particula n
        B = [Bx(x(n,i),y(n,i+1),z(n,i)) + Bp(1);
             By(x(n,i),y(n,i+1),z(n,i)) + Bp(2); 
             0];
       
        ddz(n,i+1) = (q/m)*(q*E(3) ... %Nueva aceleracion en Z
                    + dx(n,i)*B(2)... 
                    - dy(n,i+1)*B(1) );
        dz(n,i+1) = dz(n,i) + dt*ddz(n,i+1); %Nueva velocidad en y
        z(n,i+1) = z(n,i) + dt*dz(n,i+1); %Nueva posicion en y
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ COORDENADA X ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
        E = [0;0;0]; %Campo electrico total en la posicion de la particula n
        Bp = [0;0;0]; %Campo magnetico debido a otras particulas en la posicion de la particula n
        for p = NN(NN~=n) %Para los calculos en x, utilizando la nueva z
            E(1) = E(1) + Ex(x(n,i),y(n,i),z(n,i+1),x(p,i),y(p,i),z(p,i+1));
            Bp(2) = Bp(2) + Bpy(x(n,i),y(n,i),z(n,i+1),x(p,i),y(p,i),z(p,i+1),dx(p,i),dy(p,i),dz(p,i+1));
            Bp(3) = Bp(3) + Bpz(x(n,i),y(n,i),z(n,i+1),x(p,i),y(p,i),z(p,i+1),dx(p,i),dy(p,i),dz(p,i+1));
        end
        %Campo magnetico total en la posicion de la particula n
        B = [0;
             By(x(n,i),y(n,i),z(n,i+1)) + Bp(2); 
             Bz(x(n,i),y(n,i),z(n,i+1)) + Bp(3)]; 
        
        ddx(n,i+1) = (q/m)*(q*E(1) ... %Nueva aceleracion en X
                    + dy(n,i)*B(3) ...
                    - dz(n,i+1)*B(2) );
        dx(n,i+1) = dx(n,i) + dt*ddx(n,i+1); %Nueva aceleracion en X
        x(n,i+1) = x(n,i) + dt*dx(n,i+1); %Nueva aceleracion en X
    end
end

%Visualizacion
if animation
    figure('Color',[0,0,0])
    BX = gca;
    box; BX.BoxStyle = 'full'; BX.Color = 'k';
    if freeze
        axis('equal',[-a-c*.4 a+c*.4 -.8*a-c*.2 .8*a+c*.2 -.8*a-c*.2 .8*a+c*.2]) 
    else
        axis('equal','vis3d')
    end
    BX.XTickLabel = [];BX.YTickLabel = [];BX.ZTickLabel = [];
    BX.XTick = [];BX.YTick = [];BX.ZTick = [];
    BX.XColor = [1,1,1];BX.YColor = [1,1,1];BX.ZColor = [1,1,1];

    L = cell(N,1);
    if N > 7
        C = jet;
        for j = 1:N
            L{j} = animatedline('LineStyle','none','MaximumNumPoint',tail,'MarkerFaceColor',C(round(length(C)*j/N),:),'Marker','o','MarkerSize',2,'MarkerEdgeColor','none'); 
        end
    else
        C = ['r';'b';'g';'y';'w';'m';'c'];
        for j = 1:N
            L{j} = animatedline('LineStyle','none','MaximumNumPoint',tail,'MarkerFaceColor',C(j),'Marker','o','MarkerSize',2,'MarkerEdgeColor','none');
        end
    end

    ROT = 20;
    view(ROT,15)
    if ~exist('last_i')
        last_i = iter+1;
    end
    for k = 1:last_i
        for kk = 1:N
           addpoints(L{kk,1},x(kk,k),y(kk,k),z(kk,k));
        end
        if mod(k,speed)==0
            %camorbit(.1,0,'data')
            ROT = ROT+rotSpeed*.1; view(ROT,15);
            CAM = campos.*1.2;campos(CAM)
            drawnow update
        end
    end 
else
   if testAxis
       subplot(2,3,1)
       plot(x)
       axis square
       subplot(2,3,2)
       plot(y)
       axis square
       subplot(2,3,3)
       plot(z)
       axis square
       subplot(2,3,[4,5,6])
       plot3(x,y,z)
   else
       for k = 1:N
           if k == 2
               hold on
           end
           plot3(x(k,:),y(k,:),z(k,:))
       end
       axis equal vis3d
   end
end
