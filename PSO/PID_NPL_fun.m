function wskaznik_jakosci = PID_NPL_fun(parameters)


global alfa1 alfa2 beta1 beta2;
global pso_scale;
%rzad
na = 2;
nb = 5;

%czas trwania symulacji
global kk;

%wybor algorytmu
global algorytm;
%algorytm='PID'; %GPC NPL PID oraz NO, ktory nie dziala

if strcmp(algorytm, 'PID')
    %parametry ciaglego PID
    T = 1;
    K = parameters(1);
    Ti = parameters(2);
    Td = parameters(3);
    %parametry dyskretnego PID
    r0 = K*(1+(T/(2*Ti)) + Td/T);
    r1 = K*((T/(2*Ti)) - ((2*Td)/T) - 1);
    r2 = (K*Td)/T;
elseif strcmp(algorytm, 'NPL')
    N = round(pso_scale*parameters(1));
    Nu = round(pso_scale*parameters(2));
    lambda = parameters(3);
    if N<=2 || Nu<=2 || lambda <=0
        wskaznik_jakosci = 10^3;
        return
    end
end


%uchyb PID
e = zeros(1,kk);

%parametry regulatora GPC i NPL
%global N;
% N = 15;
% %global Nu;
% Nu = 5;
% %global lambda;
% lambda = 1;

%inicjalizacja
global u; u = zeros(1,kk);        %sterowanie
global y; y = zeros(1,kk);        %wyjscie obiektu
%wartosc zadana
global yzad;

ymod     = zeros(1,kk); %wyjscie modelu
ymod_lin = zeros(1,kk); %wyjscie modelu zlinearyzowanego

x1 = zeros(1,kk);       %x1
x2 = zeros(1,kk);       %x2

if strcmp(algorytm, 'NPL')
    y_zero   = zeros(N,1);  %trajektoria swobodna
    s  = zeros(1,N);        %odpowiedz skokowa
    I  = eye(Nu);
    Yzad=yzad.*ones(N,1);
end

%ograniczenia
u_min = -1;
u_max = 1;
% umin_opt = u_min*ones(Nu,1);
% umax_opt = u_max*ones(Nu,1);


global w10 w1 w20 w2
delta = 10^(-7);    %delta do linearyzacji modelu


%%%--- GPC - obliczenia offline ---%%%
if strcmp(algorytm, 'GPC')
    %%%--- model liniowy ---%%%
    w(1)=-0.0604346804947776;
    w(2)=0.0484656613758900;
    w(3)=1.89416871378879;
    w(4)=-0.9007809402858;
    b = [0 0 0 w(1) w(2)];
    a = [w(3) w(4)];
    %%%--- GPC - odpowiedz skokowa ---%%%
    for j=1:N
        %suma b
        suma_b = 0;
        i_max  = min([j, nb]);
        for i=1:i_max
            suma_b = b(i) + suma_b;
        end
        
        %suma a
        suma_a = 0;
        i_max  = min([j-1, na]);
        for i=1:i_max
            suma_a = a(i)*s(j-i) + suma_a;
        end
        s(j) = suma_b - suma_a;
    end
    
    %%%--- GPC - macierz dynamiczna ---%%%
    M=zeros(N, Nu);
    i=0;
    for j=1:Nu
        M(j:N,j)=s(1:N-i).';
        i=i+1;
    end
    K=inv(M.'*M+lambda*I)*M.';
end

%%%%%%%%%%%%%%%%-------------- PETLA REGULACJI --------------%%%%%%%%%%%%%%%%
%global k;
%global d;
for k=6:kk
    %%%---pomiar wyjscia obiektu---%%%
    x1(k) = -alfa1*x1(k-1)+x2(k-1)+beta1*g1(u(k-4));
    x2(k) = -alfa2*x1(k-1)+beta2*g1(u(k-4));
    y(k)  = g2(x1(k));
    
    %uchyb
    e(k) = yzad(k)-y(k);
    
    if strcmp(algorytm, 'PID')
        
        
        %sygnal sterujacy regulatora PID
        u(k) = r2*e(k-2) + r1*e(k-1) + r0*e(k) + u(k-1);
    end
    
    if strcmp(algorytm, 'NPL') || strcmp(algorytm, 'GPC') || strcmp(algorytm, 'NO')
        %%%---wyjscie modelu---%%%
        if strcmp(algorytm, 'NPL') || strcmp(algorytm, 'NO')
            wesn_arx = [u(k-4) u(k-5) y(k-1) y(k-2)]';
            ymod(k)  = w20+w2*tanh(w10+w1*wesn_arx);
        elseif strcmp(algorytm, 'GPC')
            ymod(k)  = b(4)*u(k-4) + b(5)*u(k-5) + a(1)*y(k-1) + a(2)*y(k-2);
        end
        
        %%%---pomiar zaklocenia---%%%
        d=y(k)-ymod(k);
        
        if strcmp(algorytm, 'NO')
            %             u0 = u(k-1)*ones(1,Nu);
            %             opcje = optimset('Algorithm','sqp', 'Display', 'off', 'TolFun',1e-10,'TolX',1e-10);
            %             uopt = fmincon(@fun, u0, [],[],[],[], umin_opt, umax_opt, [], opcje);
            %             drawnow;
            %             u(k)=uopt(1);
        elseif strcmp(algorytm, 'NPL')
            %%%---trajektoria swobodna---%%%
            wesn = [u(k-3) u(k-4) y(k) y(k-1)]';
            y_zero(1) = w20+w2*tanh(w10+w1*wesn) + d;
            if N>=2
            wesn = [u(k-2) u(k-3) y_zero(1) y(k)]';
            y_zero(2) = w20+w2*tanh(w10+w1*wesn) + d;
            end
            if N>=3
            wesn = [u(k-1) u(k-2) y_zero(2) y_zero(1)]';
            y_zero(3) = w20+w2*tanh(w10+w1*wesn) + d;
            end
            if N>=4
            for i=4:N
                wesn = [u(k-1) u(k-1) y_zero(i-1) y_zero(i-2)]';
                y_zero(i) = w20+w2*tanh(w10+w1*wesn) + d;
            end
            end
            
            %%%---linearyzacja modelu---%%%
            y_pp = ymod(k);
            
            %wspolczynnik b4
            wesn = [u(k-4)+delta u(k-5) y(k-1) y(k-2)]';
            y_d = w20+w2*tanh(w10+w1*wesn);
            b4 = (y_d-y_pp)/delta;
            
            %wspolczynnik b5
            wesn = [u(k-4) u(k-5)+delta y(k-1) y(k-2)]';
            y_d = w20+w2*tanh(w10+w1*wesn);
            b5  = (y_d-y_pp)/delta;
            
            %wspolczynnik a1
            wesn = [u(k-4) u(k-5) y(k-1)+delta y(k-2)]';
            y_d = w20+w2*tanh(w10+w1*wesn);
            a1  = -(y_d-y_pp)/delta;
            
            %wspolczynnik a2
            wesn = [u(k-4) u(k-5) y(k-1) y(k-2)+delta]';
            y_d = w20+w2*tanh(w10+w1*wesn);
            a2  = -(y_d-y_pp)/delta;
            
            %model zlinearyzwoany
            ymod_lin(k) = b4*u(k-4) + b5*u(k-5) + a1*y(k-1) + a2*y(k-2);
            
            %%%---odpowiedz skokowa---%%%
            
            a = [a1 a2];
            b = [0 0 0 b4 b5];
            for j=1:N
                %suma b
                suma_b = 0;
                i_max  = min([j, nb]);
                for i=1:i_max
                    suma_b = b(i) + suma_b;
                end
                
                %suma a
                suma_a = 0;
                i_max  = min([j-1, na]);
                for i=1:i_max
                    suma_a = a(i)*s(j-i) + suma_a; %tu mo�e byc blad
                end
                s(j) = suma_b - suma_a;
            end
            
            %%%---macierz dynamiczna---%%%
            M=zeros(N, Nu);
            i=0;
            for j=1:Nu  %wypelnianie macierzy trojkatnej dolnej M
                M(j:N,j)=s(1:N-i)';
                i=i+1;
            end
            K=inv(M.'*M+lambda*I)*M';
        elseif strcmp(algorytm, 'GPC')
            %%%-----------GPC------------%%%
            %%%---trajektoria swobodna---%%%
            y_zero(1) = b(4)*u(k-3) + b(5)*u(k-4) + a(1)*y(k) + a(2)*y(k-1) + d;
            if N>=2
                y_zero(2) = b(4)*u(k-2) + b(5)*u(k-3) + a(1)*y_zero(1) + a(2)*y(k) + d;
            end
            if N>=3
                y_zero(3) = b(4)*u(k-1) + b(5)*u(k-2) + a(1)*y_zero(2) + a(2)*y_zero(1) + d;
            end
            if N>=4
                for i=4:N
                    y_zero(i) = b(4)*u(k-1) + b(5)*u(k-2) + a(1)*y_zero(i-1) + a(2)*y_zero(i-2) + d;
                end
            end
        end
        %%%---obliczenie sterowania---%%%
        if size(Yzad(:,k),1)~=size(y_zero,1)
            disp('Pierwszy wymiar')
            fprintf('N=%d, Nu=%f, lambda=%f',N, Nu, lambda)
            
            disp('Yzad'); disp(size(Yzad));
            disp('y_zero'); disp(size(y_zero));
        end
%         if size(Yzad(:,k),2)~=size(y_zero,2)
%             disp('Drugi wymiar')
%             fprintf('N=%d, Nu=%f, lambda=%f',N, Nu, lambda)
%             
%             disp('Yzad'); disp(size(Yzad));
%             disp('y_zero'); disp(size(y_zero));
%         end
        if size(K,2)~=size((Yzad(:,k)-y_zero),1)
            fprintf('N=%d, Nu=%f, lambda=%f',N, Nu, lambda)
            disp('K'); disp(K);
            disp('Yzad'); disp(Yzad);
            disp('y_zero'); disp(y_zero);
            %fprintf('=%d, Nu=%f, lambda=%f',N, Nu, lambda)
        end
        deltaU=K*(Yzad(:,k)-y_zero);
        if size(deltaU,1) < 1
            fprintf('N=%d, Nu=%f, lambda=%f',N, Nu, lambda)
            disp(deltaU)
            
        end
        delta_u=deltaU(1);
        u(k)=u(k-1)+delta_u;
    end
    
    %ograniczenie sygna�u steruj�cego
    if not(strcmp(algorytm, 'NO'))
        if u(k)>u_max
            u(k)=u_max;
        elseif u(k)<u_min
            u(k)=u_min;
        end
    end
end
wskaznik_jakosci=sum(e.^2);
%save('temp_to_del.mat','U','Y','yzad')
end