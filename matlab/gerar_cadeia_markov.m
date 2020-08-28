% Experiência com fonte em árvore - exemplo do artigo do willems
clear
theta1=0.1;
theta10=0.3;
theta00=0.5;

% Cálculo das probabilidades estacionárias considerando a matriz de
% transição

pi1=theta00/(theta00+(1-theta1)*theta00+(1-theta1)*(1-theta10));
pi10=pi1*(1-theta1);
pi00=pi1*(1-theta1)*(1-theta10)/theta00;

%% Calculando DI analítica

clc
p=0.1; % BSC com parâmetro p
pcond(1)=(1-p)*(1-theta00)/((1-p)*(1-theta00)+p*theta00);
pcond(2)=(1-p)*(1-theta10)/((1-p)*(1-theta10)+p*theta10);
pcond(3)=(1-p)*(1-theta1)/((1-p)*(1-theta1)+p*theta1);
pcond(4)=(p)*(1-theta00)/((p)*(1-theta00)+(1-p)*theta00);
pcond(5)=(p)*(1-theta10)/((p)*(1-theta10)+(1-p)*theta10);
pcond(6)=(p)*(1-theta1)/((p)*(1-theta1)+(1-p)*theta1);
pcond(7)=(p)*(theta00)/((p*theta00)+(1-p)*(1-theta00));
pcond(8)=(p)*(theta10)/((p*theta10)+(1-p)*(1-theta10));
pcond(9)=(p)*(theta1)/((p*theta1)+(1-p)*(1-theta1));
pcond(10)=(1-p)*(theta00)/((1-p)*(theta00)+p*(1-theta00));
pcond(11)=(1-p)*(theta10)/((1-p)*(theta10)+p*(1-theta10));
pcond(12)=(1-p)*(theta1)/((1-p)*(theta1)+p*(1-theta1));

pconj(1)=(1-theta00)*pi00*(1-p)+theta00*pi00*p;
pconj(2)=(1-theta10)*pi10*(1-p)+theta10*pi10*p;
pconj(3)=(1-theta1)*pi1*(1-p)+theta1*pi1*p;
pconj(4)=(1-theta00)*pi00*p+theta00*pi00*(1-p);
pconj(5)=(1-theta10)*pi10*p+theta10*pi10*(1-p);
pconj(6)=(1-theta1)*pi1*p+theta1*pi1*(1-p);

h_causal_cond= -pconj(1)*(pcond(1)*log2(pcond(1))+pcond(7)*log2(pcond(7)))...
    -pconj(2)*(pcond(2)*log2(pcond(2))+pcond(8)*log2(pcond(8)))...
    -pconj(3)*(pcond(3)*log2(pcond(3))+pcond(9)*log2(pcond(9)))...
    -pconj(4)*(pcond(4)*log2(pcond(4))+pcond(10)*log2(pcond(10)))...
    -pconj(5)*(pcond(5)*log2(pcond(5))+pcond(11)*log2(pcond(11)))...
    -pconj(6)*(pcond(6)*log2(pcond(6))+pcond(12)*log2(pcond(12)));
h_X=-pi00*(theta00*log2(theta00)+(1-theta00)*log2(1-theta00))...
    -pi10*(theta10*log2(theta10)+(1-theta10)*log2(1-theta10))...
    -pi1*(theta1*log2(theta1)+(1-theta1)*log2(1-theta1));

di_analytical=h_X-h_causal_cond;

%% Gerando dados

N=100000; %tamanho da sequência
theta1=0.1;
theta10=0.3;
theta00=0.5;
D=2;
DIrate=zeros(50,N-D);


for trial=1:50
    
    fprintf('\n generating data, trial %d \n', trial)
    tic
    x= ones(N,1);
    y=x;
        
    for i=2:length(x)        
        if x(i-1)==1
            x(i)=binornd(1,theta1,1);
        else if x(i-2)==0
                x(i)=binornd(1,theta00,1);
            else x(i)=binornd(1,theta10,1);
            end
        end                
    end
    
    for i=1:N
        if x(i)==0
            y(i)=rand(1,1) <= p;
        else
            y(i)=rand(1,1) <=(1-p);
        end
    end
    toc
    di(trial) = DI_plugin_estimate(x,y,D)
end

%% Chi Squared

estimation = 2*N.*di;
l = 2; %tamanho do alfabeto de X
m = 2; %tamanho do alfabeto de Y
k = D; %tamanho da memoria
dof = (l^k)*(m^(k+1)-1)*(l-1); % graus de liberdade
percentil = chi2inv(0.95,dof);
if(estimation(1) > percentil)
    disp('Podemos rejeitar a hipótese nula - Há causalidade!')
else
    disp('Não há causalidade.')
end
figure(1)
histfit(di', 10, 'Normal')
title('Histograma DI_n com \mu = 0.434 e \sigma = 0.002')
xlabel('Informação Direcional')
ylabel('# de Ocorrências')
