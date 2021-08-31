function [betas] = fallossupuestos(beta,N,iter,supuesto)
%-----------------------------------------------
% PROPOSITO: Muestra los efectos de los fallos de los supuestos sobre la
%            distribucion del estimador OLS de beta, con un beta de 
%            dimension 2x1, con el elemento (1,1) siendo la constante.
%-----------------------------------------------
% INSUMOS  : beta     : 2x1 vector de parametros poblacionales
%            N        : 1x1 largo de la muestra a usar 
%            iter     : 1x1 numero de iteraciones a mostrar
%            supuesto : 1x1 : supuesto que se viola donde 
%                             0   No falla ningun supuesto
%                             1 = Linealidad
%                             2 = Exogeneidad estricta (error con media diferente de cero)
%                             3 = Exogeneidad estricta (regresores y error no ortogonales)
%                             4 = Errores esfericos (heteroscedesticidad)
%-----------------------------------------------
% OUTPUT   : Betas   : iterxK vector de parametros en todas las iteracion
%            La funcion muestra el grafico de la distribucion del parametro
%            que no es la constante.
%-----------------------------------------------

K       = size(beta,1)-1;
X       = 100*randn(N,K);          
X       = [ones(N,1), X];         
sigma   = 10;
mu      = 2;
gamma1  = 0.75;

if supuesto == 0
    betas = zeros(iter,K+1);
    for i=1:iter
        epsilon    = sigma*randn(N,1);    
        y          = X*beta + epsilon;
        beta_hat   = (X'*X)^(-1)*X'*y;
        betas(i,:) = beta_hat';
    end
elseif supuesto == 1
    betas = zeros(iter,K+1);
    for i=1:iter
        epsilon    = sigma*randn(N,1);    
        y          = (X(:,1).*beta(1,1)).*(X(:,2).*beta(2,1))  + epsilon;
        beta_hat   = (X'*X)^(-1)*X'*y;
        betas(i,:) = beta_hat';
    end
elseif supuesto == 2
    betas = zeros(iter,K+1);
    for i=1:iter
        epsilon    = mu + sigma*randn(N,1);    
        y          = X*beta + epsilon;
        beta_hat   = (X'*X)^(-1)*X'*y;
        betas(i,:) = beta_hat';
    end
elseif supuesto == 3
    betas = zeros(iter,K+1);
    for i=1:iter
        M       = mu + sigma*randn(N,2);
        R       = [1 gamma1; gamma1 1];
        L       = chol(R);
        M       = M*L;
        X(:,2)  = M(:,1); 
        epsilon = M(:,2);
        y          = X*beta + epsilon;
        beta_hat   = (X'*X)^(-1)*X'*y;
        betas(i,:) = beta_hat';
    end
elseif supuesto == 4
    betas = zeros(iter,K+1);
    MU = zeros(N,1);subplot(2,2,2)
hist(betas(:,2),80), hold on
plot(beta(2), betas(:,2),'r'), ylim([0 0.05*iter]),  title({'Distribución \beta_1'}),  hold off
sgtitle('Distribución de los Parametros Estimados por MCO') 

    V  = eye(N);
    for j=1:N
        V(j,j) = j;
    end 
    for i=1:iter
        epsilon    =mvnrnd(MU,V,1)'; 
        y          = X*beta + epsilon;
        beta_hat   = (X'*X)^(-1)*X'*y;
        betas(i,:) = beta_hat';
    end
else
     error('supuesto solo puedo tomar los valores {0,1,2,3,4}')
end 


figure(1)
subplot(2,2,1)
hist(betas(:,1),80), hold on
plot(beta(1), betas(:,1),'r'), ylim([0 0.05*iter]), title({'Distribución \beta_0'}), hold off


subplot(2,2,2)
hist(betas(:,2),80), hold on
plot(beta(2), betas(:,2),'r'), ylim([0 0.05*iter]),  title({'Distribución \beta_1'}),  hold off
suptitle('Distribución de los Parámetros Estimados por MCO') 


end 