function [output]= MPC_Sam_2obj(input_values)
x=input_values(1:3); % état x(k)
g=input_values(4,1); % Ptolhaut: contrainte et objectif
ref_SoC=30; %Référence pour le SoC
global Np Ad Bd Cc 
% Calculation of Phi_i, Pi and Psi_i H2, h1, A_ineg and B_ineg
m=size(Bd,2); %Quantité d'entrées
nobj=2; %Nombre de références
Phi_i=zeros(Np*size(Ad,1),size(Ad,1));
Phi_i(1:size(Ad,1),1:size(Ad,2))=Ad;
Psi_i=zeros(Np*size(Ad,1),Np);
psi_i=Bd;
Psi_i(1:size(Ad,1),:)=psi_i*[eye(m), zeros(m,Np-1)];
    for k=1:Np+1
        Phi_i(k*size(Ad,1)+1:(k+1)*size(Ag,1),:)=Phi_i((k-1)*size(Ad,1)+1:k*size(Ad,1),:)*Ad;            % Phi_i=A^i
    end
% Phi_i=[A; A*A; A*A*A ...]
    for k=1:Np-1
        psi_i=[Ad*psi_i(1:size(Ad,1),m),psi_i]; %psi_i=[A^(i-1)B, ... AB, B]
        Psi_i((k)*size(Ad,1)+1:(k+1)*size(Ad,1),:)=psi_i*[eye((k+1)*m) zeros((k+1)*m,(Np-k-1)*m)];  % [A^(i-1)*B, ... ,A*B,B]*Pi(i)
    end
ncons=3; % 3 contraintes: Ptolhaut, Socmin et Socmax
A_ineg=zeros(Np*ncons,Np);  
    for k=1:Np
        A_ineg((k-1)*ncons+1:k*ncons,:)=[Cc(1,:); Cc(2,:); -Cc(2,:)]*Psi_i((k-1)*size(Ad,1)+1:k*size(Ad,1),:);
        %A_ineg((k-1)*ncons+1:k*ncons,:)=[Cc(2,:); -Cc(2,:)]*Psi_i((k-1)*size(Ad,1)+1:k*size(Ad,1),:);
        %A_ineg((k-1)*ncons+1:k*ncons,:)=[Cc(1,:); Cc(2,:)]*Psi_i((k-1)*size(Ad,1)+1:k*size(Ad,1),:);
        % Première ligne: Pbat, deuxième ligne SoC, troisième ligne -SOC.
    end
B_ineg=zeros(ncons*Np,1); 
gamma1_max=1*g; % Ptol_haut
SoC_max=80; % SoC max
SoC_min=20; % Pinj and SoC
    for k=1:Np
      B_ineg((k-1)*ncons+1:k*ncons,:)=[gamma1_max; SoC_max; -SoC_min]-[Cc(1,:); Cc(2,:); -Cc(2,:)]*Phi_i((k-1)*size(Ad,1)+1:k*size(Ad,1),:)*x; %SoC and Pg
      %B_ineg((k-1)*ncons+1:k*ncons,:)=[SoC_max; -SoC_min]-[Cc(2,:); -Cc(2,:)]*Phi_i((k-1)*size(Ad,1)+1:k*size(Ad,1),:)*x;
    end
Q=diag([1 1]); R=g;
H2=0;
F1=zeros(Np,size(Ad,1));
F2=zeros(Np,Np*nobj);
Pi=zeros(1,Np);
Pi_yr=zeros(nobj,Np*nobj);
    for k=1:Np
        omega_i1=Cc*Psi_i((k-1)*size(Ad,1)+1:k*size(Ad,1),:);
        Pi(1,k)=1;
        Pi_yr(:,k+1:k+nobj)=eye(nobj); 
        H2=H2+omega_i1'*Q*omega_i1+Pi'*R*Pi; %
        C_phi_i=Cc*Phi_i((k-1)*size(Ad,1)+1:k*size(Ad,1),:);
        F1=F1+omega_i1'*Q*(C_phi_i);
        F2=F2+omega_i1'*Q*Pi_yr;
    end
% Solution of the problem min 0.5*x’*H*x + f’*x subject to: A*x <= b:   
Yref=repmat(eye(nobj),Np,1)*[g; ref_SoC]; %Yref=[g, ..., g];           % modifier
h1=F1*x-F2*Yref;
[X,fval]=quadprog(H2/2,h1,A_ineg,B_ineg);
output=[X(1,1),fval];
end