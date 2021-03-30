clc
clear
close all
%% Dynamic unicycle agent plant
st=1;
d=0.2;
vt=-1:st:1;
Bg=[0 0
    0 0
    1 0
    0 d];
Cg=[1 0 0 0
    0 1 0 0];
Dg=[0 0
    0 0];
for i=1:length(vt)
    Ag(:,:,i)=[0 vt(i)/d 1 0
               -vt(i)/d 0 0 1
               0 0 0 0
               0 0 0 0];
    CM(:,:,i)=ss(Ag(:,:,i),Bg,Cg,Dg);
end

umax = [50 15];
rmax = [10 10];
Du = diag(umax);
Dr = diag(rmax);
SCM = Dr\CM*Du; % Scaled Continous Model
Ag=SCM.A;
Bg=SCM.B(:,:,1);
Cg=SCM.C(:,:,1);
Dg=SCM.D(:,:,1);
nxg=size(Ag,1);
nu=size(Bg,2);
nw=size(Cg,1);
%% Shaping filters
WS1=tf([0 0.5e1],[1 0.5e-2]);
WS  = ss(mdiag(WS1,WS1));
As=WS.A;
Bs=WS.B;
Cs=WS.C;
Ds=WS.D;
nxs=size(As,1);
nzs=size(Cs,1);
nus=size(Ds,2);
%
WK1=tf([1 100],[1 100000]);
WK  = ss(mdiag(WK1,WK1));
Ak=WK.A;
Bk=WK.B;
Ck=WK.C;
Dk=WK.D;
nxk=size(Ak,1);
nzk=size(Ck,1);
nuk=size(Dk,2);
Wsf=tf([0 1e-2],[1 1e-4]); 
WSF  = ss(mdiag(Wsf,Wsf));                   
systemnames = 'WSF';
inputvar = '[w{2};r{2};u{2}]';
outputvar = '[r-u;WSF;w+r-u]';
input_to_WSF = '[w+r-u]';
cleanupsysic = 'yes';
PF = sysic;
[F,CL,GAM,INFO] = hinfsyn(PF,2,2,'method','lmi');
Af=F.A;
Bf=F.B;
Cf=F.C;
Df=F.D;
%% Generalized Plant
for i=1:length(vt)
    A_rho(:,:,i)=[Ag(:,:,i) zeros(nxg,nxs) zeros(nxg,nxk)
                  -Bs*Cg As zeros(nxs,nxk)
                  zeros(nxk,nxg) zeros(nxk,nxs) Ak];
end
Bp=[zeros(nxg,nw)
    Bs
    zeros(nxk,nw)];
Bu=[Bg
    zeros(nxs,nu)
    Bk];
Cp=[-Ds*Cg Cs zeros(nzs,nxk)
    zeros(nzk,nxg) zeros(nzk,nxs) Ck];
Cy=[-Cg zeros(nw,nxs) zeros(nw,nxk)];
Dpp=[Ds
     zeros(nzk,nus)];
Dpu=[zeros(nzs,nuk)
     Dk];
Dyp=eye(nw);
Dyu=zeros(nw,nu);
%% LFT representation matrices
A_=[0 0 1 0
    0 0 0 1
    0 0 0 0
    0 0 0 0];
C_=[0 1/d 0 0
    -1/d 0 0 0];
A=[A_ zeros(nxg,nxs) zeros(nxg,nxk)
   -Bs*Cg As zeros(nxs,nxk)
   zeros(nxk,nxg) zeros(nxk,nxs) Ak];
Bd=[[eye(2);zeros(nxg-2,nu)]
    zeros(nxs,2)
    zeros(nxk,2)];

Cd=[C_ zeros(2,nxs) zeros(2,nxk)];

Ddd=zeros(2,2);
Ddp=zeros(2,nw);
Ddu=zeros(2,nu);

Dpd=zeros(nzs+nzk,2);

Dyd=zeros(nw,2);
%% Matrices preparation for existance condition
G_R11=Ddd;
G_R12=[Cd Ddp];
G_R21=[zeros(size(A,1),size(Bd,2))
       Bd
       zeros(size(Bp,2),size(Bd,2))
       Dpd];
G_R22=[eye(size(A,1)) zeros(size(A,1),size(Bp,2))
       A Bp
       zeros(size(Bp,2),size(A,2)) eye(size(Bp,2))
       Cp Dpp];
GR=[G_R11 G_R12
    eye(size(G_R11,2)) zeros(size(G_R11,2),size(G_R12,2))
    G_R21 G_R22];
%
G_S11=-Ddd';
G_S12=[-Bd' -Dpd'];
G_S21=[-Cd'
       zeros(size(A,1),size(Cd',2))
       -Ddp'
       zeros(size(Dpp',2),size(Ddp',2))];
G_S22=[-A' -Cp'
       eye(size(A',1)) zeros(size(A',1),size(Cp',2))
       -Bp' -Dpp'
       zeros(size(Dpp',2),size(Bp',2)) eye(size(Dpp',2))];
GS=[eye(size(G_S11,2)) zeros(size(G_S11,2),size(G_S12,2))
    G_S11 G_S12
    G_S21 G_S22];
%
N_R=null([Dyd Cy Dyp],'r');
N_S=null([Ddu' Bu' Dpu'],'r');
%% Yalmip for solving for existance condition
% Yalimp Parameters
Zgap = 1e-100;
clear yalmip
LMIsolver = 'sdpt3';
SDPoptions = sdpsettings('savesolveroutput',0,'verbose',1,'debug', 1,'solver',LMIsolver);
SDPoptions.sdpt3.maxit = 1500;
SDPoptions.sdpt3.steptol= 1e-8;

R=sdpvar(nxg+nxs+nxk,nxg+nxs+nxk,'symmetric');
S=sdpvar(nxg+nxs+nxk,nxg+nxs+nxk,'symmetric');
gm= sdpvar(1,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D/G scalings
M11=sdpvar(2,2,'symmetric');
M22=sdpvar(2,2,'symmetric');
M12=sdpvar(2,2,'skew');
M=[M11 M12;
   M12' -1*M11];
N11=sdpvar(2,2,'symmetric');
N22=sdpvar(2,2,'symmetric');
N12=sdpvar(2,2,'skew');
N=[N11 N12;
   N12' -1*N11];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LRM=N_R'*GR'*[M zeros(4,nxg+nxs+nxk) zeros(4,nxg+nxs+nxk) zeros(4,nw+nzs+nzk)
              zeros(nxg+nxs+nxk,4) zeros(nxg+nxs+nxk) R zeros(nxg+nxs+nxk,nw+nzs+nzk)
              zeros(nxg+nxs+nxk,4) R zeros(nxg+nxs+nxk) zeros(nxg+nxs+nxk,nw+nzs+nzk)
              zeros(nw+nzs+nzk,4) zeros(nw+nzs+nzk,nxg+nxs+nxk) zeros(nw+nzs+nzk,nxg+nxs+nxk) [-gm*eye(nw) zeros(nw,nzs+nzk);zeros(nzs+nzk,nw) (1/gm)*eye(nzs+nzk)]]*GR*N_R;
Y=GR*N_R;
Wy1=Y(1:(4+2*(nxg+nxs+nxk)),:);
Wy2=Y((4+2*(nxg+nxs+nxk))+1:(4+2*(nxg+nxs+nxk))+2,:);
Wy3=Y((4+2*(nxg+nxs+nxk))+3:end,:);
Wy0=[M zeros(4,nxg+nxs+nxk) zeros(4,nxg+nxs+nxk)
    zeros(nxg+nxs+nxk,4) zeros(nxg+nxs+nxk) R
    zeros(nxg+nxs+nxk,4) R zeros(nxg+nxs+nxk)];
LRMlmi=[Wy1'*Wy0*Wy1-gm*(Wy2'*Wy2) Wy3'
        Wy3 -gm*eye(nzs+nzk)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
LSN=N_S'*GS'*[N zeros(4,nxg+nxs+nxk) zeros(4,nxg+nxs+nxk) zeros(4,nw+nzs+nzk)
              zeros(nxg+nxs+nxk,4) zeros(nxg+nxs+nxk) S zeros(nxg+nxs+nxk,nw+nzs+nzk)
              zeros(nxg+nxs+nxk,4) S zeros(nxg+nxs+nxk) zeros(nxg+nxs+nxk,nw+nzs+nzk)
              zeros(nw+nzs+nzk,4) zeros(nw+nzs+nzk,nxg+nxs+nxk) zeros(nw+nzs+nzk,nxg+nxs+nxk) [(-1/gm)*eye(nw) zeros(nw,nzs+nzk);zeros(nzs+nzk,nw) gm*eye(nzs+nzk)]]*GS*N_S;
X=GS*N_S;
W1=X(1:(4+2*(nxg+nxs+nxk)),:);
W2=X((4+2*(nxg+nxs+nxk))+1:(4+2*(nxg+nxs+nxk))+2,:);
W3=X((4+2*(nxg+nxs+nxk))+3:end,:);
W0=[N zeros(4,nxg+nxs+nxk) zeros(4,nxg+nxs+nxk)
    zeros(nxg+nxs+nxk,4) zeros(nxg+nxs+nxk) S
    zeros(nxg+nxs+nxk,4) S zeros(nxg+nxs+nxk)];

LSNlmi=[W1'*W0*W1+gm*(W3'*W3) W2'
        W2 gm*eye(nw)];    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
LRS=[R eye(nxg+nxs+nxk)
     eye(nxg+nxs+nxk) S];

Constraints =  LRMlmi <= -Zgap*eye(size(LRMlmi));
Constraints = [Constraints , LSNlmi >= Zgap*eye(size(LSNlmi))];
Constraints = [Constraints , LRS >= Zgap*eye(size(LRS))];

Constraints = [Constraints ,M11 >= Zgap*eye(size(M11))];
Constraints = [Constraints ,N11 >= Zgap*eye(size(N11))];
Constraints = [Constraints ,(gm <= 1)];
Constraints

optimize(Constraints,[],SDPoptions)
gam = value(gm)
M=value(M);
N=value(N);
S=value(S);
R=value(R);
%% Controller Construction based on p.59-60
MNT=eye(nxg+nxs+nxk)-S*R;
[u,s,v]=svd(MNT);
Mx=u*s;
Nx=v;
F=-inv(Dpu'*Dpu)*(gam*Bu'*inv(S)+Dpu'*Cp);
L=-(gam*inv(R)*Cy'+Bp*Dyp')*inv(Dyp*Dyp');
for i=1:length(vt)
    A_k(:,:,i)=-inv(Nx)*(A_rho(:,:,i)'+R*(A_rho(:,:,i)+Bu*F+L*Cy)*S+inv(gam)*R*(Bp+L*Dyp)*Bp'+inv(gam)*Cp'*(Cp+Dpu*F)*S)*(inv(Mx))';
end
B_k=inv(Nx)*R*L;
C_k=F*S*(inv(Mx))';
D_k=zeros(size(C_k,1),size(B_k,2));
Bg=[0 0
    0 0
    1 0
    0 d];
Cg=[1 0 0 0
    0 1 0 0];
Dg=[0 0
    0 0];
for i=1:length(vt)
    Ag(:,:,i)=[0 vt(i)/d 1 0
               -vt(i)/d 0 0 1
               0 0 0 0
               0 0 0 0];
    G(:,:,i)=ss(Ag(:,:,i),Bg,Cg,Dg);
    kd(:,:,i)=ss(A_k(:,:,i),B_k,C_k,D_k);
end
Kd=Du*kd*inv(Dr);

%% Simulate in simulink
sim('IFF_Heterogenous_simu')
%% Plot
figure(5)
for k=1:1000:length(out.tout)
    plot(out.ref.Data(:,1),out.ref.Data(:,2),'c','LineWidth',2)
    hold on
    for nn=1:7
        if nn==1
            line([out.XY.Data(1,nn,k)-d*cos(out.The.Data(1,nn,k))-0.15*cos(out.The.Data(1,nn,k)) out.XY.Data(1,nn,k)],[out.XY.Data(2,nn,k)-d*sin(out.The.Data(1,nn,k))-0.15*sin(out.The.Data(1,nn,k)) out.XY.Data(2,nn,k)],'Color','b','LineWidth',5)
            hold on
            plot(out.XY.Data(1,nn,k)-d*cos(out.The.Data(1,nn,k)),out.XY.Data(2,nn,k)-d*sin(out.The.Data(1,nn,k)),'.y','Markersize',10)
            hold on
            plot(out.XY.Data(1,nn,k),out.XY.Data(2,nn,k),'.r','Markersize',10)
            hold on
        else
            xxx=out.XY.Data(1,nn,1:k);
            yyy=out.XY.Data(2,nn,1:k);
            plot(xxx(:),yyy(:),'--c')
            hold on
            line([out.XY.Data(1,nn,k)-d*cos(out.The.Data(1,nn,k))-0.15*cos(out.The.Data(1,nn,k)) out.XY.Data(1,nn,k)],[out.XY.Data(2,nn,k)-d*sin(out.The.Data(1,nn,k))-0.15*sin(out.The.Data(1,nn,k)) out.XY.Data(2,nn,k)],'Color','k','LineWidth',5)
            hold on
            plot(out.XY.Data(1,nn,k)-d*cos(out.The.Data(1,nn,k)),out.XY.Data(2,nn,k)-d*sin(out.The.Data(1,nn,k)),'.y','Markersize',10)
            hold on
            plot(out.XY.Data(1,nn,k),out.XY.Data(2,nn,k),'.r','Markersize',10)
            hold on
        end
    xlim([-6 12])
    ylim([-7 7])
    xlabel('X [m]')
    ylabel('Y [m]')
    grid on
    end
    hold off
    pause(0.01)
end
%%
figure(6)
for nn=1:7
        ttt=out.The.Data(1,nn,1:1000:length(out.tout));
        plot(1e-3.*[1:1000:length(out.tout)],ttt(:),'b')
        hold on
end
grid on
xlim([0 100])
xlabel('T [sec]')
ylim([-6 2])
ylabel('\phi [rad]')
%%
figure(7)
xxxxx=out.efa.Data(1:1000:length(out.tout));
 plot(1e-3.*[1:1000:length(out.tout)],xxxxx(:),'b','Linewidth',2)
grid on
xlim([0 100])
xlabel('T [sec]')
ylim([0 4])
ylabel('$\|e_f\|_{{L}_2}~[rad]$','interpreter','latex')