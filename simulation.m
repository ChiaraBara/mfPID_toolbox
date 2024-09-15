clear all; clc; close all;

%%%% simulation of Boolean logic AND-gate with 2 continuous source variables

addpath([pwd,'\functions\']);

a = -0.5; b = -a;
d = [0.05:0.05:0.95]; c = d-1;

N = 300;
k = 5;
iy = 1; ix1 = 2; ix2 = 3;

nsurr = 100;
alpha = 0.05;

%%%% theoretical values
states = [0 0 0; 1 0 0; 0 1 0; 1 1 0; 0 0 1; 1 0 1; 0 1 1; 1 1 1]; % [Y X1 X2]
for i = 1:length(c)
    joint_prob = [(a*c(i))/((b-a)*(d(i)-c(i))); 0; (-b*c(i))/((b-a)*(d(i)-c(i))); 0;...
        (-a*d(i))/((b-a)*(d(i)-c(i))); 0; 0; (b*d(i))/((b-a)*(d(i)-c(i)))];

    outth = mfPID_2sources_th(joint_prob,states);
    Ith(i) = outth.I;
    I1th(i) = outth.I1;
    I2th(i) = outth.I2;
    U1th(i) = outth.U1;
    U2th(i) = outth.U2;
    Sth(i) = outth.S;
    Rth(i) = outth.R;
    pYth(i,:) = outth.pY;
    Ispecth(i,:,:) = outth.Ispec;
end

%%% realization
X1 = a + (b-a)*rand(N,1);
for i = 1:length(c)
    X2 = c(i) + (d(i)-c(i)).*rand(N,1);
    Y = Heaviside(X1).*Heaviside(X2); 
    Zn(i,:,:) = [Y,X1,X2];
end

for i = 1:length(c) 

    Zn_tmp = squeeze(Zn(i,:,:));

    %%% estimated value
    out = mfPID_2sources_mixed_mex(Zn_tmp,iy,ix1,ix2,k,1);
    I(i) = out.I;
    U1(i) = out.U1;
    U2(i) = out.U2;
    S(i) = out.S;
    R(i) = out.R;
    pY(i,:) = out.pY;
    Dspec(i,:,:) = out.Dspec;

    %%% significance test
    Q = size(pY,2);
    th = alpha/(Q*4);
    pos_pY0 = ones(1,Q); pos_pY0(pY(i,:)==0) = 0;
    for isurr = 1:nsurr
        Zn_surr = Zn_tmp;
        Zn_surr(:,1) = surrshuf(Zn_tmp(:,1));
        out_surr = mfPID_2sources_mixed_mex(Zn_surr,iy,ix1,ix2,k);
        Dspec_surr(isurr,:,:) = out_surr.Dspec;
    end
    Dthr = prctile(Dspec_surr,(1-th)*100,1);
    Dsign_tmp = squeeze(Dspec(i,:,:))>squeeze(Dthr);
    pos_pY0_tmp = repmat(pos_pY0,size(Dsign_tmp,1),1); Dsign_tmp = Dsign_tmp.*pos_pY0_tmp;
    Dsign = (sum(Dsign_tmp,2) ~= 0);

    Isign(i) = (sum(Dsign) ~= 0);
    Rsign(i) = (Dsign(1) ~= 0);
    U1sign(i) = (Dsign(2) ~= 0);
    U2sign(i) = (Dsign(3) ~= 0);
    Ssign(i) = (Dsign(4) ~= 0);
end

%% plot

col = [109 89 122;38 70 83;42 157 143;231 111 81;244 162 97]./255;
thMeas = [Ith;U1th;U2th;Rth;Sth];
Meas = [I;U1;U2;R;S];
signMeas = [Isign;U1sign;U2sign;Rsign;Ssign];
legend_label = {'I(Y;X_1,X_2)','U(Y;X_1)','U(Y;X_2)','R(Y;X_1,X_2)','S(Y;X_1,X_2)'};

figure;
hold on;
for imeas = 1:size(thMeas,1)
    plot(d,Meas(imeas,:),'Color',col(imeas,:),'LineWidth',2,'DisplayName',legend_label{imeas});
    for id = 1:length(d)
        if signMeas(imeas,id) == 1
            scatter(d(id),Meas(imeas,id),'filled','Marker','o','MarkerFaceColor',col(imeas,:),'HandleVisibility','off');
        end
    end
    plot(d,thMeas(imeas,:),':','Color',col(imeas,:),'LineWidth',2,'HandleVisibility','off');
end
legend;
xlabel('d');
ylabel('[nats]');