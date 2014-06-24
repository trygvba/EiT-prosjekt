F_const = -Fmat_up*uz_up-Fmat_low*uz_low;
FL_const = -2*(Fmat_low*uz_low);
F_cur = @(t) FL_const - Fmat_up*(2*uz_up + (plateDisp(t)+plateDisp(t-dt)))-(plateAcc(t)+plateAcc(t-dt))*F_acc;
Utemp_last = Amod\F_const;
u0 = putDirichletBackV2(Utemp_last,lowerNodes,upperNodes,uz_low,uz_up,x_plates,y_plates);
%U = full(u0);
U = zeros(szU,1);
disp('Starting time integration')
%------------------------------------------------
%       TIME INTEGRATION (NEWMARK)
%------------------------------------------------
n=1;
Uz=zeros(szU/3,steps);
Ux=Uz;
Uy=Uz;
%h = waitbar(0,'Pictures taken');
tic
for i=1:steps
     if i==1||(floor(NumberOfPics*i/steps)>floor(NumberOfPics*(i-1)/steps))
         State_to_vtk(output_folder,vtktitle,n,szU,tetr(:,1:4),p,U);
         n=n+1;
     end
    t =T0+(i-1)*dt;
    Utemp_cur =  K2*(2*Utemp_last + dt*vel + (0.25*dt^2)*(Mmod_inv*(F_cur(t))))-Utemp_last;
    utemp = putDirichletBackV2(Utemp_cur,lowerNodes,upperNodes,uz_low,uz_up+plateDisp(t),x_plates,y_plates);
    U = utemp;
    Ux(:,i) = U(1:3:end-2); % For transmisjonsregning

    Uy(:,i) = U(2:3:end-1); % For transmisjonsregning

    Uz(:,i) = U(3:3:end); % For transmisjonsregning
    vel = vel+ 0.5*dt*((Mmod_inv*(F_cur(t)))-(K1*(Utemp_cur+Utemp_last)));
    Utemp_last = Utemp_cur;
    %waitbar(i/steps);
end
toc
TranNode = max(abs(Uz(MarkerNode,:)))/OLT;% Transmisjonskoeffisient

TranDist = abs(Uz(find((abs(p(:,3))<0.1)),:))/OLT;
TranMax = 1:size(TranDist,1);
for i = 1:length(TranMax)
    TranMax(i) = max(TranDist(i,:));
end
TranAvg = sum(TranMax)/length(TranMax);

[TranNode; TranAvg; 2*sqrt(rho(1))/(sqrt(rho(1))+sqrt(rho(2)))]
plot(Uz(MarkerNode,:))
% tic
% for i=1:1
% %     if i==1||(floor(NumberOfPics*i/steps)>floor(NumberOfPics*(i-1)/steps))
% %         State_to_vtk(output_folder,vtktitle,n,szU,tetr(:,1:4),p,U);
% %         n=n+1;
% %     end
% tic;    t =T0+(i-1)*dt; toc;
%     F_cur = -Fmat_up*(uz_up+plateDisp(t)) - Fmat_low*uz_low - plateAcc(t)*F_acc;  %Fjerd tregest
%     Utemp_cur = K2*(Utemp_last+dt*vel-0.25*dt^2*((K1*Utemp_last)-(Mmod_inv*(F_cur+F_last))));  % Tregest
%     utemp = putDirichletBackV2(Utemp_cur,lowerNodes,upperNodes,uz_low,uz_up+plateDisp(t),x_plates,y_plates); %Tredj tregest
%     U = utemp;
%     %Uz(:,i) = U(3:3:end);
%     vel = vel+ 0.5*dt*((Mmod_inv*(F_cur+F_last))-(K1*(Utemp_cur+Utemp_last)));  % Nest tregest
%     Utemp_last = Utemp_cur;
%     F_last = F_cur;
%     %waitbar(i/steps);
% end
% toc

