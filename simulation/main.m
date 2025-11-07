clc;
clear;
clf;
par=parpool('local',parcluster('local').NumWorkers, 'IdleTimeout', 60*24);



all=tic;
[L1_fenbu_cankaol,L2_fenbu_cankaol,L3_fenbu_cankaol,L4_fenbu_cankaol,P4_shang_cankaol,P4_xia_cankaol,xiao_cankaol,xiao_tuiyanl,Px_cankaol,Px_tuiyanl]...
         =zhuchengxu(2,21000,800e3);
tend=toc(all);
tend=tend/3600;
disp(tend)


delete(par)
% for i=1:imax
%     [L1_fenbu_cankaol,L2_fenbu_cankaol,L3_fenbu_cankaol,L4_fenbu_cankaol,P4_shang_cankaol,P4_xia_cankaol,xiao_cankaol,xiao_tuiyanl,Px_cankaol,Px_tuiyanl]...
%         =zhuchengxu(i,4,700e3);
%     if(i==1)
%         L1_fenbu_cankao=zeros([size(L1_fenbu_cankaol),imax]);
%         L2_fenbu_cankao=zeros([size(L2_fenbu_cankaol),imax]);
%         L3_fenbu_cankao=zeros([size(L3_fenbu_cankaol),imax]);
%         L4_fenbu_cankao=zeros([size(L4_fenbu_cankaol),imax]);
%         P4_shang_cankao=zeros([size(P4_shang_cankaol),imax]);
%         P4_xia_cankao=zeros([size(P4_xia_cankaol),imax]);
%         xiao_cankao=zeros([size(xiao_cankaol),imax]);
%         xiao_tuiyan=zeros([size(xiao_tuiyanl),imax]);
%         Px_cankao=zeros([size(Px_cankaol),imax]);
%         Px_tuiyan=zeros([size(Px_tuiyanl),imax]);
%     end
%     L1_fenbu_cankao(:,:,:,i)=L1_fenbu_cankaol;
%     L2_fenbu_cankao(:,:,:,i)=L2_fenbu_cankaol;
%     L3_fenbu_cankao(:,:,:,i)=L3_fenbu_cankaol;
%     L4_fenbu_cankao(:,:,:,i)=L4_fenbu_cankaol;
%     P4_shang_cankao(:,:,i)=P4_shang_cankaol;
%     P4_xia_cankao(:,:,i)=P4_xia_cankaol;
%     xiao_cankao(:,:,:,i)=xiao_cankaol;
%     xiao_tuiyan(:,:,:,i)=xiao_tuiyanl;
%     Px_cankao(:,:,i)=Px_cankaol;
%     Px_tuiyan(:,:,i)=Px_tuiyanl;
% end
% L1_fenbu_cankao=mean(L1_fenbu_cankao,4);
% L2_fenbu_cankao=mean(L2_fenbu_cankao,4);
% L3_fenbu_cankao=mean(L3_fenbu_cankao,4);
% L4_fenbu_cankao=mean(L4_fenbu_cankao,4);
% P4_shang_cankao=mean(P4_shang_cankao,3);
% P4_xia_cankao=mean(P4_xia_cankao,3);
% xiao_cankao=mean(xiao_cankao,4);
% xiao_tuiyan=mean(xiao_tuiyan,4);
% Px_cankao=mean(Px_cankao,3);
% Px_tuiyan=mean(Px_tuiyan,3);
% %绘制平均值的图像
% %参考值
% flag_cankao=0;
% pic_num=1;
% ns=3000;
% h0=600e3;
% dR_fenbu=1e3;
% huitu_jiange=max(1,floor(ns/400));
% lo_fenbu=6378e3+200e3;
% hi_fenbu=6378e3+1000e3;
% Ts=2*pi*sqrt((6378e3+h0)^3/398600e9);
% while 1
%     if(pic_num>ns+1)
%         break;
%     end
%     f_tu=figure(1);
%     if(pic_num==1)
%         clf;
%     end
%     plot(L1_fenbu_cankao(1,:,pic_num),L1_fenbu_cankao(2,:,pic_num),L2_fenbu_cankao(1,:,pic_num),...
%         L2_fenbu_cankao(2,:,pic_num),L3_fenbu_cankao(1,:,pic_num),L3_fenbu_cankao(2,:,pic_num),L4_fenbu_cankao(1,:,pic_num),L4_fenbu_cankao(2,:,pic_num));
%     ylabel('轨道半径（km）');
%     xlabel('位于对应轨道半径区间的概率');
%     legend('碎片特征长度大于10cm的参考值',...
%         '碎片特征长度在3cm到10cm之间的参考值',...
%         '碎片特征长度在1.8cm到3cm之间的参考值',...
%         '碎片特征长度在1cm到1.8cm之间的参考值');
%     title({['演化时间=',num2str((pic_num-1)*Ts/86400),'天'];...
%         ['碎片特征长度大于10cm的参考值：',num2str(hi_fenbu/1e3),'km之上的概率=',num2str(P4_shang_cankao(1,pic_num)),'；',num2str(lo_fenbu/1e3),'km之下的概率=',num2str(P4_xia_cankao(1,pic_num))];...
%         ['碎片特征长度在3cm到10cm之间的参考值：',num2str(hi_fenbu/1e3),'km之上的概率=',num2str(P4_shang_cankao(2,pic_num)),'；',num2str(lo_fenbu/1e3),'km之下的概率=',num2str(P4_xia_cankao(2,pic_num))];...
%         ['碎片特征长度在1.8cm到3cm之间的参考值：',num2str(hi_fenbu/1e3),'km之上的概率=',num2str(P4_shang_cankao(3,pic_num)),'；',num2str(lo_fenbu/1e3),'km之下的概率=',num2str(P4_xia_cankao(3,pic_num))];...
%         ['碎片特征长度在1cm到1.8cm之间的参考值：',num2str(hi_fenbu/1e3),'km之上的概率=',num2str(P4_shang_cankao(4,pic_num)),'；',num2str(lo_fenbu/1e3),'km之下的概率=',num2str(P4_xia_cankao(4,pic_num))]});
%     set(gca,'FontSize',16);
%     if(flag_cankao==0)
%         axis([0 0.03 lo_fenbu/1e3 hi_fenbu/1e3]);
%     else
%         axis([0 0.01 lo_fenbu/1e3 hi_fenbu/1e3]);
%     end
%     set(f_tu, 'unit', 'normalized', 'position', [0,0,1,1]);
%     drawnow;
%     F_tu=getframe(gcf);
%     I_tu=frame2im(F_tu);
%     [I_tu,map]=rgb2ind(I_tu,256);
%     if(pic_num == 1)
%         imwrite(I_tu,map,['pingjun',',cankao_duo.gif'],'gif', 'Loopcount',inf,'DelayTime',0.05);
%         clf;
%     else
%         if(pic_num+huitu_jiange>ns+1)
%             if(flag_cankao==0)
%                 imwrite(I_tu,map,['pingjun',',cankao_duo.gif'],'gif','WriteMode','append','DelayTime',2);
%                 clf;
%             else
%                 imwrite(I_tu,map,['pingjun',',cankao_duo.gif'],'gif','WriteMode','append','DelayTime',5);
%             end
%         else
%             imwrite(I_tu,map,['pingjun',',cankao_duo.gif'],'gif','WriteMode','append','DelayTime',0.05);
%             clf;
%         end
%     end
%     if(flag_cankao==0&&pic_num+huitu_jiange>ns+1)
%         flag_cankao=1;
%         pic_num=pic_num-huitu_jiange;%特殊情况，不是错误，是因为需要重复绘制最后一幅图
%     end
%     pic_num=pic_num+huitu_jiange;
% end
% %参考值和推演值的对比
% KL_sanduPQ=zeros(1,ns+1);
% KL_sanduQP=zeros(1,ns+1);
% JS_sandu=zeros(1,ns+1);
% for i=1:ns+1
%     P_KL=xiao_cankao(1,:,i);
%     Q_KL=xiao_tuiyan(1,:,i);
%     for j=1:(hi_fenbu-lo_fenbu)/dR_fenbu
%         if(P_KL(j)~=0)
%             KL_sanduPQ(i)=KL_sanduPQ(i)+P_KL(j)*log(P_KL(j)/Q_KL(j));
%             JS_sandu(i)=JS_sandu(i)+0.5*P_KL(j)*log(P_KL(j)*2/(P_KL(j)+Q_KL(j)));
%         end
%         if(Q_KL(j)~=0)
%             KL_sanduQP(i)=KL_sanduQP(i)+Q_KL(j)*log(Q_KL(j)/P_KL(j));
%             JS_sandu(i)=JS_sandu(i)+0.5*Q_KL(j)*log(Q_KL(j)*2/(P_KL(j)+Q_KL(j)));
%         end
%     end
%     if(Px_cankao(i,1)~=0)
%         KL_sanduPQ(i)=KL_sanduPQ(i)+Px_cankao(i,1)*log(Px_cankao(i,1)/Px_tuiyan(i,1));
%         JS_sandu(i)=JS_sandu(i)+0.5*Px_cankao(i,1)*log(Px_cankao(i,1)*2/(Px_cankao(i,1)+Px_tuiyan(i,1)));
%     end
%     if(Px_cankao(i,2)~=0)
%         KL_sanduPQ(i)=KL_sanduPQ(i)+Px_cankao(i,2)*log(Px_cankao(i,2)/Px_tuiyan(i,2));
%         JS_sandu(i)=JS_sandu(i)+0.5*Px_cankao(i,2)*log(Px_cankao(i,2)*2/(Px_cankao(i,2)+Px_tuiyan(i,2)));
%     end
%     if(Px_tuiyan(i,1)~=0)
%         KL_sanduQP(i)=KL_sanduQP(i)+Px_tuiyan(i,1)*log(Px_tuiyan(i,1)/Px_cankao(i,1));
%         JS_sandu(i)=JS_sandu(i)+0.5*Px_tuiyan(i,1)*log(Px_tuiyan(i,1)*2/(Px_cankao(i,1)+Px_tuiyan(i,1)));
%     end
%     if(Px_tuiyan(i,2)~=0)
%         KL_sanduQP(i)=KL_sanduQP(i)+Px_tuiyan(i,2)*log(Px_tuiyan(i,2)/Px_cankao(i,2));
%         JS_sandu(i)=JS_sandu(i)+0.5*Px_cankao(i,2)*log(Px_cankao(i,2)*2/(Px_cankao(i,2)+Px_tuiyan(i,2)))+0.5*Px_tuiyan(i,2)*log(Px_tuiyan(i,2)*2/(Px_cankao(i,2)+Px_tuiyan(i,2)));
%     end
% end
% flag_tuiyan=0;
% pic_num=1;
% while 1
%     if(pic_num>ns+1)
%         break;
%     end
%     f_tu=figure(2);
%     if(pic_num==1)
%         clf;
%     end
%     plot(xiao_cankao(1,:,pic_num),xiao_cankao(2,:,pic_num),'blue',xiao_tuiyan(1,:,pic_num),xiao_tuiyan(2,:,pic_num),'magenta');
%     ylabel('轨道半径（km）');
%     xlabel('位于对应轨道半径区间的概率');
%     legend('小碎片（特征长度小于10cm）的参考值','小碎片（特征长度小于10cm）的推演值');
%     title({['演化时间=',num2str((pic_num-1)*Ts/86400),'天'];['JS散度=',num2str(JS_sandu(pic_num))];...
%         ['参考值：',num2str(hi_fenbu/1e3),'km之上的概率=',num2str(Px_cankao(pic_num,1)),'；',num2str(lo_fenbu/1e3),'km之下的概率=',num2str(Px_cankao(pic_num,2))];...
%         ['推演值：',num2str(hi_fenbu/1e3),'km之上的概率=',num2str(Px_tuiyan(pic_num,1)),'；',num2str(lo_fenbu/1e3),'km之下的概率=',num2str(Px_tuiyan(pic_num,2))]});
%     set(gca,'FontSize',16);
%     if(flag_tuiyan==0)
%         axis([0 0.03 lo_fenbu/1e3 hi_fenbu/1e3]);
%     else
%         axis([0 0.01 lo_fenbu/1e3 hi_fenbu/1e3]);
%     end
%     set(f_tu, 'unit', 'normalized', 'position', [0,0,1,1]);
%     drawnow;
%     F_tu=getframe(gcf);
%     I_tu=frame2im(F_tu);
%     [I_tu,map]=rgb2ind(I_tu,256);
%     if(pic_num == 1)
%         imwrite(I_tu,map,['pingjun','.gif'],'gif', 'Loopcount',inf,'DelayTime',0.05);
%         clf;
%     else
%         if(pic_num+huitu_jiange>ns+1)
%             if(flag_tuiyan==0)
%                 imwrite(I_tu,map,['pingjun','.gif'],'gif','WriteMode','append','DelayTime',2);
%                 clf;
%             else
%                 imwrite(I_tu,map,['pingjun','.gif'],'gif','WriteMode','append','DelayTime',5);
%             end
%         else
%             imwrite(I_tu,map,['pingjun','.gif'],'WriteMode','append','DelayTime',0.05);
%             clf;
%         end
%     end
%     if(flag_tuiyan==0&&pic_num+huitu_jiange>ns+1)
%         flag_tuiyan=1;
%         pic_num=pic_num-huitu_jiange;%特殊情况，不是错误，是因为需要重复绘制最后一幅图
%     end
%     pic_num=pic_num+huitu_jiange;
% end
% save('pingjun');
% delete(par);