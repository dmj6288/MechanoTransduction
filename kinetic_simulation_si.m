%Kinetic simulation of the talin rod
%R7R8 modeled as a whole
close all;
clear all;
save_folder=('');
fileprefix='testdata';
domain_size=5;%nm
Rod_contour_length=[69.2 52.4 49.6 52.4 64 60.4 91.6 66.8 63.2 66.4 62.8 ];% in units of nm, calculated by 0.4 nm per residue and number of residues in each talin alpha-helix bundle
Rod_K0_unfold=[4.2E-6 1.7E-8 0.018 4.2E-6 2.5E-5 2.5E-5 4.2E-06 4.2E-6 2.5E-5 2.5E-5 1.7E-8];% in units of  s^-1, zero force unfolding rates from Table 1
Rod_K0_refold=[0.11 0.019 22.2 0.46 1.0 1.0 0.39 0.93 0.93 0.93 0.93 ]; % in units of s^-1, zero force refolding rates from Table 1
Rod_D_unfold=[3.1 3.40 5.7 3.1 4.1 4.1 3.1 3.1 4.1 4.1 3.40];% in units of nm, the unfolding transition distance from Table 1
Rod_D_refold=[18.2 12.5 15.5 4.4 15.7 15.7 13.3 14.5 14.5 14.5 14.5 ]; % in units of nm, the dimension of transition state in refolding from Table 1

n_run=2; %number of repeated simulations to run

%% Define the initial states of talin and 
%1 -> Folded
%2 -> Unfolded
Rod_state_start=[1 1 1 1 1 1 1 1 1 1 1 ];

k_f_bells = @(force_kBT, k_0, x_c) k_0.*exp(force_kBT.*x_c);


%Time=linspace(0,40,600);
%Extension=[linspace(30,550) linspace(550,60) linspace(60,350) linspace(350,150) linspace(150,450) linspace(450,30)];
%% prepare refolding function

file_data_raw=importdata(sprintf('%s%s',fileprefix,'.txt'));
file_data=file_data_raw.data;
Time=file_data(:,1);
Extension=file_data(:,2);

time_extension_spline=spline(Time,Extension);
extension_vs_t=@(tt) ppval(time_extension_spline,tt);
time_speed_spline=fnder(time_extension_spline);
extension_speed_vs_t=@(tt) abs(ppval(time_speed_spline,tt));

ff=[0.01:0.1:10 10:1:100]; %force range to compute force dependent rates
K_2_1_f_spline=cell(size(Rod_K0_refold));
parfor k=1:length(ff)
K_2_1_f_values(k,:)=arrayfun(@talin_domain_kf_f1,ones(size(Rod_K0_refold)).*double(ff(k)),Rod_K0_refold,Rod_D_refold,Rod_contour_length);
end
parfor j=1:length(Rod_K0_refold)
K_2_1_f_spline{j}=spline(ff,K_2_1_f_values(:,j));
end
K_2_1_f_f=@(f) cellfun(@ppval,K_2_1_f_spline,num2cell(ones(size(Rod_K0_refold)).*(f)));
%% Gillespie
min_dis_updateforce=0.1;%maximum step to update force;

%parpool(8);
trajectory=cell(n_run,1);
%h = waitbar(0,'Simulation Progress');
trajectory_detail=cell(n_run,1);
parfor i=1:n_run

Rod_state_current=Rod_state_start;
t=min(Time)+0.5;%start at 0.5s to avoid boundary effects of splined extension traces.
dt_min=0.1; %maximum step size


while t<= max(Time)
    current_force=force_vs_contour_states(extension_vs_t(t), Rod_state_current);
    trajectory{i}=[trajectory{i}; [t double(current_force) sum(Rod_state_current==2)]];
    trajectory_detail{i}=[trajectory_detail{i}; [t double(current_force) Rod_state_current]];

    K_1_2_f=(k_f_bells(current_force./4.1,Rod_K0_unfold,Rod_D_unfold));
    K_2_1_f=abs(K_2_1_f_f(double(current_force)));
    K_sum_current=sum(K_1_2_f(Rod_state_current==1))+sum(K_2_1_f(Rod_state_current==2));
    dt_move=min_dis_updateforce/extension_speed_vs_t(t);
    if dt_move>dt_min
        dt_move=dt_min;
    end
    dt_event=double(exprnd(1/double(K_sum_current)));
    if dt_event>dt_move
        t=t+dt_move;
    else
        
        random_event_no = gendist(double([K_1_2_f.*(Rod_state_current==1) K_2_1_f.*(Rod_state_current==2)]./K_sum_current),1,1);
        Rod_state_current_temp=Rod_state_current;
        if random_event_no <=11 
            Rod_state_current(random_event_no)=2;
        else
            Rod_state_current(random_event_no-11)=1;
        end
        next_force_temp=force_vs_contour_states(extension_vs_t(t+dt_event), Rod_state_current);
        K_1_2_f_next_temp=(k_f_bells(next_force_temp./4.1,Rod_K0_unfold,Rod_D_unfold));
        K_2_1_f_next_temp=abs(K_2_1_f_f(double(next_force_temp)));
        if random_event_no <=11 
            K_next_temp=K_2_1_f_next_temp(random_event_no);
        else
            K_next_temp=K_1_2_f_next_temp(random_event_no-11);
        end
        % Ignore really transient transitions for the sake of simulaiton speed
        if next_force_temp==100 
            Rod_state_current=Rod_state_current_temp;
            t=t+dt_move;
        else
            if K_next_temp>5000 
                Rod_state_current=Rod_state_current_temp;
            end
            t=t+dt_event;
        end
    end
end

end


close all;
for k=1:n_run
    figure();
    set(gcf,'PaperUnits','centimeters')
    set(gcf,'Units','centimeters')
    xSize = 8; ySize = 10;
    xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
    set(gcf,'Position',[xLeft yTop xSize ySize])
    set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
    
    %set(gcf,'OuterPosition',[5 5 xSize ySize])
    subplot('Position',[0.2000    0.7350    0.7000    0.15]);
    %hold on;
    stairs(trajectory{k}(:,1),trajectory{k}(:,2),'k-','linewidth',1);
    set(gca,'LineWidth',1.2,'FontSize',10);
    ylabel('Force (pN)','FontSize',10);
    ylabh = get(gca,'YLabel');
    set(ylabh, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
    xlim([1 max(trajectory{k}(:,1))]);
    title(sprintf('%s%0.2f%s%0.2f%s','Mean Force ( \pm std): ',mean(trajectory{k}(:,2)),' \pm',std(trajectory{k}(:,2)),' pN'),'FontSize',10);
    subplot('Position',[0.2000    0.5650    0.7000    0.15]);
    vbs_mask=[1 2 2 0 0 1 1 1 0 1 1 0];
    trajectory_detail_here=[trajectory_detail{k}(:,3:9) trajectory_detail{k}(:,9) trajectory_detail{k}(:,10:end)]-1;
    num_exposed_vbs=[trajectory_detail{k}(:,1) sum(trajectory_detail_here.*repmat(vbs_mask,length(trajectory_detail_here),1),2)];
    stairs(trajectory{k}(:,1),num_exposed_vbs(:,2),'k-','linewidth',1);
    ylabel('N_{VBS}','FontSize',10);
    ylabh = get(gca,'YLabel');
    set(ylabh, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
    xlim([1 max(trajectory{k}(:,1))]);
    subplot('Position',[0.2000    0.3950    0.7000    0.15]);
    plot(Time,Extension,'k-','LineWidth',1);
    set(gca,'LineWidth',1.2,'FontSize',10);
    ylabel('Ext (nm)','FontSize',10);
    ylabh = get(gca,'YLabel');
    set(ylabh, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
    xlim([1 max(trajectory{k}(:,1))]);
    
    subplot('Position',[0.2000    0.2250    0.7000    0.15]);
    relax_trajectory_sample=[trajectory_detail{k}(:,1) trajectory_detail{k}(:,3:end)-1];
    smooth_local=@(a) smooth(a,50);
    relax_trajectory_heatmap=[];
    trajectory_sample_size=size(relax_trajectory_sample);
    for i=1:trajectory_sample_size(2) %take account of the cooperate unfolding of R7 and R8
        if i<8
            relax_trajectory_heatmap(:,i)=smooth_local(relax_trajectory_sample(:,i));
        elseif i>8
            relax_trajectory_heatmap(:,i+1)=smooth_local(relax_trajectory_sample(:,i));
        elseif i==8 
            relax_trajectory_heatmap(:,i)=smooth_local(relax_trajectory_sample(:,i));
            relax_trajectory_heatmap(:,i+1)=smooth_local(relax_trajectory_sample(:,i));
        end
    end
    colormap('jet');
    imagesc(relax_trajectory_heatmap(:,1),1:(trajectory_sample_size(2)-1),relax_trajectory_heatmap(:,2:end)');
    set(gca,'LineWidth',1.2,'FontSize',10);
    ylabel('Domain No.','FontSize',10);
    ylabh = get(gca,'YLabel');
    set(ylabh, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
    xlim([1 max(relax_trajectory_heatmap(:,1))]);
    xlabel('Time (s)','FontSize',10);
    samexaxis();
    print(sprintf('%s%s%s%d',save_folder,fileprefix,'_kinetic_trajectory_nonbell_11domain_global_boot',k),'-dpdf');
end


