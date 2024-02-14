alpha_list_bd0      = coco_bd_col(bd0, {'alpha'});
beta_list_bd0       = coco_bd_col(bd0, {'beta'});
lambda_list_bd0     = coco_bd_col(bd0, {'lambda_3'});

b2_list_bd0     = coco_bd_col(bd0, {'b2'});
b3_list_bd0     = coco_bd_col(bd0, {'b3'});
T_list_bd0      = coco_bd_col(bd0, {'T'});

par_list = [alpha_list_bd0; beta_list_bd0; lambda_list_bd0; b2_list_bd0; b3_list_bd0; T_list_bd0];
monitor  = [];
for i = 1 : size(par_list,2)
    monitor_tmp =  LCO_multiplier_3D_case(par_list(:,i));
    monitor = [monitor, monitor_tmp];
end
monitor(5,:)
monitor(6,:)
