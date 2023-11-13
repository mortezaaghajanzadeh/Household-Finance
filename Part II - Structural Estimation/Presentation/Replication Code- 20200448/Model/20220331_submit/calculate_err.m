function [err,mmt0] = calculate_err(Par,s0,W0,mmt_actual,s_additional)

    if nargin == 5
        Par.a0 = s_additional;
    end
    
    Par.beta = s0(1);
    Par.W0 = s0(2);
    Par.os = s0(3);
    Par.a1 = s0(4);
    Par.cb0 = s0(5);
    Par.cd0 = s0(6);
    Par.fix = s0(7);

    Par = Grid_exp(Par); 
    shock0 = genshocks(Par); save('Data\shocks_6_500000_1_0_0','shock0');
    solution0 = solve_model_MP(Par,'YES DMP','YES LMP'); 
    [firm0, simu_con0, simu_un0, rec0] = SimulatePanel_short(Par, solution0, shock0);
    
    mmt0= [simu_un0.meanC2V, ...
                 simu_un0.meanN2D, ...
                 simu_un0.stdN2D, ...
                 simu_un0.spreadrd, ...
                 simu_un0.spreadr, ...
                 simu_un0.meanD2A, ...
                 simu_un0.meanFIX1, ...
                 simu_un0.meanLEV, ...
                 simu_un0.meanM2B, ...
                 simu_un0.reg_T_agg, ...
                 simu_un0.reg_B_agg];
    
    for i = 1:length(mmt_actual)
        if mmt_actual(i) == 0
           mmt0(i) = 0;
        end
    end

    err = (mmt0-mmt_actual)*W0*(mmt0'-mmt_actual');
    save('Data\Par_est','Par');
end
    
