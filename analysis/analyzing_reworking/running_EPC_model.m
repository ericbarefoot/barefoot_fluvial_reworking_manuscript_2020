

ii = 4;    
avul = [0 1 2];
IR = 0.4:0.1:0.9;
vertagg_rate = 0.005:0.02:0.765;
% do another run with vertagg rate from 0.005 to 0.045 by
% 0.01
fp_lat = 15;
walk_mag = 20;
latmob = 4:1:12;
rel_window_mult = 10;

runmat = combvec(avul, IR, vertagg_rate, fp_lat, walk_mag, latmob, rel_window_mult);

[~, jj] = size(runmat);

k = 1;
jk = 1;

success_params = nan(size(runmat));
failed_params = nan(size(runmat));

qs = [1 floor(jj/4) floor(jj/2) floor(jj*3/4) jj];

f = waitbar(0,'running many models.');
for jjj = 6053:qs(ii+1)
    in_str.nits = 1000;
    in_str.mw = 500;
    in_str.avul = runmat(1,jjj); % 1 = Compensational, 2 = Random, 0 = clustered
    in_str.ch_z = 1;
    in_str.ch_w_i = 11;
    in_str.IR = runmat(2,jjj);
    in_str.randint = 20;
    in_str.vertagg_rate = runmat(3,jjj);
    in_str.fp_lat = runmat(4,jjj);
    in_str.walk_mag = runmat(5,jjj);
    in_str.wwide = 0.5;
    in_str.wthick = 0.5;
    in_str.plot_bool = false;
    in_str.latmob = runmat(6,jjj);
    in_str.fp_filling = 'local';
    in_str.rel_window_mult = runmat(7,jjj);
    in_str.store_topo = false;
    
    try
        out_str(k) = StickModel_BarPres(in_str);
        success_params(:,k) = runmat(:,jjj);
        k = k+1;
    catch
        warning('The model failed for a run. Skipping and moving on')
        failed_params(:,jk) = runmat(:,jjj);
        jk = jk+1;
    end
    waitbar((jjj-qs(ii))/(qs(ii+1)-qs(ii)))
end
s_drop = all(isnan(success_params));
f_drop = all(isnan(failed_params));
success_params(:,s_drop) = [];
failed_params(:,f_drop) = [];
close(f)
save(sprintf('test_results_20191202_%db.mat', ii), 'out_str', 'runmat', 'success_params', 'failed_params', '-v7.3')



% cla
% plot_strat_str(out_str, true, 100)
