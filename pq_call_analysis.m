function pq_call_analysis
% PQ_CALL_ANALYSIS  Used for the analysis presented in 'Telephone training
% to improve ECG quality in remote screening for atrial fibrillation' (2024).
%
%   # Inputs
%   * rec_data_anon.mat - two of these files, one for the SAFER Feas2 dataset, 
%   and one for the Trial dataset.
%
%   # Outputs
%   * Text outputs in command window, and box plots saved to plot folder.
%
%   # Exemplary usage
%
%      pq_call_analysis  % runs the analysis
%
%   # Author
%   Peter H. Charlton, University of Cambridge, August 2024
%
%   # License - MIT
%      Copyright (c) 2024 Peter H. Charlton
%      Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
%      The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
%      THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

% setup
up.paths.rec_data_anon_trial = '/Volumes/T7/SAFER_data/Trial_orig_20231207/rec_data_anon.mat';
up.paths.rec_data_anon_feas2 = '/Volumes/T7/SAFER_data/Feas2_orig_20231102/rec_data_anon.mat';
up.paths.plot_savefolder = '/Users/petercharlton/Library/CloudStorage/GoogleDrive-peterhcharlton@gmail.com/My Drive/Work/Publications/In Preparation/2024 ECG telephone training/202408 PMea Resubmission/figures/';

% do analysis for each phase of SAFER
phases = {'Combined', 'Feas2', 'Trial'};
var_names = {'precall', 'postcall', 'group'};
for phase_no = 1 : length(phases)
    curr_phase = phases{phase_no};

    % load data for current phase
    if strcmp(curr_phase, 'Trial')
        filename = {up.paths.rec_data_anon_trial};
    elseif strcmp(curr_phase, 'Feas2')
        filename = {up.paths.rec_data_anon_feas2};
    elseif strcmp(curr_phase, 'Combined')
        filename = {up.paths.rec_data_anon_trial, ...
            up.paths.rec_data_anon_feas2};
    end
    
    % rename loaded variable
    for file_no = 1 : length(filename)
        eval(['temp_' num2str(file_no) ' = load(filename{file_no}, ''rec_data_anon'');']);
    end
    
    % combine data for combined Feas2 and Trial analysis.
    rec_data_anon = temp_1.rec_data_anon;
    if length(filename)>1
        temp_2.rec_data_anon.ptID = temp_2.rec_data_anon.ptID+max(temp_1.rec_data_anon.ptID);
        temp_2.rec_data_anon.gpID = temp_2.rec_data_anon.gpID+max(temp_1.rec_data_anon.gpID);
        for file_no = 2 : length(filename)
            eval(['rec_data_anon = [rec_data_anon; temp_' num2str(file_no) '.rec_data_anon];']);
        end
    end
    
    % Perform analysis for this phase
    fprintf(['\n~~~~ ', curr_phase ' Analysis ~~~~'])
    use_true_time_of_call = false;
    [precall, postcall, group] = compare_before_and_after_call(rec_data_anon, 0, use_true_time_of_call, curr_phase);

    % make boxplots of % poor quality ECGs per participant
    if strcmp(curr_phase, 'Combined')
        make_boxplots(precall, postcall, group, up);
    end

    % store results
    for var_name = var_names
        eval([curr_phase '.' var_name{1,1} ' = ' var_name{1,1} ';'])
    end

end

do_remainder = 0;

if do_remainder

    % (I think this is now redundant, as the necessary part of the combined analysis is performed above)

    fprintf('\n\n ~~~~ Combined Analysis ~~~~')

    % do analysis for both phases
    for var_name = var_names
        eval([var_name{1,1} ' = [Feas2.' var_name{1,1} '; Trial.' var_name{1,1} '];']);
    end
    test_using_ancova(precall, postcall, group, no_pts);

    % do additional analysis of prop poor quality over time
    make_plot_of_prop_low_quality_over_time(rec_data_anon);

end

fprintf('\n\n --- Finished analysis ---\n')

end

function make_plot_of_prop_low_quality_over_time(a)

pts = unique(a.ptID);

% exclude participants with < 80 ECGs or who received a poor quality call
exc_pts = 0;
if exc_pts
    rows_to_exc = false(height(a),1);
    for pt_no = 1 : length(pts)
        if rem(pt_no,1000) == 0
            fprintf('\n pt no %d', pt_no)
        end
        
        curr_pt = pts(pt_no);
        pt_rows = a.ptID==curr_pt;
        rec_call = max(a.totalNoPQcalls(pt_rows))>0;
        if sum(pt_rows)< 80 || rec_call
            rows_to_exc(pt_rows) = true;
        end
    end
    a = a(~rows_to_exc,:);
    pts = unique(a.ptID);
end

% make bar chart of proportion low quality (days)
days = 0:20;
[prop_lq, no_pts] = deal(nan(length(days),1));
for day_no = 1 : length(days)
    curr_day = days(day_no);
    rel_rows = a.measRelTime>= curr_day & a.measRelTime<= curr_day+1;
    prop_lq(day_no) = 100*mean(a.tag_orig_Poor_Quality(rel_rows));
    no_pts(day_no) = length(unique(a.ptID(rel_rows)));
end
bar(days,prop_lq)
ftsize = 14;
ylabel('Proportion of ECGs which were low-quality (%)', 'FontSize', ftsize)
xlabel('Day', 'FontSize', ftsize)
set(gca, 'FontSize', ftsize)
box off
grid on

% make bar chart of proportion low quality (ecgs)
figure
ecgs = 0:80;
[prop_lq, no_pts] = deal(nan(length(ecgs),1));
for ecg_no = 1 : length(ecgs)
    curr_ecg = ecgs(ecg_no);
    rel_rows = a.measNo == curr_ecg;
    prop_lq(ecg_no) = 100*mean(a.tag_orig_Poor_Quality(rel_rows));
    no_pts(ecg_no) = length(unique(a.ptID(rel_rows)));
end
bar(ecgs,prop_lq)
ftsize = 14;
ylabel('Proportion of ECGs which were low-quality (%)', 'FontSize', ftsize)
xlabel('ECG no.', 'FontSize', ftsize)
set(gca, 'FontSize', ftsize)
box off
grid on


end

function [precall, postcall, group] = compare_before_and_after_call(a, verbose_log, use_true_time_of_call, curr_phase)

pts = unique(a.ptID);

fprintf('\n\n Outputting initial statistics for %s dataset (for Table 1)', curr_phase)

% - prepare stats
no_pts_at_least1_pq = 0;
no_pts_at_least25pc_pq = 0;
no_female = 0;
for pt_no = 1 : length(pts)
    curr_rows = a.ptID == pts(pt_no);
    no_pq_pt = sum(a.tag_orig_Poor_Quality(curr_rows));
    if no_pq_pt>0
        no_pts_at_least1_pq = no_pts_at_least1_pq + 1;
    end
    if no_pq_pt>0.25*sum(curr_rows)
        no_pts_at_least25pc_pq = no_pts_at_least25pc_pq + 1;
    end
    if sum(a.gender(curr_rows) == 1)
        no_female = no_female+1;
    end
end
ages = a.age(~isnan(a.age));

% - output stats
fprintf('\n - No. participants: %d', length(pts))
fprintf('\n - Age, mean (SD): %.1f (%.1f)', mean(ages), std(ages))
fprintf('\n - Gender, no (%% female): %d (%.1f%%)', no_female, 100*no_female/length(pts))
fprintf('\n - No. Practices: %d', length(unique(a.gpID)))
fprintf('\n - Total no. ECGs: %d', length(a.measNo))
no_pq = length(a.measNo(a.tag_orig_Poor_Quality));
fprintf('\n - Total no. poor-quality ECGs (%%): %d (%.1f%%)', no_pq, 100*no_pq/length(a.measNo))
fprintf('\n - No. participants with at at least one poor-quality ECG: %d (%.1f%%)', no_pts_at_least1_pq, 100*no_pts_at_least1_pq/length(pts))
fprintf('\n - No. participants with at at least 25%% poor-quality ECGs: %d (%.1f%%)', no_pts_at_least25pc_pq, 100*no_pts_at_least25pc_pq/length(pts))

fprintf('\n\n Comparing proportions of low quality ECGs before and after call')

thresh_breach_t = nan(length(pts),1);
rec_call = false(length(pts),1);
[prop_lq_before_pt, prop_lq_after_pt] = deal(nan(length(pts),1));

do_new_definition = true; % implemented on 23-July-2024 in response to A Dymond's email.

% cycle through each participant
outputted_ref_times = false;
for pt_no = 1 : length(pts)
    curr_pt = pts(pt_no);

    % find ECGs for this participant
    rel_rows = a.ptID == curr_pt;

    % skip this participant if there are less than three poor quality ECGs for this participant (because they won't breach the threshold)
    if sum(a.tag_orig_Poor_Quality(rel_rows))<3
        continue
    end

    rel_rows = find(rel_rows);

    if do_new_definition
%         % A dymond definition
%         day_4_start_RelTime = 3 + (1-a.measTime(rel_rows(1)));
%         day_11_start_RelTime = 10 + (1-a.measTime(rel_rows(1)));
        % J Brimicombe definition
        day_4_start_RelTime = 2 + (1-a.measTime(rel_rows(1)));
        day_11_start_RelTime = 9 + (1-a.measTime(rel_rows(1)));
    else
        day_4_start_RelTime = 3;
        day_11_start_RelTime = 10;
    end
    if ~outputted_ref_times
        fprintf('\n - Start of day 4 to end of day 10 defined as (e.g.) time %.1f to %.1f', day_4_start_RelTime, day_11_start_RelTime)
        outputted_ref_times = true;
    end

    % setup to find out the proportion of poor quality ECGs before and after
    no_lq = 0;

    % go through each ECG for this participant, looking for the time at which they received a call or breached the threshold
    for row_no = 1 : length(rel_rows)
        % skip if this ECG was measured after the end of the threshold period
        if a.measRelTime(rel_rows(row_no)) > day_11_start_RelTime
            break
        end
        % count the number of poor quality ECGs so far
        no_lq = no_lq + a.tag_orig_Poor_Quality(rel_rows(row_no));
        no_lq_temp = sum(a.tag_orig_Poor_Quality(rel_rows(1:row_no)));
        if no_lq_temp~=no_lq
            error('these should match')
        end
        % calculate the proportion poor quality so far
        prop_lq_before = 100*no_lq/row_no;
        t = a.measRelTime(rel_rows(row_no));
        % see whether the participant received a call or if the threshold was breached within the relevant period
        rec_call(pt_no) = max(a.totalNoPQcalls(rel_rows))>0;
        if rec_call(pt_no) || (prop_lq_before>=25 && t < day_11_start_RelTime && t>=day_4_start_RelTime)
            if verbose_log
                fprintf('\n Participant %d ', curr_pt);
            end
            %  adjust the time if this participant received a call
            if rec_call(pt_no)
                if verbose_log
                    fprintf('RECEIVED call, and ')
                end
                % adjust row no to be the ECG immediately before the call was received
                row_no = find(a.noPQcalls(rel_rows),1)-1;
                t = a.measRelTime(rel_rows(row_no));
                % recalculate prop poor quality for this new time
                no_lq = sum(a.tag_orig_Poor_Quality(rel_rows(1:row_no)));
                prop_lq_before = 100*no_lq/row_no;
            end
            % note down the time at which the threshold was breached
            thresh_breach_t(pt_no) = t;
            if verbose_log
                fprintf('breached at %.2f days with prop %.1f%%.', thresh_breach_t(pt_no), prop_lq_before)
            end
            
            prop_lq_after_pt(pt_no) = 100*sum(a.tag_orig_Poor_Quality(rel_rows(row_no+1:end)))/(length(rel_rows)-row_no);
            prop_lq_before_pt(pt_no) = prop_lq_before;
            if ~use_true_time_of_call
                curr_rel_rows = rel_rows(a.measRelTime(rel_rows)<day_4_start_RelTime);
                prop_lq_before_pt(pt_no) = 100*sum(a.tag_orig_Poor_Quality(curr_rel_rows))/length(curr_rel_rows);
                curr_rel_rows = rel_rows(a.measRelTime(rel_rows)>=day_11_start_RelTime);
                prop_lq_after_pt(pt_no) = 100*sum(a.tag_orig_Poor_Quality(curr_rel_rows))/length(curr_rel_rows);
            end
            
            if verbose_log
                fprintf(' Prop before, after: %.1f%%, %.1f%%', prop_lq_before_pt(pt_no), prop_lq_after_pt(pt_no))
            end
            break
        end
    end
end

% test using signed rank test
test_using_signed_rank_test(rec_call, prop_lq_before_pt, prop_lq_after_pt, thresh_breach_t);

% test using ANOCVA
rec_call_els = rec_call;
didnt_rec_call_els = ~rec_call & ~isnan(thresh_breach_t) & ~isnan(prop_lq_after_pt);
precall = [prop_lq_before_pt(rec_call_els); prop_lq_before_pt(didnt_rec_call_els)];
postcall = [prop_lq_after_pt(rec_call_els); prop_lq_after_pt(didnt_rec_call_els)];
group = [repmat({'received call'}, sum(rec_call_els),1); repmat({'didnt receive call'}, sum(didnt_rec_call_els), 1)];
group = categorical(group);
no_pts = length(rec_call_els);
test_using_ancova(precall, postcall, group, no_pts);

end

function test_using_signed_rank_test(rec_call, prop_lq_before_pt, prop_lq_after_pt, thresh_breach_t)

% output results
fprintf('\n\n --- Results (using Wilcoxon signed rank test, in Table 2) ---')
sig_level = 0.05;
rel_els = rec_call;
fprintf('\n No. received call and breached thresh:\n    %d (%.1f%%) with median (IQR) of %.1f (%.1f - %.1f) %% of ECGs poor quality before, and %.1f (%.1f - %.1f) %% after', sum(rel_els), 100*sum(rel_els)/length(rel_els), median(prop_lq_before_pt(rel_els)), quantile(prop_lq_before_pt(rel_els),0.25), quantile(prop_lq_before_pt(rel_els),0.75), median(prop_lq_after_pt(rel_els)), quantile(prop_lq_after_pt(rel_els),0.25), quantile(prop_lq_after_pt(rel_els),0.75))
[p,h,stats] = signrank(prop_lq_before_pt(rel_els), prop_lq_after_pt(rel_els));
if p<sig_level
    fprintf('\n    - SIGNIFICANT (p = %.3f)', p)
else
    fprintf('\n    - non-significant (p = %.3f)', p)
end
rel_els = ~rec_call & ~isnan(thresh_breach_t) & ~isnan(prop_lq_after_pt);
fprintf('\n No. did not receive call and breached thresh:\n    %d (%.1f%%) with %.1f (%.1f - %.1f) %% of ECGs poor quality before, and %.1f (%.1f - %.1f) %% after', sum(rel_els), 100*sum(rel_els)/length(rel_els), median(prop_lq_before_pt(rel_els)), quantile(prop_lq_before_pt(rel_els),0.25), quantile(prop_lq_before_pt(rel_els),0.75), median(prop_lq_after_pt(rel_els)), quantile(prop_lq_after_pt(rel_els),0.25), quantile(prop_lq_after_pt(rel_els),0.75))
[p,h,stats] = signrank(prop_lq_before_pt(rel_els), prop_lq_after_pt(rel_els));
if p<sig_level
    fprintf('\n    - SIGNIFICANT (p = %.3f)', p)
else
    fprintf('\n    - non-significant (p = %.3f)', p)
end

end

function test_using_ancova(precall, postcall, group, no_pts)

fprintf('\n\n --- Results (using ANCOVA, in Table 3) ---\n')

sig_level = 0.05;

% report no. subjects
fprintf('No. subjects: %d (%.1f%%)', length(precall), 100*length(precall)/no_pts);
cats = categories(group);
cat_counts = countcats(group);
for cat_no = 1 : length(cats)
    fprintf('\n   - %d (%.1f%%) %s', cat_counts(cat_no), 100*cat_counts(cat_no)/no_pts, cats{cat_no})
end

% report means pre and post call
for cat_no = 1 : length(cats)
    curr_cat = cats{cat_no};
    fprintf('\n For those who %s:', curr_cat)
    rel_els = group == curr_cat;
    rel_vals = precall(rel_els);
    mean_val = mean(rel_vals);
    se = std(rel_vals)/sqrt(length(rel_vals));
    upper_ci = mean_val + (1.96*se);
    lower_ci = mean_val - (1.96*se);
    fprintf('\n   - Mean prop poor quality pre-call: %.1f (%.1f-%.1f) %%', mean_val, lower_ci, upper_ci)
    rel_vals = postcall(rel_els);
    mean_val = mean(rel_vals);
    se = std(rel_vals)/sqrt(length(rel_vals));
    upper_ci = mean_val + (1.96*se);
    lower_ci = mean_val - (1.96*se);
    fprintf('\n   - Mean prop poor quality post-call: %.1f (%.1f-%.1f) %%', mean_val, lower_ci, upper_ci)
    rel_vals = precall(rel_els)-postcall(rel_els);
    mean_val = mean(rel_vals);
    se = std(rel_vals)/sqrt(length(rel_vals));
    upper_ci = mean_val + (1.96*se);
    lower_ci = mean_val - (1.96*se);
    fprintf('\n   - Mean difference: %.1f (%.1f-%.1f) %%', mean_val, lower_ci, upper_ci)
    fprintf('\n   - Prop reduction: %.1f%%', 100*mean(precall(rel_els)-postcall(rel_els))/mean(precall(rel_els)))
end

% do stats test
[h,atab,ctab,stats] = aoctool(precall,postcall,group, sig_level, '', '', '', 'on', 'parallel lines');

% report results of stats test
fprintf('\n Equation is:   post-call prop = %.0f + %.2f x pre-call prop + %.1f x group', min(stats.intercepts), ctab{end,2}, 2*ctab{3,2});
se = ctab{3,3}*2;
mean_val = 2*ctab{3,2};
half_96_ci = 1.96* se;
fprintf('\n Difference between means is %.1f, with a 95%% CI of %.1f - %.1f', mean_val, mean_val - half_96_ci, mean_val + half_96_ci);
z = mean_val/se;
p = exp(-0.717*z - 0.416*z^2);
fprintf('\n P-value is: %.3f\n', p)

end

function make_boxplots(precall, postcall, group, up)

ftsize = 20;

groups = {'Received training', 'Did not receive training'};

for group_no = 1 : length(groups)

    curr_group = groups{group_no};

    % - extract relevant data
    if strcmp(curr_group, 'Received training')
        rel_els = group == 'received call';
    else
        rel_els = group == 'didnt receive call';
    end
    precall_data = precall(rel_els);
    postcall_data = postcall(rel_els);
    all_data = [precall_data; postcall_data];
    all_groups = [repmat(categorical({'Days 1-3'}), [length(precall_data),1]); repmat(categorical({'Days 11-21'}), [length(postcall_data),1])];

    % make box plot
    figure('Position', [20, 20, 600,600])
    h = boxchart(all_groups, all_data);
    h.BoxWidth = 0.75;

    % tidy up box plot
    set(gca, 'FontSize', ftsize, 'YGrid', 'on');
    xlabel('Time period during screening')
    ylabel('% poor quality ECGs per participant')
    box off

    % add title
    title(curr_group, 'FontSize', ftsize)

    % save box plot
    filepath = [up.paths.plot_savefolder, 'boxplot_', strrep(curr_group, ' ', '_')];
    print(filepath, '-dpng')

end


end






