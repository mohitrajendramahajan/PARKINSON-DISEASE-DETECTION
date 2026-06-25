/*=============================================================
    Parkinson's Disease Prediction Project
    Advanced Version (WOE + K-Fold + ROC + FULL OUTPUTS)
=============================================================*/

ods pdf file="/home/u64477338/sasuser.v94/sas project ss/parkinsons_full_output.pdf";
ods graphics on;

/*=============================================================
    [0] Import Data
=============================================================*/
proc import datafile="/home/u64477338/sasuser.v94/sas project ss/Parkinsson disease.csv"
    out=parkinsons_raw
    dbms=csv
    replace;
    guessingrows=max;
run;

proc print data=parkinsons_raw(obs=10);
    title "Raw Data Sample";
run;


/*=============================================================
    [1] Rename Variables
=============================================================*/
data parkinsons_renamed;
    set parkinsons_raw;

    rename
        'MDVP:Fo(Hz)'n = avg_fre
        'MDVP:Fhi(Hz)'n = max_fre
        'MDVP:Flo(Hz)'n = min_fre
        'MDVP:Jitter(%)'n = var_fre1
        'MDVP:Shimmer'n = var_amp1;
run;

proc contents data=parkinsons_renamed;
    title "Dataset After Renaming";
run;


/*=============================================================
    [2] Data Cleaning
=============================================================*/
data parkinsons_clean;
    set parkinsons_renamed;

    if missing(status) then delete;
run;

proc freq data=parkinsons_clean;
    title "Status Distribution";
    tables status;
run;


/*=============================================================
    [3] Binning for WOE
=============================================================*/
data parkinsons_binned;
    set parkinsons_clean;

    length avgfre_bin $20;

    if avg_fre < 120 then avgfre_bin = "Low";
    else if avg_fre < 180 then avgfre_bin = "Medium";
    else avgfre_bin = "High";
run;

proc freq data=parkinsons_binned;
    title "Binning vs Status";
    tables avgfre_bin*status;
run;


/*=============================================================
    [4] WOE Calculation
=============================================================*/
proc sql;
    create table woe_table as
    select avgfre_bin,
        sum(case when status=0 then 1 else 0 end) as Good,
        sum(case when status=1 then 1 else 0 end) as Bad
    from parkinsons_binned
    group by avgfre_bin;
quit;

proc print data=woe_table;
    title "WOE Raw Counts";
run;

proc sql;
    create table totals as
    select sum(Good) as total_good,
           sum(Bad) as total_bad
    from woe_table;
quit;

data woe_final;
    if _n_=1 then set totals;
    set woe_table;

    good_dist = Good/total_good;
    bad_dist = Bad/total_bad;

    if good_dist>0 and bad_dist>0 then WOE = log(good_dist/bad_dist);
run;

proc print data=woe_final;
    title "WOE Final Table";
run;


/*=============================================================
    [5] Join WOE
=============================================================*/
proc sql;
    create table parkinsons_woe as
    select a.*, b.WOE
    from parkinsons_binned a
    left join woe_final b
    on a.avgfre_bin=b.avgfre_bin;
quit;

proc print data=parkinsons_woe(obs=10);
    title "WOE Joined Dataset";
run;


/*=============================================================
    [6] Normalization
=============================================================*/
proc means data=parkinsons_woe noprint;
    var avg_fre max_fre;
    output out=norm min=min1 min2 max=max1 max2;
run;

data parkinsons_final;
    if _n_=1 then set norm;
    set parkinsons_woe;

    norm_avg = (avg_fre-min1)/(max1-min1);
    norm_max = (max_fre-min2)/(max2-min2);
run;

proc print data=parkinsons_final(obs=10);
    title "Final Dataset Sample";
run;


/*=============================================================
    [7] Train/Test Split
=============================================================*/
proc surveyselect data=parkinsons_final
    out=train_cv
    samprate=0.8
    seed=123
    outall;
run;

data train_data test_data;
    set train_cv;

    if Selected=1 then output train_data;
    else output test_data;
run;

proc freq data=train_data;
    title "Training Data Distribution";
    tables status;
run;

proc freq data=test_data;
    title "Test Data Distribution";
    tables status;
run;


/*=============================================================
    [8] K-Fold Setup
=============================================================*/
data train_cv;
    set train_data;
    Fold = mod(_n_-1,5)+1;
run;

proc freq data=train_cv;
    title "Fold Distribution";
    tables Fold;
run;


/*=============================================================
    [9] K-Fold Cross Validation
=============================================================*/
%macro kfold;

%do k=1 %to 5;

    data cv_train cv_valid;
        set train_cv;

        if Fold=&k then output cv_valid;
        else output cv_train;
    run;

    proc logistic data=cv_train descending;
        title "Fold &k Training Model";
        model status = norm_avg norm_max WOE;
    run;

    proc logistic data=cv_train;
        score data=cv_valid out=cv_score_&k;
    run;

    proc print data=cv_score_&k(obs=5);
        title "Fold &k Validation Predictions";
    run;

%end;

%mend;

%kfold;


/*=============================================================
    [10] Final Model
=============================================================*/
ods output ParameterEstimates=final_coef;
proc logistic data=train_data descending outmodel=model_final;
    title "Final Logistic Model";
    model status = norm_avg norm_max WOE;
run;

proc print data=final_coef;
    title "Final Model Coefficients";
run;


/*=============================================================
    [11] Test Scoring
=============================================================*/
proc logistic inmodel=model_final;
    score data=test_data out=test_scored;
run;

proc print data=test_scored(obs=10);
    title "Test Predictions";
run;


/*=============================================================
    [12] Evaluation
=============================================================*/
data eval;
    set test_scored;

    if P_1>=0.5 then pred=1;
    else pred=0;
run;

proc freq data=eval;
    title "Confusion Matrix";
    tables status*pred / norow nocol nopercent;
run;

/* Metrics */
proc sql;
    create table metrics as
    select
        sum(case when status=1 and pred=1 then 1 else 0 end) as TP,
        sum(case when status=0 and pred=0 then 1 else 0 end) as TN,
        sum(case when status=0 and pred=1 then 1 else 0 end) as FP,
        sum(case when status=1 and pred=0 then 1 else 0 end) as FN
    from eval;
quit;

data metrics;
    set metrics;

    Accuracy = (TP+TN)/(TP+TN+FP+FN);
    Precision = TP/(TP+FP);
    Recall = TP/(TP+FN);
    F1 = 2*(Precision*Recall)/(Precision+Recall);
run;

proc print data=metrics;
    title "Evaluation Metrics";
run;


/*=============================================================
    [13] ROC Curve + AUC
=============================================================*/
ods output ROCAssociation=auc_data;
proc logistic data=train_data plots=roc;
    title "ROC Curve";
    model status = norm_avg norm_max WOE;
    roc;
run;

proc print data=auc_data;
    title "AUC Values";
run;


/*=============================================================
    [14] New Prediction
=============================================================*/
data new_patient;
    avg_fre=130;
    max_fre=200;
    avgfre_bin="Medium";
    WOE=0.2;
run;

data new_patient_norm;
    if _n_=1 then set norm;
    set new_patient;

    norm_avg=(avg_fre-min1)/(max1-min1);
    norm_max=(max_fre-min2)/(max2-min2);
run;

proc logistic inmodel=model_final;
    score data=new_patient_norm out=prediction;
run;

proc print data=prediction;
    title "New Patient Prediction";
run;


/*=============================================================
    [15] Export Key Results
=============================================================*/
ods excel file="/home/u64477338/sasuser.v94/sas project ss/parkinsons_results.xlsx"
    options(sheet_interval="proc");

proc print data=woe_final;
run;

proc print data=metrics;
run;

proc print data=auc_data;
run;

proc print data=prediction;
run;

proc export data=parkinsons_final
    outfile="/home/u64477338/sasuser.v94/sas project ss/parkinsons_final_dataset.xlsx"
    dbms=xlsx
    replace;
run;

ods excel close;

ods pdf close;
ods graphics off;

/*=============================================================
    END PROJECT
=============================================================*/
