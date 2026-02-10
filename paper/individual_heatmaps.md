# Heatmaps from individual rats
The following pages show heatmaps and associated analyses for each rat. Plots are limited to the high concentration (0.45 M) in deplete conditions to try to see which bodyparts, behavioural mesaurements, and analysis choices give the best discrimination between appetitive and aversive trials.

## Using Dionne's tracking
The 'best' bodyparts from visual inspection appear to be [r_ear, l_ear, head_base]. Some others could also be included although unclear whether this will provide better quantification. Could consider an unbiased way to select bodyparts, maybe by inputting them into a classifer and using other conditions as positive and negative training sets.
![Figure 1](../results/figure_snips_allbodyparts_deplete_45NaCl_individualrats.png)

<div style="page-break-before: always;"></div>

## Using movement metric
The following figure uses a movement metric that includes [r_ear, l_ear, head_base] (I think, need to double check). Gaussian smoothing with windowsize=10 is used. ~~Each bodypart is zscored to itself~~ and then summed. No baseline zscoring is done.
Each row is a different rat. Col. 1 is heatmap, col. 2 is median of infusion time, col.3 is time moving during infusion (threshold of 0.02), col.4 is histogram of movement metric for all bins/trials. Black line shows threshold for movement (0.02). The red triangle/dashed lines show transition point from dopamine data.

![Figure 2](../results/figure_snips_movement_deplete_45NaCl_individualrats.png)

show average heatmap here

<div style="page-break-before: always;"></div>

## Using movement metric with z-scoring
Z-scoring doesn't seem useful for movement metric (which makes sense, i.e. zero movement should be an absolute)

![Figure 3](../results/figure_snips_movement_z_deplete_45NaCl_individualrats.png)

## Show angular velocity
