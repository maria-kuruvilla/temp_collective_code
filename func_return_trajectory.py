# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 12:54:54 2020

@author: Maria Kuruvilla
"""

import trajectorytools as tt

      def trajectory(i, j , k):
    if j == 1:
        trajectories_file_path = 'G:/My Drive/CollectiveBehavior_Thermal_Experiments/Tracked/'+str(i)+'/' +str(j)+'/session_GS_'+str(j)+'_T_'+str(i)+'_'+str(k+1)+'/trajectories/trajectories.npy'
    else:
        trajectories_file_path = 'G:/My Drive/CollectiveBehavior_Thermal_Experiments/Tracked/'+str(i)+'/' +str(j)+'/session_GS_'+str(j)+'_T_'+str(i)+'_'+str(k+1)+'/trajectories_wo_gaps/trajectories_wo_gaps.npy'
    sigma_values = 1.5 #smoothing parameter
    tr = tt.Trajectories.from_idtrackerai(trajectories_file_path, center=True, smooth_params={'sigma': sigma_values}).normalise_by('body_length') # normalizing by body length
    tr.new_time_unit(tr.params['frame_rate'], 'seconds') # changing time unit to seconds
    return(tr) 