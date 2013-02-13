function [ ratio ] = right_hand_side_ratio(s, fs, P_LF, X_LF, P_DG, X_DG)
    %prob_dg = interp1(pdg_s,pdg,s,'spline');
    P_LF_at_fs = interp1(X_LF, P_LF, fs, 'spline');
    
    P_DG_at_s = interp1(X_DG, P_DG, s, 'spline');
    
    ratio = P_DG_at_s / P_LF_at_fs;
end