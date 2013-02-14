function [ ratio ] = right_hand_side_ratio(s, f, P_LF, X_LF, P_DG, X_DG)
    P_LF_at_fs = interp1(X_LF, P_LF, f, 'spline');
    
    P_DG_at_s = interp1(X_DG, P_DG, s, 'spline');
    
    ratio = P_DG_at_s / P_LF_at_fs;
end