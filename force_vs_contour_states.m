function current_force=force_vs_contour_states(current_contour, current_state)
syms F_kT;
offset=10;%offset the stiff part(gfp label... not important)
linker_size=40;% the contour length of the flexible talin neck domain in nm (not important)
Rod_contour_length=[69.2 52.4 49.6 52.4 64 60.4 91.6 66.8 63.2 66.4 62.8];
uniform_rigid_size_unit=5;% mean domain size of folded domains
total_wormlike_contour=sum(Rod_contour_length(current_state~=1))+linker_size; 
total_N_rigid=sum(current_state==1);
if current_contour>(total_wormlike_contour+total_N_rigid*uniform_rigid_size_unit+offset)
    current_force=100;
else
wormlike_x_current=current_contour-offset-total_N_rigid.*(uniform_rigid_size_unit.*coth(F_kT.*uniform_rigid_size_unit)-1./F_kT);
eqn= F_kT*0.8==0.25.*(1-wormlike_x_current./total_wormlike_contour)^(-2)+wormlike_x_current./total_wormlike_contour-0.25;
current_force = vpasolve(eqn,F_kT).*4.1;
end
end