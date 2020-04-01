[good_solution,good_t,good_x] = beam_eq_projection_dav(161,-1,0,1,0,0.1,1,1,2,1); 
[solution81,t81,x81] = beam_eq_projection_dav(81,-1,0,1,0,0.1,1,1,2,1);
% [solution41,t41,x41] = beam_eq_projection_dav(41,-1,0,1,0,0.1,1,1,2,1);
% [solution21,t21,x21] = beam_eq_projection_dav(21,-1,0,1,0,0.1,1,1,2,1);

error81 = 0;
for(i=1:320/160:length(good_x))
    if i==1
        error81 = error81 + good_solution(i)-solution81(i);
    else
        error81 = error81 + good_solution(i)-solution81(i*(160/320));
    end
end

