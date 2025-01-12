
function [dot_bkg_L,dot_ret_L,dot_head_L]=stim_map(dot_x,dot_y,bkg_x,bkg_y,fix_x,fix_y,trial_spike,params)
dot_start=params.dot_start;
dot_end=params.dot_end;
stim_center=params.stim_center;
stim_dist=[81,61];
dot_bkg_L=zeros(max(trial_spike),stim_dist(1),stim_dist(2));
dot_ret_L=zeros(max(trial_spike),stim_dist(1),stim_dist(2));
dot_head_L=zeros(max(trial_spike),stim_dist(1),stim_dist(2));

for i=trial_spike
    %sti_frame
    te_frame=zeros(stim_dist(1),stim_dist(2));
    te_retino=zeros(stim_dist(1),stim_dist(2));
    te_head=zeros(stim_dist(1),stim_dist(2));
    for l=dot_start:dot_end        %change here for later version
        
        dot_bkg_pos_x=stim_center(1)+dot_x(i,l)-bkg_x(i);
        dot_bkg_pos_y=stim_center(2)+dot_y(i,l)-bkg_y(i);
        
        dot_ret_pos_x=stim_center(1)+dot_x(i,l)-fix_x(i);
        dot_ret_pos_y=stim_center(2)+dot_y(i,l)-fix_y(i);
        
        dot_head_pos_x=stim_center(1)+dot_x(i,l);
        dot_head_pos_y=stim_center(2)+dot_y(i,l);
        
        if ((dot_bkg_pos_x>0)&&(dot_bkg_pos_x<=stim_dist(1)))
            if ((dot_bkg_pos_y>0)&&(dot_bkg_pos_y<=stim_dist(2)))
                te_frame(dot_bkg_pos_x,dot_bkg_pos_y)=1;
            end
        end
        if ((dot_ret_pos_x>0)&&(dot_ret_pos_x<=stim_dist(1)))
            if ((dot_ret_pos_y>0)&&(dot_ret_pos_y<=stim_dist(2)))
                te_retino(dot_ret_pos_x,dot_ret_pos_y)=1;
            end
        end
        if ((dot_head_pos_x>0)&&(dot_head_pos_x<=stim_dist(1)))
            if ((dot_head_pos_y>0)&&(dot_head_pos_y<=stim_dist(2)))
                te_head(dot_head_pos_x,dot_head_pos_y)=1;
            end
        end
    end
    dot_bkg_L(i,:,:)= te_frame;
    dot_ret_L(i,:,:)= te_retino;
    dot_head_L(i,:,:)=te_head;
end
