# This is a piece of code in propensities, separated here for clarity
# We take the ASVs with significant propensities towards the different
# trajectories (convergent to class 1 at starting/final point, to class 2
# or divergent) and represent the sum of their relative abundances at
# starting and final times.

setwd(dirOut)
# --- extract the intersections
inters_list = venn(ASV_sig_prop.list,show.plot = F)

inters_names =  names(attributes(inters_list)$intersections)
#inters_name = inters_names[24] # debug

for(inters_name in inters_names){
  intersection = unlist(attributes(inters_list)$intersections[inters_name])
  if((length(intersection) < 20) & 
     (inters_name != "final_class1:final_class2")){next} # represent only those with at least 30 ASVs
  for(exp in experiment){
    if(exp == "starting"){
      ASV.table.freq = ASV.table.0D.freq
      tmp1_id = parent1_id
      tmp2_id = parent2_id
      tmpmix_id = parentmix_id
    }else{
      ASV.table.freq = ASV.table.7D.freq
      tmp1_id = class1_id
      tmp2_id = class2_id
      tmpmix_id = classmix_id
    }
    matched = match(colnames(ASV.table.freq), intersection)
    
    freq_core_in_class1 = rowSums(ASV.table.freq[tmp1_id,!is.na(matched)])
    names(freq_core_in_class1) = rownames(ASV.table.freq)[tmp1_id]
    freq_core_in_class2 = rowSums(ASV.table.freq[tmp2_id,!is.na(matched)])
    names(freq_core_in_class2) = rownames(ASV.table.freq)[tmp2_id]
    freq_core_in_classmix = rowSums(ASV.table.freq[tmpmix_id,!is.na(matched)])
    names(freq_core_in_classmix) = rownames(ASV.table.freq)[tmpmix_id]
    
    if(exp == "final"){
      names_parents = child_to_parents[names(freq_core_in_class1)]
      freq_core = aggregate(freq_core_in_class1, list(names_parents), FUN=mean) 
      freq_core_in_class1 = freq_core$x; names(freq_core_in_class1) = freq_core$Group.1
      names_parents = child_to_parents[names(freq_core_in_class2)]
      freq_core = aggregate(freq_core_in_class2, list(names_parents), FUN=mean) 
      freq_core_in_class2 = freq_core$x; names(freq_core_in_class2) = freq_core$Group.1
      names_parents = child_to_parents[names(freq_core_in_classmix)]
      freq_core = aggregate(freq_core_in_classmix, list(names_parents), FUN=mean) 
      freq_core_in_classmix = freq_core$x; names(freq_core_in_classmix) = freq_core$Group.1
    }
    
    times = length(freq_core_in_class1)
    freq_core_in_class1 = data.frame(freq_core_in_class1,
                                     rep("convergent, class 1", times = times),
                                     rep(exp, times = times))
    times = length(freq_core_in_class2)
    colnames(freq_core_in_class1) = c("freq","Trajectory","Experiment")
    freq_core_in_class2 = data.frame(freq_core_in_class2,
                                     rep("convergent, class 2", times = times),
                                     rep(exp, times = times))
    colnames(freq_core_in_class2) = c("freq","Trajectory","Experiment")
    times = length(freq_core_in_classmix)
    freq_core_in_classmix = data.frame(freq_core_in_classmix,
                                       rep("divergent", times = times),
                                       rep(exp, times = times))
    colnames(freq_core_in_classmix) = c("freq","Trajectory","Experiment")
    freq_core.df = rbind(freq_core_in_class1,freq_core_in_class2,freq_core_in_classmix)
    if(exp == "starting"){
      freq_core.df.long = freq_core.df
      colnames(freq_core.df) = paste(colnames(freq_core.df),exp,sep = "_")
      freq_core.df.wide = freq_core.df
    }else{
      freq_core.df.long = rbind(freq_core.df.long, freq_core.df)
      colnames(freq_core.df) = paste(colnames(freq_core.df),exp,sep = "_")
      freq_core.df.wide = cbind(freq_core.df.wide,freq_core.df)
    }
  }
  
  # ... create plot
  my_palette = c("convergent, class 1"="chocolate4",
                 "convergent, class 2" = "chartreuse", "divergent"="blue")
  title = gsub(inters_name, pattern = ":", replacement = " and ")
  title = gsub(title, pattern = "_", replacement = " ") 
  plot_title = gsub(inters_name, pattern = ":", replacement = "Vs")
  plotOut = paste0("Plot_Frequency_",plot_title,"_vs_trajectory.pdf")
  pdf(plotOut, width = 10)
  gg = ggplot(freq_core.df.wide,
              aes(x = freq_starting, y = freq_final, 
                  color = Trajectory_starting))+
    geom_point(alpha = 0.5, size = 4) + #, 
    #           position=position_jitter(height=.3, width=.3))+
    #geom_hline(yintercept = 0.1)+
    #geom_line()+
    geom_abline(slope = 1, intercept = 0,linetype='dotted')+
    #scale_y_log10()+
    xlab("Sum of ASVs rel. abundances (starting)")+
    ylab("<Sum of ASVs rel. abundances> (final)")+
    ggtitle("ASVs with significant positive propensities at:",
            subtitle = title)+
    scale_color_manual(values = my_palette)+
    labs(color = "Trajectories")+
    theme_bw() +
    theme(title = element_text(size = 18),
          axis.title = element_text(size = 22),
          axis.text =  element_text(size = 18),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 22))#,
          #legend.position = "inside",
          #legend.position.inside = c(0.8,0.8))
  if(inters_name == "starting_class1:final_class1"){
    gg = gg + geom_hline(yintercept = 0.125, color = "black", linetype = "longdash") +
      geom_vline(xintercept = 0.125, linetype = "longdash", color = "black")
    #cat("debug1 \n")
    freq_core.df.wide.caseA = freq_core.df.wide
    colnames(freq_core.df.wide.caseA)=paste0(colnames(freq_core.df.wide),"A")
  }else if(inters_name == "starting_class1:starting_class2:starting_classmix:final_class1:final_class2:final_classmix"){
    gg = gg + geom_hline(yintercept = 0.2, color = "black", linetype = "longdash") +
      geom_vline(xintercept = 0.2, linetype = "longdash", color = "black")
    freq_core.df.wide.caseB = freq_core.df.wide
    colnames(freq_core.df.wide.caseB)=paste0(colnames(freq_core.df.wide),"B")
  }
  print(gg)
  dev.off()
}

# --- Plot two specific groups, those with propensity with trajectories
#     converging to class 1 and those generalists.
freq_core.df.wide.case = cbind(freq_core.df.wide.caseA, freq_core.df.wide.caseB)

xlab0 = "Sum of ASVs rel. abundances (S&F1)"
ylab0 = "Sum of ASVs rel. abundances (cosmopolitan)"
for(exp in experiment){
  my_palette = c("convergent, class 1"="chocolate4",
                 "convergent, class 2" = "chartreuse", "divergent"="blue")
  plotOut = paste0("Plot_Frequency_",exp,"_PropTrajectory1VsCosmopolitan_vs_trajectory.pdf")
  pdf(plotOut, width = 10)
  if(exp == "starting"){
    gg = ggplot(freq_core.df.wide.case,
                aes(x = freq_startingA, y = freq_startingB, 
                    color = Trajectory_startingA))
    title = "Starting communities"
    xlab = xlab0
    ylab = ylab0
  }else{
    gg = ggplot(freq_core.df.wide.case,
                aes(x = freq_finalA, y = freq_finalB, 
                    color = Trajectory_finalA))
    title = ""#Final communities, mean across replicates"
    xlab = paste0("<",xlab0,">")
    ylab = paste0("<",ylab0,">")
  }
  
  gg = gg + geom_point(alpha = 0.5, size = 4) + #, 
    geom_abline(slope = 1, intercept = 0,linetype='dotted')+
    xlab(xlab)+ ylab(ylab)+ ggtitle(title)+
    scale_color_manual(values = my_palette)+
    labs(color = "Trajectories")+
    theme_bw() +
    theme(title = element_text(size = 18),
          axis.title = element_text(size = 22),
          axis.text =  element_text(size = 18),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 22))
  print(gg)
  dev.off()
}
