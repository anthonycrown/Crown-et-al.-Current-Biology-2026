
#inputs a list and a smooting window, adds values for smoothed angular velocity for each item
smooth_velocity_values_list = function(n,m){
  library(smoother)
  
  input.list = n
  smoothing.frames = m
  
  for(i in c(1:length(input.list)))
  {
    curr.exp = input.list[[i]]
    
    curr.vel = curr.exp$Velocity
    curr.ang.vel = curr.exp$`Angular Velocity`
    
    #sometimes we have to transpose angular/trans vel matrices.
    #guess they get loaded in wrong sometimes?
    if(dim(curr.vel)[1] > dim(curr.vel)[2])
    {
      curr.vel = t(curr.vel)
      curr.ang.vel = t(curr.ang.vel)
      
      curr.exp$Velocity = curr.vel
      curr.exp$`Angular Velocity` = curr.ang.vel
    }
    
    n.flies = nrow(curr.vel)
    n.frames = ncol(curr.vel)
    
    #angular veloctiy:
    smoothed.ang = as.data.frame(matrix(NA,nrow = n.flies,ncol = n.frames))
    
    for(j in c(1:n.flies))
    {
      smoothing.vals = as.numeric(curr.vel[j,])
      smoothed.ang[j,] = smth.gaussian(smoothing.vals,window=smoothing.frames)
    }
    
    curr.exp$`Smoothed Translational Velocity` = smoothed.ang
    
    #translational velocity:
    
    smoothed.trans = as.data.frame(matrix(NA,nrow = n.flies,ncol = n.frames))
    
    for(j in c(1:n.flies))
    {
      smoothing.vals = as.numeric(curr.ang.vel[j,])
      smoothed.trans[j,] = smth.gaussian(smoothing.vals,window = smoothing.frames)
    }
    
    curr.exp$`Smoothed Angular Velocity` = smoothed.trans
    
    input.list[[i]] = curr.exp
    
  }
  
  return(input.list)
  
}



extract_walk_bouts_list = function(n,m){
  input.list = n
  bout.threshold = m
  
  
  for(i in c(1:length(input.list)))
  {
    
    curr.exp = input.list[[i]]
    
    curr.vel = curr.exp$`Smoothed Translational Velocity`
    curr.ang.vel = curr.exp$`Smoothed Angular Velocity`
    
    n.flies = nrow(curr.vel)
    n.frames = ncol(curr.vel)
    
    bouts.df = as.data.frame(matrix(NA,nrow=n.flies,ncol=n.frames))
    
    for (j in c(1:n.flies)) {
      
      current.bout = 1
      tracking.bout = FALSE
      tracking.cutoff = bout.threshold
      
      trial.trans = as.numeric(curr.vel[j,])
      
      bouts <- c()
      
      for (k in c(1:length(trial.trans))) {
        curr.trans <- trial.trans[k]
        
        bouts[k] <- NA
        if(is.na(curr.trans))
        {
          curr.trans = -1
        }
        
        if(curr.trans > tracking.cutoff) {
          tracking.bout <- TRUE
          bouts[k] <- current.bout
        } else{
          if(tracking.bout == TRUE)
          {
            
            tracking.bout = FALSE
            current.bout = current.bout + 1
          }
        }
        
      }
      
      bouts.df[j,] = as.character(bouts)
      
    }
    
    curr.exp$`Walking bouts` = bouts.df
    
    input.list[[i]] = curr.exp
  }
  
  return(input.list)
}


#takes a list of kinematics experiments and a vector describing what variables are in the names
#outputs a stats matrix for use in ggplot
#only used for non-opto data
stats_matrix_from_list = function(n,m,l){
  input.list = n
  exp.levels = m
  frames.to.analyze = l
  
  n.variables = length(exp.levels)
  
  input.list.levels = strsplit(names(input.list)," ")
  
  trans.vels = c()
  #trans.vels.smoothed = c()
  ang.vels = c()
  pearson.cor = c()
  n.bouts = c()
  bout.length = c()
  bout.distance = c()
  bout.trans.vel = c()
  bout.ang.vel = c()
  bout.pearson.cor = c()
  
  n.flies.total = c()
  
  for(i in c(1:length(input.list))){
    curr.exp = input.list[[i]]
    
    curr.vel = curr.exp$`Smoothed Translational Velocity`
    curr.ang.vel = curr.exp$`Smoothed Angular Velocity`
    curr.bouts = curr.exp$`Walking bouts`
    
    n.flies = nrow(curr.bouts)
    n.flies.total[i] = n.flies
    
    for(j in c(1:n.flies)){
      fly.trans.vel = as.numeric(curr.vel[j,c(1:frames.to.analyze)])
      fly.ang.vel = as.numeric(curr.ang.vel[j,c(1:frames.to.analyze)])
      fly.bouts = as.numeric(curr.bouts[j,c(1:frames.to.analyze)])
      
      
      trans.vels = append(trans.vels, mean(fly.trans.vel,na.rm=TRUE))
      ang.vels = append(ang.vels,mean(fly.ang.vel,na.rm=TRUE))
      pearson.cor = append(pearson.cor,cor(fly.trans.vel,fly.ang.vel,method = "pearson",use = "complete.obs"))
      n.bouts = append(n.bouts,max(fly.bouts,na.rm = TRUE))
      bout.length = append(bout.length,mean(table(fly.bouts)))
      
      fly.trans.vel.bout = fly.trans.vel[-which(is.na(fly.bouts))]
      fly.ang.vel.bout = fly.ang.vel[-which(is.na(fly.bouts))]
      
      bout.trans.vel = append(bout.trans.vel, mean(fly.trans.vel.bout))
      bout.ang.vel = append(bout.ang.vel, mean(fly.ang.vel.bout))
      #bout.pearson.cor = append(bout.pearson.cor,cor(fly.trans.vel,fly.ang.vel,method = "pearson", use = "complete"))
      bout.pearson.cor = append(bout.pearson.cor,cov(fly.trans.vel,fly.ang.vel,method = "pearson", use = "complete"))
      
    }
    
    
  }
  
  final.matrix = data.frame(
    mean_translational_velocity = trans.vels,
    mean_angular_velocity = ang.vels,
    pearson_cor = pearson.cor,
    n_bouts = n.bouts,
    bout_duration = bout.length,
    mean_bout_translational_velocity = bout.trans.vel,
    mean_bout_angular_velocity = bout.ang.vel,
    bout_pearson_cor = bout.pearson.cor
  )
  
  #adding in varibles
  start.col = ncol(final.matrix)+1
  end.col = start.col + n.variables-1
  cols = c(start.col:end.col)
  
  for(i in c(1:n.variables)){
    final.matrix[,cols[i]] = rep(  unlist(lapply(input.list.levels, function(n){n[i]})), n.flies.total)
  }
  
  colnames(final.matrix)[cols] = exp.levels
  
  return(final.matrix)
  
}



#takes a given exp (list format) and generates df with timecourse data
#specify which fly (use 1 if single fly trials)
timecourse.data.from.list =  function(experiment,fly,variables){
  curr.exp = experiment
  
  
  x.curr = curr.exp$`X mm`
  if(dim(x.curr)[1] > dim(x.curr)[2]) x.curr = t(x.curr)
  x.curr = as.numeric(x.curr[fly,])
  
  y.curr = curr.exp$`Y mm`
  if(dim(y.curr)[1] > dim(y.curr)[2]) y.curr = t(y.curr)
  y.curr = as.numeric(y.curr[fly,])
  
  ori.curr = curr.exp$Ori
  if(dim(ori.curr)[1] > dim(ori.curr)[2]) ori.curr = t(ori.curr)
  ori.curr = as.numeric(ori.curr[fly,])
  
  ang.vel = as.numeric(curr.exp$`Smoothed Angular Velocity`[fly,])
  #ang.vel = as.numeric(curr.exp$`Angular Velocity`[fly,])
  ang.vel = c(NA,ang.vel)
  
  trans.vel = as.numeric(curr.exp$`Smoothed Translational Velocity`[fly,])
  trans.vel = c(NA,trans.vel)
  
  walk.bouts = as.character( curr.exp$`Walking bouts`[fly,])
  walk.bouts = c(NA,walk.bouts)
  
  time = c(1:length(x.curr))/curr.exp$FPS
  
  final.data = data.frame(x = x.curr,
                          y = y.curr,
                          orientation = ori.curr,
                          time = time,
                          angular_velocity = ang.vel,
                          translational_velocity = trans.vel,
                          walking_bout_number = walk.bouts)
  
  #adding in variables
  start.col = ncol(final.data)
  end.col = start.col + length(variables) - 1
  
  replace.cols = c(start.col:end.col)
  
  #curr.variable.values = strsplit(names(curr.exp)," ")[[1]]
  
  for(i in c(1:length(variables)))
  {
    final.data[,replace.cols[i]] = rep( variables[i],nrow(final.data))
  }
  
  colnames(final.data)[replace.cols] = variables 
  
  return(final.data)
  
}





#same as non-opto version but relies on "Phase"
#Assuming single fly trials
#outputs a stats matrix for use in ggplot
stats_matrix_from_list_opto_single = function(n,m){
  input.list = n
  exp.levels = m
  
  n.variables = length(exp.levels)
  
  input.list.levels = strsplit(names(input.list), " ")
  
  trans.vels.on = c()
  trans.vels.off = c()
  delta.trans.vels = c()
  trans.vels.p.values = c()
  
  ang.vels.on = c()
  ang.vels.off = c()
  delta.ang.vels = c()
  ang.vel.p.values = c()
  
  abs.ang.vels.on = c()
  abs.ang.vels.off = c()
  abs.delta.ang.vels = c()
  abs.ang.vel.p.values = c()
  
  pearson.cor.on = c()
  pearson.cor.off = c()
  delta.pearson.cor = c()
  
  for (i in c(1:length(input.list))) {
    curr.exp = input.list[[i]]
    
    curr.vel =  as.numeric(curr.exp$`Smoothed Translational Velocity`[1, ])
    curr.vel = curr.vel / max(curr.vel, na.rm = TRUE)
    curr.ang.vel = as.numeric(curr.exp$`Smoothed Angular Velocity`[1, ])
    curr.phase = as.numeric(curr.exp$Phase[, 1])
    curr.phase = curr.phase[-1]
    
    curr.on.frames = curr.exp$`On frames`
    curr.off.frames = curr.exp$`Off frames`
    
    #trimming beginning and endings of experiments
    curr.stim.bout.length = curr.off.frames[1, 1] - curr.on.frames[1, 1]
    
    start.frame = curr.on.frames[1, 1] - curr.stim.bout.length
    if (start.frame <= 0)
      start.frame = 1
    end.frame = curr.off.frames[3, 1] + curr.stim.bout.length
    
    curr.vel = curr.vel[start.frame:end.frame]
    curr.ang.vel = curr.ang.vel[start.frame:end.frame]
    curr.ang.vel.abs = abs(curr.ang.vel)
    curr.phase = curr.phase[start.frame:end.frame]
    
    #splitting by stimulus phase
    curr.vel.off = curr.vel[which(curr.phase == 0)]
    curr.vel.on = curr.vel[which(curr.phase == 1)]
    
    curr.ang.vel.off = curr.ang.vel[which(curr.phase == 0)]
    curr.ang.vel.on = curr.ang.vel[which(curr.phase == 1)]
    
    curr.ang.abs.vel.off = curr.ang.vel.abs[which(curr.phase == 0)]
    curr.ang.abs.vel.on = curr.ang.vel.abs[which(curr.phase == 1)]
    
    #calculating stats
    vel.on.mean = mean(curr.vel.on, na.rm = TRUE)
    vel.off.mean = mean(curr.vel.off, na.rm = TRUE)
    vel.test = wilcox.test(curr.vel.on, curr.vel.off)
    
    ang.vel.on.mean = mean(curr.ang.vel.on, na.rm = TRUE)
    ang.vel.off.mean = mean(curr.ang.vel.off, na.rm = TRUE)
    ang.vel.test = wilcox.test(curr.ang.vel.off, curr.ang.vel.on)
    
    ang.vel.abs.on.mean = mean(curr.ang.abs.vel.on, na.rm = TRUE)
    ang.vel.abs.off.mean = mean(curr.ang.abs.vel.off, na.rm = TRUE)
    ang.vel.abs.test = wilcox.test(curr.ang.abs.vel.off, curr.ang.abs.vel.on)
    
    on.cor = cor(curr.vel.on,
                 curr.ang.vel.on,
                 method = "pearson",
                 use = "complete.obs")
    off.cor = cor(curr.vel.off,
                  curr.ang.vel.off,
                  method = "pearson",
                  use = "complete.obs")
    
    #filling data columns with above values
    trans.vels.on = append(trans.vels.on, vel.on.mean)
    trans.vels.off = append(trans.vels.off, vel.off.mean)
    delta.trans.vels = append(delta.trans.vels, (vel.on.mean - vel.off.mean))
    trans.vels.p.values = append(trans.vels.p.values, vel.test$p.value)
    
    ang.vels.on = append(ang.vels.on, ang.vel.on.mean)
    ang.vels.off = append(ang.vels.off, ang.vel.off.mean)
    delta.ang.vels = append(delta.ang.vels, (ang.vel.on.mean - ang.vel.off.mean))
    ang.vel.p.values = append(ang.vel.p.values, ang.vel.test$p.value)
    
    abs.ang.vels.on = append(abs.ang.vels.on, ang.vel.abs.on.mean)
    abs.ang.vels.off = append(abs.ang.vels.off, ang.vel.abs.off.mean)
    abs.delta.ang.vels = append(abs.delta.ang.vels,
                                (ang.vel.abs.on.mean - ang.vel.abs.off.mean))
    abs.ang.vel.p.values = append(abs.ang.vel.p.values, ang.vel.abs.test$p.value)
    
    pearson.cor.on = append(pearson.cor.on, on.cor)
    pearson.cor.off = append(pearson.cor.off, off.cor)
    delta.pearson.cor = append(delta.pearson.cor, (on.cor - off.cor))
  }
  
  
  final.matrix = data.frame(
    off_translational_velocity = trans.vels.off,
    on_translational_velocity = trans.vels.on,
    delta_translational_velocity = delta.trans.vels,
    translational_velocity_p.value = trans.vels.p.values,
    off_angular_velocity = ang.vels.off,
    on_angular_velocity = ang.vels.on,
    delta_angular_velocity = delta.ang.vels,
    angular_velocity_p_value = ang.vel.p.values,
    off_absolute_angular_velocity = abs.ang.vels.off,
    on_absolute_angular_velocity = abs.ang.vels.on,
    delta_absolute_angular_velocity = abs.delta.ang.vels,
    absolute_angular_velocity_p_value = abs.ang.vel.p.values,
    off_pearson_correlation = pearson.cor.off,
    on_pearson_correlation = pearson.cor.on,
    delta_pearson_correlation = delta.pearson.cor
  )
  
  #adding in varibles
  start.col = ncol(final.matrix) + 1
  end.col = start.col + n.variables - 1
  cols = c(start.col:end.col)
  
  for (i in c(1:n.variables)) {
    final.matrix[, cols[i]] = unlist(lapply(input.list.levels, function(n) {
      n[i]
    }))
  }
  
  colnames(final.matrix)[cols] = exp.levels
  
  return(final.matrix)
  
}

#gets the angle between x and y positions from successive time points
#used for angular motion that doesn't involve rotating the body
get_displacement_angles = function(input.list){
  
  for(i in c(1:length(input.list))){
    
    curr.exp = input.list[[i]]
    
    curr.x = curr.exp$X
    curr.y = curr.exp$Y
    curr.ori = curr.exp$Ori
    
    disp.ang = c(NA)
    
    for(j in c(2:(length(curr.x))))
    {
      x1 = as.numeric(curr.x[,(j-1)])
      y1 = as.numeric(curr.y[,(j-1)])
      x2 = as.numeric(curr.x[,j])
      y2 = as.numeric(curr.y[,j])
      
      dX = x2 - x1 
      dY = y2 - y1
      
      curr.ang = atan2(dY,dX)
      
      curr.ang = curr.ang - as.numeric(curr.ori[,j-1])
      
      if(curr.ang > pi){
        curr.ang = curr.ang - 2*pi
      }
      
      if(curr.ang < -pi){
        curr.ang = curr.ang + 2*pi
      }
      
      disp.ang[j] = curr.ang
    }
    
    curr.exp$`Angle of displacement` = disp.ang
    
    input.list[[i]] = curr.exp
  }
  
  return(input.list)
}


#assuming single flies
sep_opto_exp_list_to_bouts = function(input.list,level.names)
{
  
  #n.levels = length(input.list[[1]]$`On frames`)
  list.levels = strsplit(names(input.list)," ")
  n.levels = length(list.levels[[1]])
  
  n.bouts = nrow(input.list[[1]]$`On frames`)
  bout.length = input.list[[1]]$`Off frames`[1,] - input.list[[1]]$`On frames`[1,] - 1
  
  #phase = rep(0,bout.length*3)
  #phase[(bout.length+1):(bout.length*2)] = 1
  phase = rep(c("Pre","Light on","Post"),each = bout.length)
  
  #assuming 30 FPS recording
  time = c(1:length(phase))/32.63
  
  #output.list = list()
  output.df = as.data.frame(matrix(NA,ncol = (7+n.levels)))
  colnames(output.df) = c("Translational","Angular","Absolute_angular","Orientation","Phase","Bout","Time",level.names)
  
  for(i in c(1:length(input.list)))
  {
    curr.exp = input.list[[i]]
    curr.name = names(input.list)[i]
    
    curr.on.frames = curr.exp$`On frames`
    curr.off.frames = curr.exp$`Off frames`
    
    n.flies = nrow(curr.exp$`Smoothed Translational Velocity`)
    
    for(k in c(1:n.flies))
    {
      
      curr.trans.vel = as.numeric(curr.exp$`Smoothed Translational Velocity`[k,])
      curr.trans.vel = curr.trans.vel/max(curr.trans.vel,na.rm = TRUE)
      curr.ang.vel =  as.numeric(curr.exp$`Smoothed Angular Velocity`[k,])
      curr.orientation = as.numeric(curr.exp$Ori[k,])
      for(j in c(1:n.bouts))
      {
        #exp.df = as.data.frame(matrix(NA,nrow = length(phase),ncol=3))
        #colnames(exp.df) = c("Translational velocity","Angular velocity","Phase")
        
        start = curr.on.frames[j,1] - bout.length + 1
        if(start > 0)
        {
          end = curr.on.frames[j,1] + bout.length*2 #-1
          
          trans = curr.trans.vel[start:end]
          ang = curr.ang.vel[start:end]
          abs.ang = abs(ang)
          ori = curr.orientation[start:end]
          
          exp.df = data.frame("Translational" = trans, 
                              "Angular" = ang,
                              "Absolute_angular" = abs.ang,
                              "Orientation" = ori,
                              "Phase" = phase,
                              "Bout" = rep(j,length(phase)),
                              "Time"=time)
          
          for(k in c(1:n.levels))
          {
            exp.df[,(7+k)] =  rep(list.levels[[i]][k],length(phase))
          }
          
          colnames(exp.df)[8:(7+n.levels)] = level.names
          
          #curr.index = length(output.list) + 1
          #output.list[[curr.index]] = exp.df
          #names(output.list)[curr.index] = paste(curr.name,"Bout",j,sep = " ")
          
          output.df = rbind(output.df,exp.df)
        }
        
      }
      
    }
    
  }
  
  #return(exp.df)
  #return(output.list)
  output.df = output.df[-1,]
  
  return(output.df)
  
}


remove_na_trials = function(input.list){
  
  
  for(i in c(1:length(input.list))){
    
    curr.exp = input.list[[i]]
    curr.x = curr.exp$X
    curr.y = curr.exp$Y
    curr.ori = curr.exp$Ori
    curr.trans.vel = curr.exp$`Smoothed Translational Velocity`
    curr.ang.vel = curr.exp$`Smoothed Angular Velocity`
    curr.walking.bouts = curr.exp$`Walking bouts`
    
    to.remove = c()
    for( j in c(1:nrow(curr.trans.vel)))
    {
      if(length(which(!is.na(curr.trans.vel[j,])))==0)
      {
        
        to.remove = append(to.remove,j)
      }
    }
    
    if(length(to.remove > 0))
    {
      #print(names(input.list)[i])
      #print(to.remove)
      curr.x = curr.x[-to.remove,]
      curr.y =  curr.y[-to.remove,]
      curr.ori =  curr.ori[-to.remove,]
      curr.trans.vel =  curr.trans.vel[-to.remove,]
      curr.ang.vel =  curr.ang.vel[-to.remove,]
      curr.walking.bouts =  curr.walking.bouts[-to.remove,]
    }
    
    
    
    curr.exp$X = curr.x
    curr.exp$Y = curr.y
    curr.exp$Ori = curr.ori
    curr.exp$`Smoothed Translational Velocity` = curr.trans.vel
    curr.exp$`Smoothed Angular Velocity` = curr.ang.vel
    curr.exp$`Walking bouts` = curr.walking.bouts
    
    input.list[[i]] = curr.exp
  }
  
  return(input.list)
}



#takes in file path of csv files from tracking experiments (opto only)
#assumes 3 levels to naming...change later to be more modular
#also assuming 3 light bouts
opto_csv_to_list = function(file.directory,on.off.frames){
  
  
  library(imputeTS)
  library(readr)
  library(gridExtra)
  library(readxl)
  library(knitr)
  
  setwd(file.directory)
  
  files <- list.files()
  
  files.split <- strsplit(files, " ")
  
  tracked.files <- unlist(lapply(files.split,function(x){paste(x[1],x[2],x[3],x[4],sep = " ")}))
  
  tracked.files <- table(tracked.files)
  tracked.files <- names(tracked.files)
  
  n.exps <- length(tracked.files)
  
  on.off.df = read_excel(on.off.frames)
  
  output.list <- list()
  
  list.names <- tracked.files
  
  offset <- 0
  
  #reading in data
  
  for(i in c(1:n.exps))
  {
    
    input.file = tracked.files[i]
    
    on.off.index <- which(on.off.df[,1] == input.file)
    
    
    #frame of light turning on
    light.on = as.numeric(on.off.df[on.off.index,c(2,4,6)])
    
    #frame of light turning off
    light.off =  as.numeric(on.off.df[on.off.index,c(3,5,7)])
    
    
    #recording rate in FPS:
    fps = 32.63
    
    input.file.x = paste(input.file,"x.csv")
    input.file.y = paste(input.file,"y.csv")
    input.file.ori = paste(input.file,"ori.csv")
    
    if(file.exists(input.file.x)){
      
      x.coords <- read_csv(input.file.x,col_names = FALSE,show_col_types = FALSE)
      y.coords <- read_csv(input.file.y,col_names = FALSE,show_col_types = FALSE)
      ori.vals <- read_csv(input.file.ori,col_names = FALSE,show_col_types = FALSE)
      
      total.frames <- ncol(x.coords)
      
      phase = rep(0,total.frames)
      
      for(j in c(1:length(light.on)))
      {
        
        phase[c(light.on[j]:light.off[j])] = 1
        
      }
      
      list.contents = c(fps,NA,NA,NA,NA,NA,NA)
      names(list.contents) <- c("FPS","Phase","X","Y","Ori","On frames","Off frames")
      list.temp <- as.list(list.contents)
      list.temp[[2]] <- as.data.frame(phase)
      list.temp[[3]] <- as.data.frame(x.coords)
      list.temp[[4]] <- as.data.frame(y.coords)
      list.temp[[5]] <- as.data.frame(ori.vals)
      list.temp[[6]] <- as.data.frame(light.on)
      list.temp[[7]] <- as.data.frame(light.off)
      #change this, this is the variable you'll use to later concatenate all these files
      output.list[[i-offset]] <- list.temp
      
    }else{
      cat("Error: the following exp could not be found:")
      cat(input.file)
      list.names <- list.names[-i]
      offset <- offset + 1
    }
    
    
    
  }
  
  names(output.list) <- list.names
  
  #trans and ang vel values
  
  
  for (i in c(1:length(output.list)))
  {
    list.temp <- output.list[[i]]
    
    if(!"Velocity" %in% names(list.temp))
    {
      x.temp <- list.temp[["X"]]
      x.temp <- data.frame(lapply(x.temp,function(x){gsub("NaN",NA,x)}))
      x.temp <- data.frame(sapply(x.temp,function(x)as.numeric(x)))
      y.temp <- list.temp[["Y"]]
      y.temp <- data.frame(lapply(y.temp,function(x){gsub("NaN",NA,x)}))
      y.temp <- data.frame(sapply(y.temp,function(x)as.numeric(x)))
      ori.temp <- list.temp[["Ori"]]
      ori.temp <- data.frame(lapply(ori.temp,function(x){gsub("NaN",NA,x)}))
      ori.temp <- data.frame(sapply(ori.temp,function(x)as.numeric(x)))
      
      if(dim(x.temp)[1] > dim(x.temp)[2]){
        x.temp = t(x.temp)
        y.temp = t(y.temp)
        ori.temp = t(ori.temp)
      }
      
      
      list.temp[["X"]] <-x.temp
      list.temp[["Y"]] <- y.temp
      list.temp[["Ori"]] <- ori.temp
      
      #use this later to correct jumps in angular velocity when going full circle
      ori.temp.2 <- ori.temp
      
      for(l in (c(1:nrow(ori.temp.2))))
      {
        ori.temp.2.curr <- as.numeric(ori.temp.2[l,])
        ori.temp.2.curr[which(sign(ori.temp.2.curr) == -1)] <- ori.temp.2.curr[which(sign(ori.temp.2.curr) == -1)] + 2* pi
        ori.temp.2[l,] <- ori.temp.2.curr
      }
      
      
      n.pos <- ncol(x.temp)
      n.vel <- n.pos - 1
      
      fps.temp <- as.numeric(list.temp[["FPS"]])
      
      vel.temp <- as.data.frame(matrix(NA, ncol = n.vel, nrow = nrow(x.temp)))
      
      for(j in c(2:n.pos))
      {
        #calc the distance between current pos and pos at previous frame
        
        
        dist <- sqrt(
          (x.temp[,j] - x.temp[,j-1])^2 + (y.temp[,j] - y.temp[,j-1])^2
        )
        
        #velocity is distance multiplied by frames per second:
        vel.temp[,j - 1] = dist*fps.temp
        
      }
      
      ang.vel.temp <- as.data.frame(matrix(NA, ncol = n.vel, nrow = nrow(x.temp)))
      
      
      for(k in c(2:n.pos))
      {
        
        delta.theta.1 <- ori.temp[,k] - ori.temp[,k-1]
        
        delta.theta.2 <- ori.temp.2[,k] - ori.temp.2[,k-1]
        
        jump.instances <- which(abs(delta.theta.1) > abs(delta.theta.2))
        delta.theta.1[jump.instances] <- delta.theta.2[jump.instances]
        
        ang.vel.temp[,k-1] <- delta.theta.1*fps.temp
        
        
      }
      
      list.temp[["Velocity"]] <- vel.temp
      list.temp[["Angular Velocity"]] <- ang.vel.temp
      output.list[[i]] <- list.temp
    }
  }
  
  
  return(output.list)
  
  
}


tracker_csv_to_list = function(directory){
  
  setwd(directory)
  
  exp.files <- list.files()
  
  exp.files.split <- strsplit(exp.files, " ")
  
  
  exps  <- unlist(lapply(exp.files.split, function(x){
    if(length(x) == 1) {paste(x[1],sep=" ")}
    if(length(x) == 2) {paste(x[1],x[2],sep=" ")}
    if(length(x) == 3) {paste(x[1],x[2],x[3],sep=" ")}
    if(length(x) == 4) {paste(x[1],x[2],x[3],x[4],sep=" ")}
    if(length(x) == 5) {paste(x[1],x[2],x[3],x[4],x[5],sep=" ")}
    }))
  exps <- unique(exps)
  n.exps <-length(exps)
  
  
  output.list <- list()
  
  
  list.names <- exps
  
  offset <- 0
  
  for(i in c(1:n.exps))
  {
    
    input.file = exps[i]
    
    #recording rate in FPS:
    fps = 32.63
    
    input.file.x = paste(input.file,"x.csv")
    input.file.y = paste(input.file,"y.csv")
    input.file.ori = paste(input.file,"ori.csv")
    
    if(file.exists(input.file.x)){
      
      x.coords <- read_csv(input.file.x,col_names = FALSE,show_col_types = FALSE)
      y.coords <- read_csv(input.file.y,col_names = FALSE,show_col_types = FALSE)
      ori.vals <- read_csv(input.file.ori,col_names = FALSE,show_col_types = FALSE)
      
      total.frames <- ncol(x.coords)
      
      list.contents = c(fps,NA,NA,NA)
      names(list.contents) <- c("FPS","X","Y","Ori")
      list.temp <- as.list(list.contents)
      list.temp[[2]] <- as.data.frame(x.coords)
      list.temp[[3]] <- as.data.frame(y.coords)
      list.temp[[4]] <- as.data.frame(ori.vals)
      
      #change this, this is the variable you'll use to later concatenate all these files
      output.list[[i-offset]] <- list.temp
      
    }else{
      cat("Error: the following exp could not be found:")
      cat(input.file)
      list.names <- list.names[-i]
      offset <- offset + 1
    }
    
    
    
  }
  
  names(output.list) <- list.names
  
  return(output.list)
  
}


#scaling videos to mm to compare between trials
#(accounts for variations in camera zoom that may have occurred)
#pass in csv of arena dimensions
convert_list_to_mm = function(input.lists,dimensions.path){
  
  
  dimensions = read_excel(dimensions.path)
  
  
  for (i in c(1:length(input.lists))){
    
    curr.exp = input.lists[[i]]
    curr.exp.dimensions = subset(dimensions, subset = Exp == names(input.lists)[i])
    
    if(nrow(curr.exp.dimensions) >0)
    {
      
      curr.exp$X = data.frame(lapply(curr.exp$X,function(x){gsub("NaN",NA,x)}))
      curr.exp$X = data.frame(sapply(curr.exp$X,function(x)as.numeric(x)))
      
      curr.exp$Y = data.frame(lapply(curr.exp$Y,function(x){gsub("NaN",NA,x)}))
      curr.exp$Y = data.frame(sapply(curr.exp$Y,function(x)as.numeric(x)))
      
      curr.exp$Ori = data.frame(lapply(curr.exp$Ori,function(x){gsub("NaN",NA,x)}))
      curr.exp$Ori = data.frame(sapply(curr.exp$Ori,function(x)as.numeric(x)))
      
      
      #curr.exp$`Center` = c(((curr.exp.dimensions$x1 + curr.exp.dimensions$x2)/2),
      #                      ((curr.exp.dimensions$y1 + curr.exp.dimensions$y2)/2))
      
      xConvert = abs(curr.exp.dimensions$x1 - curr.exp.dimensions$x2)/148
      yConvert = abs(curr.exp.dimensions$y1 - curr.exp.dimensions$y2)/148
      
      curr.exp$`Arena X1` = curr.exp.dimensions$x1/xConvert
      curr.exp$`Arena X2` = curr.exp.dimensions$x2/xConvert
      curr.exp$`Arena y1` = curr.exp.dimensions$y1/yConvert
      curr.exp$`Arena y2` = curr.exp.dimensions$y2/yConvert
      
      x.mm = curr.exp$X/xConvert
      y.mm = curr.exp$Y/yConvert
      
      if(dim(x.mm)[1] > dim(y.mm)[2]) x.mm = t(x.mm)
      if(dim(y.mm)[1] > dim(y.mm)[2]) y.mm = t(y.mm)
      
      curr.exp$`X mm` = x.mm
      curr.exp$`Y mm` = y.mm
      
      curr.exp$Center = c(((curr.exp.dimensions$x1 + curr.exp.dimensions$x2)/2)/xConvert,
                          ((curr.exp.dimensions$y1 + curr.exp.dimensions$y2)/2)/yConvert)
      
      input.lists[[i]] = curr.exp
      
    }
    
    
    
  }
  
  return(input.lists)
  
}

calculate_velocity_values = function(input.lists)
{
  
  for (i in c(1:length(input.lists)))
  {
    list.temp <- input.lists[[i]]
    x.temp <- list.temp[["X mm"]]
    x.temp <- data.frame(lapply(x.temp,function(x){gsub("NaN",NA,x)}))
    x.temp <- data.frame(sapply(x.temp,function(x)as.numeric(x)))
    y.temp <- list.temp[["Y mm"]]
    y.temp <- data.frame(lapply(y.temp,function(x){gsub("NaN",NA,x)}))
    y.temp <- data.frame(sapply(y.temp,function(x)as.numeric(x)))
    ori.temp <- list.temp[["Ori"]]
    ori.temp <- data.frame(lapply(ori.temp,function(x){gsub("NaN",NA,x)}))
    ori.temp <- data.frame(sapply(ori.temp,function(x)as.numeric(x)))
    
    if(dim(x.temp)[1] >  dim(x.temp)[2]) x.temp = t(x.temp)
    if(dim(y.temp)[1] >  dim(y.temp)[2]) y.temp = t(y.temp)
    if(dim(ori.temp)[1] >  dim(ori.temp)[2]) ori.temp = t(ori.temp)
    
    list.temp[["X mm"]] <-x.temp
    list.temp[["Y mm"]] <- y.temp
    list.temp[["Ori"]] <- ori.temp
    
    #use this later to correct jumps in angular velocity when going full circle
    ori.temp.2 <- ori.temp
    
    for(l in (c(1:nrow(ori.temp.2))))
    {
      ori.temp.2.curr <- as.numeric(ori.temp.2[l,])
      ori.temp.2.curr[which(sign(ori.temp.2.curr) == -1)] <- ori.temp.2.curr[which(sign(ori.temp.2.curr) == -1)] + 2* pi
      ori.temp.2[l,] <- ori.temp.2.curr
    }
    
    
    n.pos <- ncol(x.temp)
    n.vel <- n.pos - 1
    
    fps.temp <- as.numeric(list.temp[["FPS"]])
    
    vel.temp <- as.data.frame(matrix(NA, ncol = n.vel, nrow = nrow(x.temp)))
    
    for(j in c(2:n.pos))
    {
      #calc the distance between current pos and pos at previous frame
      
      
      dist <- sqrt(
        (x.temp[,j] - x.temp[,j-1])^2 + (y.temp[,j] - y.temp[,j-1])^2
      )
      
      #velocity is distance multiplied by frames per second:
      vel.temp[,j - 1] = dist*fps.temp
      
    }
    
    ang.vel.temp <- as.data.frame(matrix(NA, ncol = n.vel, nrow = nrow(x.temp)))
    
    
    for(k in c(2:n.pos))
    {
      
      delta.theta.1 <-  ori.temp[,k] - ori.temp[,k-1] 
      
      delta.theta.2 <-  ori.temp.2[,k] - ori.temp.2[,k-1]
      
      jump.instances <- which(abs(delta.theta.1) > abs(delta.theta.2))
      delta.theta.1[jump.instances] <- delta.theta.2[jump.instances]
      
      ang.vel.temp[,k-1] <- delta.theta.1*fps.temp
      
      
    }
    
    list.temp[["Velocity"]] <- vel.temp
    list.temp[["Angular Velocity"]] <- ang.vel.temp
    input.lists[[i]] <- list.temp
  }
  
  if(length(input.lists) == 1)
  {
    input.lists <- input.lists[[1]]
  }
  
  return(input.lists)
  
  
}


calculate_distance_to_center = function(input.lists){
  
  
  for(i in c(1:length(input.lists)))
  {
    
    curr.exp = input.lists[[i]]
    
    
    x.temp = curr.exp$`X mm`
    y.temp = curr.exp$`Y mm`
    
    center.x = curr.exp$`Center`[1]
    center.y = curr.exp$`Center`[2]
    
    x.comp = x.temp-center.x
    x.comp = x.comp * x.comp
    
    y.comp = y.temp - center.y
    y.comp = y.comp * y.comp
    
    dist = x.comp + y.comp
    dist = sqrt(dist)
    
    curr.exp$`Distance to center` = dist
    
    input.lists[[i]] = curr.exp
    
  }
  
  return(input.lists)
  
  
}


timecourse_df = function(input.lists,exp.levels,start.frame,end.frame){
  
  n.levels = length(exp.levels)
  
  timecourse.length = length(c(start.frame:end.frame))
  
  output.df = as.data.frame(matrix(NA,ncol=(8+n.levels)))
  colnames(output.df) = c("x","y","ori","ang_vel","trans_vel","dist_to_center","fly","time",exp.levels)
  
  for(i in c(1:length(input.lists)))
  {
    
    curr.exp = input.lists[[i]]
    
    x.temp = curr.exp$`X mm`
    if(dim(x.temp)[1] > dim(x.temp)[2]) x.temp = t(x.temp)
    
    y.temp = curr.exp$`Y mm`
    if(dim(y.temp)[1] > dim(y.temp)[2]) y.temp = t(y.temp)
    
    ori.temp = curr.exp$Ori
    if(dim(ori.temp)[1] > dim(ori.temp)[2]) ori.temp = t(ori.temp)
    
    ang.vel.temp = curr.exp$`Smoothed Angular Velocity`
    if(dim(ang.vel.temp)[1] > dim(ang.vel.temp)[2]) ang.vel.temp = t(ang.vel.temp)
    
    trans.vel.temp = curr.exp$`Smoothed Translational Velocity`
    if(dim(trans.vel.temp)[1] > dim(trans.vel.temp)[2]) trans.vel.temp = t(trans.vel.temp)
    
    dist.temp = curr.exp$`Distance to center`
    if(dim(dist.temp)[1] > dim(dist.temp)[2]) dist.temp = t(dist.temp)
    
    curr.name = names(input.lists)[i]
    curr.name = strsplit(curr.name," ")[[1]]
    
    n.flies = nrow(x.temp)
    
    if(!is.null(n.flies))
    {
      
      for(j in c(1:n.flies)){
        
        df.temp = data.frame(
          x = as.numeric(x.temp[j,c(start.frame:end.frame)]),
          y = as.numeric(x.temp[j,c(start.frame:end.frame)]),
          ori = as.numeric(ori.temp[j,c(start.frame:end.frame)]),
          ang_vel = as.numeric(ang.vel.temp[j,c(start.frame:end.frame)]),
          trans_vel = as.numeric(trans.vel.temp[j,c(start.frame:end.frame)]),
          dist_to_center = as.numeric(dist.temp[j,c(start.frame:end.frame)])
        )
        
        df.temp$fly = j
      }
      
    } else
    {
      
      df.temp = data.frame(
        x = as.numeric(x.temp)[c(start.frame:end.frame)],
        y = as.numeric(x.temp)[c(start.frame:end.frame)],
        ori = as.numeric(ori.temp)[c(start.frame:end.frame)],
        ang_vel = as.numeric(ang.vel.temp)[1,c(start.frame:end.frame)],
        trans_vel = as.numeric(trans.vel.temp)[1,c(start.frame:end.frame)],
        dist_to_center = as.numeric(dist.temp)[c(start.frame:end.frame)]
      )
      
      df.temp$fly = 1
    }
    
    df.temp$time = c(1:timecourse.length)/32.63
    
    start.col = ncol(df.temp)
    for(l in c(1:n.levels))
    {
      
      df.temp[,(start.col+l)] = rep(curr.name[l],nrow(df.temp))
      
    }
    
    colnames(df.temp)[9:ncol(df.temp)] = exp.levels
    
    output.df = rbind(output.df,df.temp)
    
    
  }
  
  output.df = output.df[-1,]
  
  return(output.df)
  
}


#using to check how our boundaries match the data

path_plot_with_bounds = function(experiment){
  
  x1 = experiment$`Arena X1`
  x2 = experiment$`Arena X2`
  y1 = experiment$`Arena y1`
  y2 = experiment$`Arena y2`
  
  centerX = experiment$Center[1]
  centerY = experiment$Center[2]
  
  x = as.numeric(experiment$`X mm`[1,])
  y = as.numeric(experiment$`Y mm`[1,])
  
  data = data.frame(x = x,
                    y = y)
  
  library(ggplot2)
  
  ggplot(data, aes(x = x, y = y))+
    geom_path()+
    geom_point(aes(x = centerX, y = centerY))+
    geom_rect(aes(xmin = x1,
                  xmax = x2,
                  ymin = y1,
                  ymax = y2),
              alpha = 0.2,
              color = "orange",
              fill = NA)+
    coord_fixed(ratio = 1)
  
  
}


normalize_velocity_values = function(input.lists){
  
  
  for (i in c(1:length(input.lists))) {
    
    curr.data = input.lists[[i]]$`Smoothed Translational Velocity`
    
    for(j in c(1:nrow(curr.data))){
      curr.data[j,] = curr.data[j,]/max(curr.data[j,],na.rm = TRUE)
    }
    
    input.lists[[i]]$`Smoothed Translational Velocity Normalized` = curr.data
  }
  
  return(input.lists)
  
}


rad_to_angle = function(x) {((x * 360/(2*pi)) +180)}

dist_heading_to_edge = function(input.list)
{
  
  
  
  x.temp = input.list$`X mm`
  y.temp = input.list$`Y mm`
  ori.temp = rad_to_angle(input.list$Ori)
  
  center.temp = input.list$Center
  dist.to.center.temp = input.list$Center
  
  x.temp = x.temp - dist.to.center.temp[1]
  y.temp = y.temp - dist.to.center.temp[2]
  
  y.temp = -y.temp
  
  n.frames = ncol(x.temp)
  n.flies = nrow(x.temp)
  
  output.df = as.data.frame(matrix(NA, ncol = n.frames,nrow = n.flies))
  
  for(j in c(1:length(n.flies))){
    
    a = rad_to_angle(atan2(y.temp[j,],x.temp[j,]))
    b = ori.temp[j,] - a
    b[which(b < 0 )] = b[which(b < 0 )] + 360
    
    c = c()
    c[which(b<90)] = 90-b
    
    
    #return()
    test.df = data.frame(x = x.temp[j,],y=y.temp[j,],color = b)
    
    
  }
  
  ggplot(test.df,aes(x=x,y=y,color=color))+
    geom_path()+
    scale_color_viridis()+
    coord_fixed(ratio = 1)
  
}


circle_plot = function(experiment){
  
  library(circular)
  
  ori = experiment$Ori[1,]
  
  data = data.frame(ori = ori)
  
  ori.circular = circular(ori,units = "radians")
  mean = mean(ori.circular)
  
  ggplot(data, aes(x = ori))+
    geom_vline(xintercept = mean)+
    geom_histogram()+
    coord_polar()
  
  #return(rayleigh.test(ori.circular))
  
}

opto.single.trial.line.plot = function(experiment,fly){
  
  on = experiment$`On frames`$light.on
  off = experiment$`Off frames`$light.off
  
  interval.length =  off[1] - on[1]
  
  start = on[1] - interval.length
  end = off[3] + interval.length
  
  on = on - start
  off = off - start
  
  trans.vel = as.numeric(experiment$`Smoothed Translational Velocity`[fly,c(start:end)])
  trans.vel = trans.vel /max(trans.vel,na.rm = TRUE)
  ang.vel = as.numeric(experiment$`Smoothed Angular Velocity`[fly,c(start:end)])
  #ang.vel = abs(ang.vel)
  ang.vel = ang.vel / max(ang.vel,na.rm = TRUE)
  time = c(1:length(trans.vel)) /32.63
  
  data = data.frame(trans = trans.vel, ang = ang.vel, time = time)
  rect = data.frame(on = on/32.63, off = off/32.63)
  
  
  plots = list()
  
  plots[[1]] = ggplot(data,aes(x= time,y=trans))+
    geom_rect(data=rect, inherit.aes=FALSE, 
              aes(xmin=on, xmax=off, ymin=0, 
                  ymax=1), 
              color="transparent", fill="#f194a3", alpha=0.3)+
    geom_line()+
    theme_minimal()+
    geom_hline(yintercept = 0,linetype = "dashed",alpha = 0.1,col= "black")+
    theme(axis.line.y = element_line(size=1,linetype = "solid"),
          axis.text.y = element_text(face="bold",size=15),
          axis.text.x = element_text(face="bold",size=15))
  
  plots[[2]] = ggplot(data,aes(x= time,y=ang))+
    geom_rect(data=rect, inherit.aes=FALSE, 
              aes(xmin=on, xmax=off, ymin=0, 
                  ymax=0.5), 
              color="transparent", fill="#f194a3", alpha=0.3)+
    geom_line()+
    theme_minimal()+
    geom_hline(yintercept = 0,linetype = "dashed",alpha = 0.1,col= "black")+
    theme(axis.line.y = element_line(size=1,linetype = "solid"),
          axis.text.y = element_text(face="bold",size=15),
          axis.text.x = element_text(face="bold",size=15))
  
  
  return(plots)
  
  
}


opto.single.trial.bout.line.plot = function(experiment,fly,stim_bout){
  
  on = experiment$`On frames`$light.on[stim_bout]
  off = experiment$`Off frames`$light.off[stim_bout]
  
  interval.length =  off - on
  
  start = on - interval.length
  end = off + interval.length
  
  on = on - start
  off = off - start
  
  trans.vel = as.numeric(experiment$`Smoothed Translational Velocity`[fly,c(start:end)])
  trans.vel = trans.vel /max(trans.vel,na.rm = TRUE)
  ang.vel = as.numeric(experiment$`Smoothed Angular Velocity`[fly,c(start:end)])
  #ang.vel = abs(ang.vel)
  #ang.vel = ang.vel / max(ang.vel,na.rm = TRUE)
  time = c(1:length(trans.vel)) /32.63
  
  data = data.frame(trans = trans.vel, ang = ang.vel, time = time)
  rect = data.frame(on = on/32.63, off = off/32.63)
  
  ymax.ang = max(abs(ang.vel),na.rm = TRUE)
  
  plots = list()
  
  plots[[1]] = ggplot(data,aes(x= time,y=trans))+
    geom_rect(data=rect, inherit.aes=FALSE, 
              aes(xmin=on, xmax=off, ymin=0, 
                  ymax=1), 
              color="transparent", fill="#f194a3", alpha=0.3)+
    geom_line()+
    theme_minimal()+
    geom_hline(yintercept = 0,linetype = "dashed",alpha = 0.1,col= "black")+
    theme(axis.line.y = element_line(size=1,linetype = "solid"),
          axis.text.y = element_text(face="bold",size=15),
          axis.text.x = element_text(face="bold",size=15))
  
  plots[[2]] = ggplot(data,aes(x= time,y=ang))+
    geom_rect(data=rect, inherit.aes=FALSE, 
              aes(xmin=on, xmax=off, ymin=-ymax.ang, 
                  ymax=ymax.ang), 
              color="transparent", fill="#f194a3", alpha=0.3)+
    geom_line()+
    theme_minimal()+
    geom_hline(yintercept = 0,linetype = "dashed",alpha = 0.1,col= "black")+
    theme(axis.line.y = element_line(size=1,linetype = "solid"),
          axis.text.y = element_text(face="bold",size=15),
          axis.text.x = element_text(face="bold",size=15))
  
  
  return(plots)
  
  
}

opto.single.trial.path.plot = function(experiment,fly){
  
  on = experiment$`On frames`$light.on
  off = experiment$`Off frames`$light.off
  
  interval.length =  off[1] - on[1]
  
  start = on[1] - interval.length
  end = off[3] + interval.length
  
  on = on - start
  off = off - start
  
  x = as.numeric(experiment$X[fly,c(start:end)])
  y = as.numeric(experiment$Y[fly,c(start:end)])
  
  phase = rep("off", length(x))
  phase[c(c(on[1]:off[1]),c(on[2]:off[2]),c(on[3]:off[3]))] = "on"
  
  data = data.frame(x= x, y= y,phase = phase)
  rect = data.frame(on = on/32.63, off = off/32.63)
  
  
  #plots = list()
  
  plot = ggplot(data,aes(x= x,y=y, color = phase))+
    geom_path()+
    theme_minimal()+
    scale_color_manual(values = c("black","#f194a3"))+
    theme(panel.grid = element_blank(),
          axis.line.y = element_line(size=1,linetype = "solid"),
          axis.text.y = element_text(face="bold",size=15),
          axis.text.x = element_text(face="bold",size=15))
  
  
  return(plot)
  
  
}

'%!in%' <- function(x,y)!('%in%'(x,y))

bouts_plot = function(experiment,fly,dist.cutoff,exp.length){
  
  curr.x = as.numeric(experiment$`X mm`[fly,])
  curr.x = curr.x[-1]
  curr.x = curr.x[c(0:exp.length)]
  curr.y = as.numeric(experiment$`Y mm`[fly,])
  curr.y = curr.y[-1]
  curr.y = curr.y[c(0:exp.length)]
  
  curr.dist = as.numeric(experiment$`Distance to center`[fly,])
  curr.dist = curr.dist[-1]
  curr.dist = curr.dist[c(0:exp.length)]
  
  curr.ang = as.numeric(experiment$`Smoothed Angular Velocity`[fly,])[c(0:exp.length)]
  #curr.ang = as.numeric(experiment$`Angular Velocity`[fly,])[c(0:exp.length)]
  curr.trans = as.numeric(experiment$`Smoothed Translational Velocity`[fly,])[c(0:exp.length)]
  
  curr.bouts = as.character(experiment$`Walking bouts`[fly,])
  curr.bouts = curr.bouts[c(0:exp.length)]
  
  time = c(1:exp.length)/32.63
  
  data.temp = data.frame(x = curr.x, y = curr.y, bouts = curr.bouts,trans = curr.trans,ang = curr.ang,dist = curr.dist,time = time)
  
  include = curr.bouts[which(curr.dist < dist.cutoff)]
  
  data.temp = subset(data.temp, subset = bouts %in% include)
  
  data.temp = subset(data.temp, subset = !is.na(bouts))
  
  final.bouts = names(table(data.temp$bouts))
  
  for(i in final.bouts)
  {
    
    curr.inds = which(data.temp$bouts == i)
    
    x.transformed = data.temp$x[curr.inds]
    y.transformed = data.temp$y[curr.inds]
    
    x.transformed = x.transformed - x.transformed[1]
    y.transformed = y.transformed - y.transformed[1]
    
    data.temp$x[curr.inds] = x.transformed
    data.temp$y[curr.inds] = y.transformed
    
    data.temp$time[curr.inds] = c(1:length(curr.inds))/32.63
    
  }
  
  plots = list()
  
  plots[[1]] = ggplot(data.temp, aes(x = x,y = y, color = bouts))+
    geom_path(show.legend = FALSE)+
    theme_minimal()+
    coord_fixed(ratio = 1)+
    theme(panel.grid = element_blank(),
          axis.line.y = element_line(size=1,linetype = "solid"),
          axis.line.x = element_line(size=1,linetype = "solid"),
          axis.ticks = element_line(size=1,linetype = "solid"),
          axis.text.y = element_text(face="bold",size=15),
          axis.text.x = element_text(face="bold",size=15))
  
  plots[[2]] = ggplot(data.temp, aes(x = ang, y = trans, color = bouts))+
    geom_vline(xintercept = (2 * pi/3), linetype = "dashed", color = "black", size = 1)+
    geom_point(size = 1, alpha = 0.3)+
    #geom_density_2d()+
    theme_minimal()+
    #stat_smooth(method = lm, aes(x = ang_vel, y = trans_vel))+
    theme(#panel.grid = element_blank(),
      axis.line.y = element_line(linewidth=1,linetype = "solid"),
      axis.line.x = element_line(linewidth=1,linetype = "solid"),
      axis.ticks = element_line(linewidth=1,linetype = "solid"),
      legend.title = element_blank(),
      legend.position = "none")
  
  plots[[3]] = ggplot(data.temp, aes(x = time,y=trans,color = bouts))+
    geom_line(size = 1)+
    theme_minimal()+
    theme(#panel.grid = element_blank(),
      axis.line.y = element_line(linewidth=1,linetype = "solid"),
      axis.line.x = element_line(linewidth=1,linetype = "solid"),
      axis.ticks = element_line(linewidth=1,linetype = "solid"),
      legend.title = element_blank(),
      legend.position = "none")
  
  
  
  plots[[4]] = ggplot(data.temp, aes(x = time,y=ang,color = bouts))+
    geom_hline(size = 1,yintercept = 2*pi/3,linetype="dashed",color = "black")+
    geom_line(size = 1)+
    theme_minimal()+
    theme(#panel.grid = element_blank(),
      axis.line.y = element_line(linewidth=1,linetype = "solid"),
      axis.line.x = element_line(linewidth=1,linetype = "solid"),
      axis.ticks = element_line(linewidth=1,linetype = "solid"),
      legend.title = element_blank(),
      legend.position = "none")+
    ylim(c(-10,10))
  
  #plots[[5]] = ggplot(data.temp, aes(x = time,y=dist,color = bouts))+
  #  geom_line(size = 1)+
  #  theme_minimal()+
  #  theme(#panel.grid = element_blank(),
  #    axis.line.y = element_line(linewidth=1,linetype = "solid"),
  #    axis.line.x = element_line(linewidth=1,linetype = "solid"),
  #    axis.ticks = element_line(linewidth=1,linetype = "solid"),
  #    legend.title = element_blank(),
  #    legend.position = "none")
  
  #plots[[6]] = ggplot(data.temp, aes(x = dist,y=ang,color = bouts))+
  #  geom_line(size = 1)+
  #  theme_minimal()+
  #  theme(#panel.grid = element_blank(),
  #    axis.line.y = element_line(linewidth=1,linetype = "solid"),
  #    axis.line.x = element_line(linewidth=1,linetype = "solid"),
  #    axis.ticks = element_line(linewidth=1,linetype = "solid"),
  #    legend.title = element_blank(),
  #    legend.position = "none")
  
  
  return(plots)  
  
}

timecourse.data.from.list.02 =  function(experiment,fly,variables){
  curr.exp = experiment
  
  
  x.curr = curr.exp$`X`
  if(dim(x.curr)[1] > dim(x.curr)[2]) x.curr = t(x.curr)
  x.curr = as.numeric(x.curr[fly,])
  
  y.curr = curr.exp$`Y`
  if(dim(y.curr)[1] > dim(y.curr)[2]) y.curr = t(y.curr)
  y.curr = as.numeric(y.curr[fly,])
  
  ori.curr = curr.exp$Ori
  if(dim(ori.curr)[1] > dim(ori.curr)[2]) ori.curr = t(ori.curr)
  ori.curr = as.numeric(ori.curr[fly,])
  
  ang.vel = as.numeric(curr.exp$`Smoothed Angular Velocity`[fly,])
  #ang.vel = as.numeric(curr.exp$`Angular Velocity`[fly,])
  ang.vel = c(NA,ang.vel)
  
  trans.vel = as.numeric(curr.exp$`Smoothed Translational Velocity`[fly,])
  trans.vel = c(NA,trans.vel)
  
  walk.bouts = as.character( curr.exp$`Walking bouts`[fly,])
  walk.bouts = c(NA,walk.bouts)
  
  time = c(1:length(x.curr))/curr.exp$FPS
  
  final.data = data.frame(x = x.curr,
                          y = y.curr,
                          orientation = ori.curr,
                          time = time,
                          angular_velocity = ang.vel,
                          translational_velocity = trans.vel,
                          walking_bout_number = walk.bouts)
  
  #adding in variables
  start.col = ncol(final.data)
  end.col = start.col + length(variables) - 1
  
  replace.cols = c(start.col:end.col)
  
  #curr.variable.values = strsplit(names(curr.exp)," ")[[1]]
  
  for(i in c(1:length(variables)))
  {
    final.data[,replace.cols[i]] = rep( variables[i],nrow(final.data))
  }
  
  colnames(final.data)[replace.cols] = variables 
  
  return(final.data)
  
}

bouts_matrix = function(input.lists,exp.levels,start.frame,end.frame,radius){
  
  library(tseries)
  
  n.levels = length(exp.levels)
  
  n.frames = length(c(start.frame:end.frame))
  
  time.temp = c(1:n.frames)/32.63
  
  plot.data = as.data.frame(matrix(NA, ncol = (9+n.levels)))
  colnames(plot.data) = c("nBouts","bout_length","bout_duration","mean_trans_vel","mean_ang_vel","curvature","fly","autocorrelation","ang_var",exp.levels)
  
  exp.names = strsplit( names(input.lists), " ")
  
  for(i in c(1:length(input.lists)))
  {
    
    curr.exp = input.lists[[i]]
    curr.name = exp.names[[i]]
    
    curr.bouts = curr.exp$`Walking bouts`
    curr.trans.vel = curr.exp$`Velocity`/32.63
    curr.trans.vel.smoothed = curr.exp$`Smoothed Translational Velocity`
    curr.ang.vel = curr.exp$`Angular Velocity`
    curr.dist = curr.exp$`Distance to center`
    
    if(!is.null(nrow(curr.bouts)))
    {
      
      for (j in c(1:nrow(curr.bouts))) {
        
        data.temp = data.frame(bouts = as.numeric(curr.bouts[j,c(start.frame:end.frame)]),
                               trans = as.numeric(curr.trans.vel[j,c(start.frame:end.frame)]),
                               smooth_trans = as.numeric(curr.trans.vel.smoothed[j,c(start.frame:end.frame)]),
                               ang = as.numeric(curr.ang.vel[j,c(start.frame:end.frame)]),
                               dist = as.numeric(curr.dist[j,c(start.frame:end.frame)]))
        
        #
        include.bouts = names(table(data.temp$bouts[which(data.temp$dist < radius)]))
        
        curr.n.bouts = length(include.bouts)
        curr.bout.duration = NA
        curr.bout.length = NA
        curr.bout.mean.trans.vel = NA
        curr.bout.mean.ang.vel = NA
        curr.bout.curvature = NA
        curr.bout.autocorrelation = NA
        curr.ang.var = NA
        
        if(curr.n.bouts > 0)
        {
          dists = c()
          trans.vels = c()
          ang.vels = c()
          curvatures = c()
          durations = c()
          autocorrelations = c()
          vars = c()
          
          for(k in include.bouts){
            
            data.temp.sub = subset(data.temp, subset = bouts == k)
            
            bout.dist = sum(data.temp.sub$trans,na.rm = TRUE)
            bout.trans.vel = mean(data.temp.sub$smooth_trans,na.rm = TRUE)
            bout.ang.vel = mean(data.temp.sub$ang,na.rm = TRUE)
            
            bout.curvature = sum(data.temp.sub$ang,na.rm = TRUE)/sum(data.temp.sub$smooth_trans,na.rm = TRUE)
            bout.curvature = mean(bout.curvature,na.rm = TRUE)
            
            bout.duration = nrow(data.temp.sub) / 32.63
            
            bout.acf = acf(as.numeric(data.temp.sub$ang),
                           plot = FALSE,
                           na.action = na.remove)[10]
            
            bout.acf = as.numeric(bout.acf$acf)
            
            bout.var = var(data.temp.sub$ang ,na.rm = TRUE)
            
            dists = append(dists,bout.dist)
            trans.vels = append(trans.vels,bout.trans.vel)
            ang.vels = append(ang.vels,bout.ang.vel)
            curvatures = append(curvatures,bout.curvature)
            durations = append(durations,bout.duration)
            autocorrelations = append(autocorrelations,bout.acf)
            vars = append(vars,bout.var)
            
          }
          
          curr.bout.length = mean(dists,na.rm = TRUE)
          curr.bout.mean.trans.vel = mean(trans.vels, na.rm = TRUE)
          curr.bout.mean.ang.vel = mean(ang.vels, na.rm = TRUE)
          curr.bout.curvature = mean(curvatures,na.rm = TRUE)
          
          curr.bout.duration = mean(durations, na.rm = TRUE)
          curr.bout.autocorrelation = mean(autocorrelations,na.rm = TRUE)
          curr.ang.var = mean(vars, na.rm = TRUE)
        }
        
        fly.data = data.frame(nBouts = curr.n.bouts,
                              bout_length = curr.bout.length,
                              bout_duration = curr.bout.duration,
                              mean_trans_vel = curr.bout.mean.trans.vel,
                              mean_ang_vel = curr.bout.mean.ang.vel,
                              curvature = curr.bout.curvature,
                              fly = j,
                              autocorrelation = curr.bout.autocorrelation,
                              ang_var = curr.ang.var)
        
        for (l in c(1:length(exp.levels))) {
          
          fly.data[,9+l] = curr.name[l]
          
        }
        
        colnames(fly.data)[c(10:(9+n.levels))] = exp.levels
        
        plot.data = rbind(plot.data,fly.data)
        
      }
      
    }
    
  }
  
  plot.data = plot.data[-1,]
  
  return(plot.data)
  
}



path_plot = function(experiment,start.frame,end.frame){
  
  x = as.numeric(experiment$`X mm`[1,])[c(start.frame:end.frame)]
  y = as.numeric(experiment$`Y mm`[1,])[c(start.frame:end.frame)]
  
  data = data.frame(x = x,
                    y = y)
  
  library(ggplot2)
  
  ggplot(data, aes(x = x, y = y))+
    geom_path()+
    coord_fixed(ratio = 1)+
    theme_minimal()+
    theme(panel.grid = element_blank(),
          axis.line.y = element_line(size=1,linetype = "solid"),
          axis.line.x = element_line(size=1,linetype = "solid"),
          axis.ticks = element_line(size=1,linetype = "solid"),
          axis.text.y = element_text(face="bold",size=15),
          axis.text.x = element_text(face="bold",size=15))
  
  
}


bouts_plot_corrected_heading = function(experiment,fly,dist.cutoff,exp.length,stim.speed){
  
  curr.x = as.numeric(experiment$`X mm`[fly,])
  curr.x = curr.x[-1]
  curr.x = curr.x[c(0:exp.length)]
  curr.y = as.numeric(experiment$`Y mm`[fly,])
  curr.y = curr.y[-1]
  curr.y = curr.y[c(0:exp.length)]
  curr.ori = as.numeric(experiment$Ori[fly,])
  curr.ori = curr.ori[-1]
  curr.ori = curr.ori[c(0:exp.length)]
  
  curr.dist = as.numeric(experiment$`Distance to center`[fly,])
  curr.dist = curr.dist[-1]
  curr.dist = curr.dist[c(0:exp.length)]
  
  curr.ang = as.numeric(experiment$`Smoothed Angular Velocity`[fly,])[c(0:exp.length)]
  #curr.ang = as.numeric(experiment$`Angular Velocity`[fly,])[c(0:exp.length)]
  #curr.trans = as.numeric(experiment$`Smoothed Translational Velocity`[fly,])[c(0:exp.length)]
  curr.trans = as.numeric(experiment$Velocity[fly,])[c(0:exp.length)]
  
  curr.bouts = as.character(experiment$`Walking bouts`[fly,])
  curr.bouts = curr.bouts[c(0:exp.length)]
  
  time = c(1:exp.length)/32.63
  
  data.temp = data.frame(x = curr.x, 
                         y = curr.y,
                         ori = curr.ori,
                         bouts = curr.bouts,
                         trans = curr.trans,
                         ang = curr.ang,
                         dist = curr.dist,
                         time = time)
  
  include = curr.bouts[which(curr.dist < dist.cutoff)]
  
  data.temp = subset(data.temp, subset = bouts %in% include)
  
  data.temp = subset(data.temp, subset = !is.na(bouts))
  
  final.bouts = names(table(data.temp$bouts))
  
  #new ang vel column which is the speed of the rotating stimulus 
  #i.e. the theoretical error between heading and goal
  
  data.temp$ang_corrected = data.temp$ang - stim.speed
  data.temp$ori_error_corrected = NA
  data.temp$x_error_corrected = NA
  data.temp$y_error_corrected = NA
  
  
  for(i in final.bouts)
  {
    
    curr.inds = which(data.temp$bouts == i)
    
    x.transformed = data.temp$x[curr.inds]
    y.transformed = data.temp$y[curr.inds]
    
    ang.corr.sub = as.numeric(data.temp$ang_corrected[curr.inds])
    #ang.corr.sub = as.numeric(data.temp$ang[curr.inds])
    trans.sub = data.temp$trans[curr.inds]
    
    x.transformed = x.transformed - x.transformed[1]
    y.transformed = y.transformed - y.transformed[1]
    
    data.temp$x[curr.inds] = x.transformed
    data.temp$y[curr.inds] = y.transformed
    
    data.temp$time[curr.inds] = c(1:length(curr.inds))/32.63
    
    #simulating walking trajectories with error-subtracted ang vel values
    
    theta.ang.corrected = c()
    x.ang.corrected = c()
    y.ang.corrected = c()
    
    starting.theta = as.numeric(data.temp$ori[curr.inds])
    starting.theta = starting.theta[1]
    
    theta.ang.corrected[1] = starting.theta
    x.ang.corrected[1] = 0
    y.ang.corrected[1] = 0
    
    for(j in c(2:length(curr.inds))){
      
      theta.ang.corrected[j] = theta.ang.corrected[j-1] + ang.corr.sub[j]
      x.ang.corrected[j] = trans.sub[j]*cos(theta.ang.corrected[j])+x.ang.corrected[j-1]
      y.ang.corrected[j] = trans.sub[j]*sin(theta.ang.corrected[j])+y.ang.corrected[j-1]
      
    }
    
    data.temp$ori_error_corrected[curr.inds] = theta.ang.corrected
    data.temp$x_error_corrected[curr.inds] = x.ang.corrected
    data.temp$y_error_corrected[curr.inds] = y.ang.corrected
    
  }
  
  plots = list()
  
  plots[[1]] = ggplot(data.temp, aes(x = x,y = y, color = bouts))+
    geom_path(show.legend = FALSE)+
    theme_minimal()+
    coord_fixed(ratio = 1)+
    theme(panel.grid = element_blank(),
          axis.line.y = element_line(size=1,linetype = "solid"),
          axis.line.x = element_line(size=1,linetype = "solid"),
          axis.ticks = element_line(size=1,linetype = "solid"),
          axis.text.y = element_text(face="bold",size=15),
          axis.text.x = element_text(face="bold",size=15))
  
  plots[[2]] = ggplot(data.temp, aes(x = x_error_corrected,y = y_error_corrected, color = bouts))+
    geom_path(show.legend = FALSE)+
    theme_minimal()+
    coord_fixed(ratio = 1)+
    theme(panel.grid = element_blank(),
          axis.line.y = element_line(size=1,linetype = "solid"),
          axis.line.x = element_line(size=1,linetype = "solid"),
          axis.ticks = element_line(size=1,linetype = "solid"),
          axis.text.y = element_text(face="bold",size=15),
          axis.text.x = element_text(face="bold",size=15))
  
  
  
  return(plots)  
  
}

per_bouts_matrix = function(input.lists,exp.levels,start.frame,end.frame,radius){
  
  library(tseries)
  
  n.levels = length(exp.levels)
  
  n.frames = length(c(start.frame:end.frame))
  
  time.temp = c(1:n.frames)/32.63
  
  plot.data = as.data.frame(matrix(NA, ncol = (9+n.levels)))
  colnames(plot.data) = c("nBouts","bout_length","bout_duration","mean_trans_vel","mean_ang_vel","curvature","fly","autocorrelation","ang_var",exp.levels)
  
  exp.names = strsplit( names(input.lists), " ")
  
  for(i in c(1:length(input.lists)))
  {
    
    curr.exp = input.lists[[i]]
    curr.name = exp.names[[i]]
    
    curr.bouts = curr.exp$`Walking bouts`
    curr.trans.vel = curr.exp$`Velocity`/32.63
    curr.trans.vel.smoothed = curr.exp$`Smoothed Translational Velocity`
    curr.ang.vel = curr.exp$`Smoothed Angular Velocity`
    curr.dist = curr.exp$`Distance to center`
    
    if(!is.null(nrow(curr.bouts)))
    {
      
      for (j in c(1:nrow(curr.bouts))) {
        
        data.temp = data.frame(bouts = as.numeric(curr.bouts[j,c(start.frame:end.frame)]),
                               trans = as.numeric(curr.trans.vel[j,c(start.frame:end.frame)]),
                               smooth_trans = as.numeric(curr.trans.vel.smoothed[j,c(start.frame:end.frame)]),
                               ang = as.numeric(curr.ang.vel[j,c(start.frame:end.frame)]),
                               dist = as.numeric(curr.dist[j,c(start.frame:end.frame)]))
        
        #
        include.bouts = names(table(data.temp$bouts[which(data.temp$dist < radius)]))
        
        curr.n.bouts = length(include.bouts)
        curr.bout.duration = NA
        curr.bout.length = NA
        curr.bout.mean.trans.vel = NA
        curr.bout.mean.ang.vel = NA
        curr.bout.curvature = NA
        curr.bout.autocorrelation = NA
        curr.bout.var = NA
        
        if(curr.n.bouts > 0)
        {
          dists = c()
          trans.vels = c()
          ang.vels = c()
          curvatures = c()
          durations = c()
          autocorrelations = c()
          ang.vars = c()
          
          for(k in include.bouts){
            
            data.temp.sub = subset(data.temp, subset = bouts == k)
            
            bout.dist = sum(data.temp.sub$trans,na.rm = TRUE)
            bout.trans.vel = mean(data.temp.sub$smooth_trans,na.rm = TRUE)
            bout.ang.vel = mean(data.temp.sub$ang,na.rm = TRUE)
            
            bout.curvature = data.temp.sub$ang/data.temp.sub$smooth_trans
            bout.curvature = mean(bout.curvature,na.rm = TRUE)
            
            bout.duration = nrow(data.temp.sub) / 32.63
            
            bout.acf = acf(as.numeric(data.temp.sub$ang),
                           plot = FALSE,
                           na.action = na.remove)[10]
            
            bout.acf = as.numeric(bout.acf$acf)
            
            bout.ang.var = var(data.temp.sub$ang,na.rm = TRUE)
            
            dists = append(dists,bout.dist)
            trans.vels = append(trans.vels,bout.trans.vel)
            ang.vels = append(ang.vels,bout.ang.vel)
            curvatures = append(curvatures,bout.curvature)
            durations = append(durations,bout.duration)
            autocorrelations = append(autocorrelations,bout.acf)
            ang.vars = append(ang.vars,bout.ang.var)
            
          }
          
          curr.bout.length = dists#mean(dists,na.rm = TRUE)
          curr.bout.mean.trans.vel = trans.vels#mean(trans.vels, na.rm = TRUE)
          curr.bout.mean.ang.vel = ang.vels#mean(ang.vels, na.rm = TRUE)
          curr.bout.curvature = curvatures#mean(curvatures,na.rm = TRUE)
          
          curr.bout.duration = durations#mean(durations, na.rm = TRUE)
          curr.bout.autocorrelation = autocorrelations#mean(autocorrelations,na.rm = TRUE)
          
          curr.bout.var = ang.vars
        }
        
        fly.data = data.frame(nBouts = curr.n.bouts,
                              bout_length = curr.bout.length,
                              bout_duration = curr.bout.duration,
                              mean_trans_vel = curr.bout.mean.trans.vel,
                              mean_ang_vel = curr.bout.mean.ang.vel,
                              curvature = curr.bout.curvature,
                              fly = j,
                              autocorrelation = curr.bout.autocorrelation,
                              ang_var = curr.bout.var)
        
        for (l in c(1:length(exp.levels))) {
          
          fly.data[,9+l] = curr.name[l]
          
        }
        
        colnames(fly.data)[c(10:(9+n.levels))] = exp.levels
        
        plot.data = rbind(plot.data,fly.data)
        
      }
      
    }
    
  }
  
  plot.data = plot.data[-1,]
  
  return(plot.data)
  
}


temporally_sparse_df = function(input.lists,exp.levels,start.frame,end.frame){
  
  n.levels = length(exp.levels)
  
  n.frames = length(c(start.frame:end.frame))
  
  time.temp = c(1:n.frames)/32.63
  
  plot.data = as.data.frame(matrix(NA, ncol = (5+n.levels)))
  colnames(plot.data) = c("fly","ang_vel","trans_vel","stimulus_moving","stimulus_visible",exp.levels)
  
  exp.names = strsplit( names(input.lists), " ")
  
  
  for(i in c(1:length(input.lists))){
    
    curr.exp = input.lists[[i]]
    
    curr.name = exp.names[[i]]
    
    curr.trans.vel = curr.exp$`Velocity`
    curr.ang.vel = curr.exp$`Angular Velocity`
    curr.stim.moving = curr.exp$`Stimulus moving`
    curr.stim.visible = curr.exp$`Stimulus visible`
    
    delay.frames =  length(curr.trans.vel[1,])- length(curr.stim.moving)
    
    if(delay.frames > 0)
    {
      curr.stim.moving = curr.stim.moving[c(delay.frames:length(curr.stim.moving))]
      curr.stim.visible = curr.stim.visible[c(delay.frames:length(curr.stim.visible))]
    } else{
      
      curr.stim.moving = c(rep(NA,abs(delay.frames)),curr.stim.moving)
      curr.stim.visible = c(rep(NA,abs(delay.frames)),curr.stim.visible)
      
    }
    
    for (j in c(1:nrow(curr.trans.vel))) {
      
      data.temp = data.frame(fly = j,
                             ang_vel = as.numeric(curr.ang.vel[j,c(start.frame:end.frame)]),
                             trans_vel = as.numeric(curr.trans.vel[j,c(start.frame:end.frame)]),
                             stimulus_moving = curr.stim.moving[c(start.frame:end.frame)],
                             stimulus_visible = curr.stim.visible[c(start.frame:end.frame)])
      
      
      for(k in  c(1:n.levels)){
        
        data.temp[,(5+k)] = curr.name[k]
        
      }
      
      colnames(data.temp)[c(6:(5+n.levels))] = exp.levels
      
      plot.data = rbind(plot.data,data.temp)  
      
    }
    
  }
  
  return(plot.data)  
  
}

timecourse_matrix_per_bout = function(input.lists,exp.levels,start.frame,end.frame,radius){
  
  n.levels = length(exp.levels)
  
  n.frames = length(c(start.frame:end.frame))
  
  time.temp = c(1:n.frames)/32.63
  
  plot.data = as.data.frame(matrix(NA, ncol = (8+n.levels)))
  colnames(plot.data) = c("x","y","trans_vel","ang_vel","dist_to_center","time","fly","bout",exp.levels)
  
  exp.names = strsplit( names(input.lists), " ")
  
  for(i in c(1:length(input.lists)))
  {
    
    curr.exp = input.lists[[i]]
    
    x.temp = curr.exp$`X mm`
    y.temp = curr.exp$`Y mm`
    
    x.temp = x.temp - curr.exp$`Arena X1`
    y.temp = y.temp - curr.exp$`Arena y1`
    
    if(dim(x.temp)[1] > dim(x.temp)[2]) x.temp = t(x.temp)
    if(dim(y.temp)[1] > dim(y.temp)[2]) y.temp = t(y.temp)
    
    trans.temp = curr.exp$`Smoothed Translational Velocity`
    ang.temp = curr.exp$`Smoothed Angular Velocity`
    
    bouts.temp = curr.exp$`Walking bouts`
    
    dist.temp = curr.exp$`Distance to center`
    
    time.temp = bouts.temp
    
    #remove later
    #error.temp = curr.exp$Flow_offset
    
    for(j in c(1:nrow(bouts.temp))){
      
      curr.fly.bouts = as.numeric(bouts.temp[j,])
      
      curr.dist = as.numeric(dist.temp[j,])
      
      exclude.bouts = as.numeric(names(table(curr.fly.bouts[which(curr.dist >= radius)])))
      
      curr.fly.bouts[which(curr.fly.bouts %in% exclude.bouts)] = NA
      
      bout.values = table(curr.fly.bouts)
      
      num.bout.values = as.numeric(names(bout.values))
      max.bout = max(num.bout.values, na.rm = TRUE)
      
      
      if(max.bout != Inf & max.bout != -Inf & !is.na(max.bout))
      {
        for(k in c(1:length(bout.values))){
          
          time.sequence = c(0:(bout.values[k]-1))/32.63
          
          curr.fly.bouts[which(curr.fly.bouts == num.bout.values[k])] = time.sequence
          
        }
        
        
        time.temp[j,] = curr.fly.bouts
      } else
      {
        time.temp[j,] = NA
      }
      
      
      
    }
    
    if(dim(trans.temp)[1] > dim(trans.temp)[2]) trans.temp = t(trans.temp)
    if(dim(ang.temp)[1] > dim(ang.temp)[2]) ang.temp = t(ang.temp)
    
    dist.temp = curr.exp$`Distance to center`
    
    if(dim(dist.temp)[1] > dim(dist.temp)[2]) dist.temp = t(dist.temp)
    
    n.flies = nrow(x.temp)
    
    experiment.name = exp.names[[i]]
    
    
    if(is.null(n.flies))
    {
      
      temp.df = data.frame(x = as.numeric(x.temp[c(start.frame:end.frame)]),
                           y = as.numeric(x.temp[c(start.frame:end.frame)]),
                           trans_vel = as.numeric(trans.temp[c(start.frame:end.frame)]),
                           ang_vel = as.numeric(ang.temp[c(start.frame:end.frame)]),
                           dist_to_center = as.numeric(dist.temp[c(start.frame:end.frame)]),
                           time = as.numeric(time.temp[j,c(start.frame:end.frame)]),
                           fly = 1,
                           bout = as.numeric(bouts.temp)[c(start.frame:end.frame)])
      
      #temp.df$fly = 1
      
      for(j in c(1:n.levels))
      {
        temp.df[,(ncol(temp.df)+1)] = experiment.name[j]
        colnames(temp.df)[ncol(temp.df)] = exp.levels[j]
      }
      
      plot.data = rbind(plot.data,temp.df)
      
    } else{
      
      for(j in c(1:n.flies))
      {
        
        temp.df = data.frame(x = as.numeric(x.temp[j,c(start.frame:end.frame)]),
                             y = as.numeric(y.temp[j,c(start.frame:end.frame)]),
                             trans_vel = as.numeric(trans.temp[j,c(start.frame:end.frame)]),
                             ang_vel = as.numeric(ang.temp[j,c(start.frame:end.frame)]),
                             dist_to_center = as.numeric(dist.temp[j,c(start.frame:end.frame)]),
                             time = as.numeric(time.temp[j,c(start.frame:end.frame)]),
                             fly = j,
                             bout = as.numeric(bouts.temp[j,c(start.frame:end.frame)]))
        
        
        #temp.df$fly = j
        
        for(k in c(1:n.levels))
        {
          temp.df[,(ncol(temp.df)+1)] = experiment.name[k]
          colnames(temp.df)[ncol(temp.df)] = exp.levels[k]
        }
        
        
        plot.data = rbind(plot.data,temp.df)
        
      }
      
    }
    
    
  }
  
  plot.data = plot.data[-1,]
  
  return(plot.data)
  
}

colored_path_plot = function(experiment,fly,exp.levels,start.frame,end.frame){
  
  library(ggplot2)
  library(scico)
  
  data = timecourse.data.from.list(experiment,fly,exp.levels)
  data = data[c(start.frame:end.frame),]
  
  #normalize vels
  max.trans = max(data$translational_velocity,na.rm = TRUE)
  max.ang = max(abs(data$angular_velocity),na.rm = TRUE)
  
  data$translational_velocity = data$translational_velocity/max.trans
  data$angular_velocity = data$angular_velocity/max.ang
  
  output = list()
  
  output[[1]] = paste("Maximum translational velocity: ",max.trans," mm/s",sep = "")
  output[[2]] = paste("Maximum angular veloctity: ",max.ang," rad./s",sep = "")
  
  output[[3]] = ggplot(data,aes(x=x,y=y,color=angular_velocity))+
    geom_path(size=1)+
    scale_color_scico(palette = "romaO",limits = c(-1,1))+
    coord_fixed(ratio = 1)+
    #ylim(c(0,ymax))+
    #xlim(c(0,xmax))+
    theme_minimal()+
    theme(panel.grid = element_blank(),
          axis.line.y = element_line(size=1,linetype = "solid"),
          axis.line.x = element_line(size=1,linetype = "solid"),
          axis.ticks = element_line(size=1,linetype = "solid"),
          legend.title = element_blank())
  
  output[[4]] = ggplot(data,aes(x=x,y=y,color=translational_velocity))+
    geom_path(size=1)+
    scale_color_scico(palette = "romaO",begin=0.5,end=1)+
    coord_fixed(ratio = 1)+ 
    #ylim(0,ymax)+
    #xlim(0,xmax)+
    theme_minimal()+
    theme(panel.grid = element_blank(),
          axis.line.y = element_line(size=1,linetype = "solid"),
          axis.line.x = element_line(size=1,linetype = "solid"),
          axis.ticks = element_line(size=1,linetype = "solid"),
          legend.title = element_blank(),
          legend.position = "none")
  
  data.lineplot = rbind(data,data)
  data.lineplot$Velocity = c(data$translational_velocity,data$angular_velocity)
  data.lineplot$Velocity_type = c(rep("Translational",nrow(data)),rep("Angular",nrow(data)))
  data.lineplot$Time = c(c(1:nrow(data))/30,c(1:nrow(data))/30)
  
  output[[5]] = ggplot(data.lineplot,aes(x=Time,y=Velocity,color=Velocity_type))+
    geom_hline(yintercept = 0,linetype = "dashed",alpha =1,col= "black",size=1)+
    geom_line(size=1)+
    theme_minimal()+
    ylim(c(-1,1))+
    scale_color_manual(values=c("#f194a3","#45ada8"))+
    theme(axis.line.y = element_line(size=1,linetype = "solid"),
          axis.text.y = element_text(face="bold",size=15),
          axis.text.x = element_text(face="bold",size=15),
          legend.position = "none")
  
  
  return(output)
  
}


flow_offset_per_xy = function(fly_x,fly_y,fly_ori,center_x,center_y,clockwise){
  
  centered_fly_x = fly_x - center_x
  centered_fly_y = center_y - fly_y
  if(!clockwise) centered_fly_y = fly_y - center_y
  
  angular_position = atan2(centered_fly_y, centered_fly_x)
  
  flow_offset = (angular_position - fly_ori + pi/2) %% pi - pi/2
  
  return(flow_offset)
}

add_error_to_list = function(exp_list){
  
  for(i in c(1:length(exp_list))){
    
    curr_exp = exp_list[[i]]
    curr_name = names(exp_list)[i]
    
    x = curr_exp$X[,1]
    y = curr_exp$Y[,1]
    ori = curr_exp$Ori[1,]
    
    center = curr_exp$Center
    
    clockwise = length(grep("Clockwise",curr_name)) > 0
    
    flow_offset = flow_offset_per_xy(x,y,ori,center[1],center[2],clockwise)
    
    curr_exp$Flow_offset = flow_offset
    exp_list[[i]] = curr_exp
  }
  
  return(exp_list)
}

rotate_around_fulcrum = function(original_x,original_y,fulcrum_x,fulcrum_y,rotation_angle){
  
  centered_x = original_x - fulcrum_x
  centered_y = original_y - fulcrum_y
  
  dist = sqrt(centered_x*centered_x+centered_y*centered_y)
  
  original_angle = atan2(centered_y,centered_x)
  
  #new_angle = (original_angle - rotation_angle + pi/2) %% pi - pi/2
  new_angle = original_angle - rotation_angle
  
  
  new_x = dist*cos(new_angle)
  new_y = dist*sin(new_angle)
  
  
  output_df = data.frame(rotated_x = new_x,rotated_y=new_y)
  
}

opto.single.trial.path.plot.shifted.heading = function(experiment,fly,delivery){
  
  on = experiment$`On frames`$light.on[delivery]
  off = experiment$`Off frames`$light.off[delivery]
  
  interval.length =  off - on
  
  start = on - round(5*32.63)
  end = off #+ interval.length
  
  on = on - start
  off = off - start
  
  x = as.numeric(experiment$`X mm`[fly,c(start:end)])
  y = as.numeric(experiment$`Y mm`[fly,c(start:end)])
  ori = as.numeric(experiment$Ori[c(start:end),fly])
  
  phase = rep("pre", length(x))
  phase[c(c(on:length(phase)))] = "on"
  #phase[c(1:(on-1))] = "pre"
  
  data = data.frame(x= x, y= y,phase = phase)
  rect = data.frame(on = on/32.63, off = off/32.63)
  
  x_initial = x[(on-1)]
  y_initial = y[(on-1)]
  ori_initial = ori[(on-1)]
  
  rotated_data = rotate_around_fulcrum(x,y,x_initial,y_initial,ori_initial)
  
  data = cbind(data,rotated_data)
  
  #plots = list()
  
  plot = ggplot(data,aes(x= x,y=y, color = phase))+
    geom_path()+
    theme_minimal()+
    coord_fixed()+
    scale_color_manual(values = c("#f194a3","black"))+
    theme(panel.grid = element_blank(),
          axis.line.y = element_line(size=1,linetype = "solid"),
          axis.text.y = element_text(face="bold",size=15),
          axis.text.x = element_text(face="bold",size=15))
  
  
  return(plot)
  
  
}
