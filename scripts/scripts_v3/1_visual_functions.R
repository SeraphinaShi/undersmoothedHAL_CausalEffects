##############################################################
# This function runs the simulation for B rounds with given:
#       - data generating function
#       - sample size n
#       - number of simulations: B

# returns the estimated ATE, empirical 95% CI, and coverage rates
##############################################################
plot_perforences_1lambda_alla <- function(df, z_para=1, est_plot_only=F, plot_list=F, add_oracal=F){
  
  df <- df %>% filter(z==z_para)
  
  color = ifelse(z_para==1, "#F8766D", "#00BFC4")
  oracal_color = ifelse(z_para==1,"#00BA38", "#F564E3")

  p_est_avg <- ggplot(data=df, aes(x=a)) +
    geom_line(aes(y=psi0), alpha = 0.5, color="darkgrey") +
    geom_errorbar(aes(ymin=ci_lwr, ymax=ci_upr, color='Empirical'), width=0.7) +
    geom_point(aes(y=psi0), color = "black") +
    geom_point(aes(y=psi_hat, color='Empirical'), shape=17, size=2) +
    labs(x="a", y="ATE", title = "Estimation") +
    scale_color_manual(name='Method',
      breaks=c('Empirical', 'Oracal'),
      values=c('Empirical'=color, 'Oracal'=oracal_color)) +
    scale_x_continuous(limits = c(0, 5), breaks = 0:5)

  if(add_oracal){
    p_est_avg <- p_est_avg + 
      geom_errorbar(aes(ymin=oracal_ci_lwr, ymax=oracal_ci_upr, color='Oracal'), width=0.7) 
  } 
  
  legend <- get_legend(p_est_avg)
  p_est_avg <- p_est_avg + theme(legend.position='none')
  
  if(est_plot_only){
    return(p_est_avg + labs(title=""))
  } else {
    p_bias <- ggplot(df, aes(x = a, y = bias)) +  
      geom_line(color=color) +
      geom_point(color=color) + 
      labs(title="|Bias|") 
    
    if(!add_oracal){
      p_se <- ggplot(df, aes(x = a, y = SE)) +  
        geom_line(color=color) +
        geom_point(color=color) + 
        labs(title="Standard Error") 
      
      p_bias_d_df <- ggplot(df, aes(x = a, y = bias_se_ratio)) +  
        geom_line(color=color) + 
        geom_point(color=color) + 
        labs(title="|Bias| / Standard Error") 
      
      p_cr <- ggplot(df, aes(x = a, y = cover_rate)) +  
        geom_line(color=color) + 
        geom_point(color=color) +
        labs(title="95% CI Coverage Rate")
    } else {
      p_se <- ggplot(df, aes(x = a)) +  
        geom_line(aes(y = SE, color='Empirical'),alpha=0.7) +
        geom_point(aes(y = SE, color='Empirical'),alpha=0.7) + 
        geom_line(aes(y = oracal_SE, color='Oracal'),alpha=0.7) +
        geom_point(aes(y = oracal_SE, color='Oracal'),alpha=0.7) + 
        scale_color_manual(name='Method',
          breaks=c('Empirical', 'Oracal'),
          values=c('Empirical'=color, 'Oracal'=oracal_color)) + 
        theme(legend.position='none') +
        labs(title="Standard Error") 
      
      
      p_bias_d_df <- ggplot(df, aes(x = a)) +  
        geom_line(aes(y = bias_se_ratio, color='Empirical'), alpha=0.7) +
        geom_point(aes(y = bias_se_ratio, color='Empirical'), alpha=0.7) + 
        geom_line(aes(y = oracal_bias_se_ratio, color='Oracal'), alpha=0.7) +
        geom_point(aes(y = oracal_bias_se_ratio, color='Oracal'), alpha=0.7) + 
        scale_color_manual(name='Method',
          breaks=c('Empirical', 'Oracal'),
          values=c('Empirical'=color, 'Oracal'=oracal_color)) + 
        theme(legend.position='none') +
        labs(title="|Bias| / Standard Error") 
      
      p_cr <- ggplot(df, aes(x = a)) +  
        geom_line(aes(y = cover_rate, color='Empirical'), alpha=0.7) +
        geom_point(aes(y = cover_rate, color='Empirical'), alpha=0.7) + 
        geom_line(aes(y = oracal_cover_rate, color='Oracal'), alpha=0.7) +
        geom_point(aes(y = oracal_cover_rate, color='Oracal'), alpha=0.7) + 
        scale_color_manual(name='Method',
          breaks=c('Empirical', 'Oracal'),
          values=c('Empirical'=color, 'Oracal'=oracal_color)) + 
        theme(legend.position='none') +
        labs(title="95% CI Coverage Rate")
    }
    
    
    if(plot_list){
      return(list(p_est_avg, p_bias, p_se, p_bias_d_df, p_cr))
    } else {
      p <- grid.arrange(p_est_avg, p_bias, p_se, p_bias_d_df, p_cr, legend,
                          layout_matrix = rbind(c(NA,1,1,NA),
                                                c(NA,1,1,6),
                                                c(2,2,3,3),
                                                c(2,2,3,3),
                                                c(4,4,5,5),
                                                c(4,4,5,5)),
                          top = textGrob(paste0("HAL-based plug in estimator performence for ATE"), 
                                         gp=gpar(fontsize=11, fontface = 'bold')))
        return(p)
    }
  }

}


estimation_qqplot <- function(results_list, z_para){
  df <- do.call("rbind", results_list) %>% 
    as.data.frame() %>% filter(z == z_para)
  
  p <- ggplot(df, aes(sample = psi_hat)) + 
    stat_qq() + stat_qq_line() +
    # theme(axis.title.x=element_blank(), 
    #       axis.title.y=element_blank()) +
    labs(x = "Theoretical Quantiles",
         y = "Sample Quantiles",
         title = paste0("Normal Q-Q Plot for Estimated ATE when z=",z_para)) +
    facet_wrap(~a, nrow=2) 
  
  return(p)
}



plot_perforences_alllambda_1a <- function(df, a_para, z_para, add_oracal=F, u_scaler=NA, avg_u_scaler=NA){
  df <- df %>% filter(a == a_para, z == z_para)
  
  p_est_avg <- ggplot(df) +  
    geom_line(aes(x = lambda_scaler, y = psi_hat), color = "grey") + 
    geom_point(aes(x = lambda_scaler, y = psi_hat)) + 
    geom_hline(aes(yintercept=psi0)) +
    labs(title="Estimation average") +
    theme()
  
  p_bias <- ggplot(df, aes(x = lambda_scaler, y = bias)) +  
    geom_line(color = "grey") + 
    geom_point() + 
    labs(title="Bias") 
  
  p_se <- ggplot(df, aes(x = lambda_scaler)) +  
    geom_line(aes(y = SE), color = "grey") + 
    geom_point(aes(y = SE, color = "Empirical")) + 
    labs(title="Standard Error") + 
    scale_color_manual(name='method',
                       breaks=c('Empirical', 'Oracal'),
                       values=c('Empirical'='black', 'Oracal'='#8B6508'))

  
  p_bias_d_df <- ggplot(df, aes(x = lambda_scaler)) + 
    geom_line(aes(y = bias_se_ratio), color = "grey") + 
    geom_point(aes(y = bias_se_ratio, color = "Empirical")) + 
    labs(title="|Bias| / Standard Error") + 
    scale_color_manual(name='method',
                       breaks=c('Empirical', 'Oracal'),
                       values=c('Empirical'='black', 'Oracal'='#8B6508')) + 
    theme(legend.position='none') 
  
  p_cr <- ggplot(df, aes(x = lambda_scaler)) +  
    geom_line(aes(y = cover_rate), color = "grey") + 
    geom_point(aes(y = cover_rate, color = "Empirical")) + 
    labs(title="Coverage rate") + 
    scale_color_manual(name='method',
                       breaks=c('Empirical', 'Oracal'),
                       values=c('Empirical'='black', 'Oracal'='#8B6508')) + 
    theme(legend.position='none')  
  
  if(!add_oracal){
    legend <- get_legend(p_se)
    p_se <- p_se + theme(legend.position='none')
  } else {
    p_se <- p_se + 
      geom_line(aes(y = oracal_SE), color = "#EEC591") + 
      geom_point(aes(y = oracal_SE, color = "Oracal")) 
    
    legend <- get_legend(p_se)
    p_se <- p_se + theme(legend.position='none')
    
    p_bias_d_df <- p_bias_d_df +  
      geom_line(aes(y = oracal_bias_se_ratio), color = "#EEC591") + 
      geom_point(aes(y = oracal_bias_se_ratio, color = "Oracal")) 
    
    p_cr <- p_cr +  
      geom_line(aes(y = oracal_cover_rate), color = "#EEC591") + 
      geom_point(aes(y = oracal_cover_rate, color = "Oracal")) 
  }
  
  if(!is.na(avg_u_scaler)){
    p_est_avg <- p_est_avg +
      geom_vline(xintercept = avg_u_scaler, lty=2)
    p_bias <- p_bias +
      geom_vline(xintercept = avg_u_scaler, lty=2)
    p_se <- p_se + 
      geom_vline(xintercept = avg_u_scaler, lty=2)
    p_bias_d_df <- p_bias_d_df +  
      geom_vline(xintercept = avg_u_scaler, lty=2)
    p_cr <- p_cr +
      geom_vline(xintercept = avg_u_scaler, lty=2)
  }
  
  p <- grid.arrange(p_est_avg, p_bias, p_se, p_bias_d_df, p_cr, legend, 
                    layout_matrix = rbind(c(NA,1,1,NA),
                                          c(NA,1,1,6),
                                          c(2,2,3,3),
                                          c(2,2,3,3),
                                          c(4,4,5,5),
                                          c(4,4,5,5)),
                    top = textGrob(paste0("HAL-based plug in estimator performence for a=", a_para, ", z=", z_para), 
                                   gp=gpar(fontsize=11, fontface = 'bold')))
  return(p)
}

