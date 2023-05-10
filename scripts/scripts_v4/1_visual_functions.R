
plot_perforences_cv_u_gl_alla <- function(df, save_plot=NA){
  
  color_cv =  "#F8766D"
  color_u_g = "#00BA38"
  color_u_l = '#619CFF'
  
  p_lambda <- ggplot(data=df, aes(x=a)) +
    geom_line(aes(y=lambda_scaler, color=method), alpha = 0.7) +
    geom_point(aes(y=lambda_scaler, color=method), shape=17, size=2, alpha= 0.7) +
    labs(x="a", y="lambda_scaler", 
         title = "Lambda scaler",
         subtitle = paste0("upon CV_lambda ", round(mean(df$lambda[df$method == 'CV']), 6)))+
    scale_x_continuous(limits = c(0, 5), breaks = 0:5) +
    scale_color_manual(name='Method',
                       breaks=c('CV', 'U_G', 'U_L'),
                       values=c('CV'=color_cv, 'U_G'=color_u_g, 'U_L'=color_u_l)) +
    scale_linetype_manual(breaks=c('Empirical', 'Oracal'),
                          values=c('Empirical'= 1, 'Oracal'=2)) +
    theme_bw()+
    theme(legend.position='none') 
  
  p_est_avg <- ggplot(data=df, aes(x=a)) +
    geom_line(aes(y=psi0), alpha = 0.5, color="darkgrey") +
    geom_errorbar(aes(ymin=ci_lwr, ymax=ci_upr, color=method, linetype='Empirical'), width=0.7) +
    geom_errorbar(aes(ymin=oracal_ci_lwr, ymax=oracal_ci_upr, color=method, linetype = "Oracal"), width=0.7) +
    geom_point(aes(y=psi0), color = "black") +
    geom_point(aes(y=y_hat, color=method), shape=17, size=2, alpha= 0.7) +
    labs(x="a", y="E[Y|a, W]", title = "Estimation") +
    scale_x_continuous(limits = c(0, 5), breaks = 0:5) +
    scale_color_manual(name='Method',
                       breaks=c('CV', 'U_G', 'U_L'),
                       values=c('CV'=color_cv, 'U_G'=color_u_g, 'U_L'=color_u_l)) +
    scale_linetype_manual(breaks=c('Empirical', 'Oracal'),
                          values=c('Empirical'= 1, 'Oracal'=2)) +
    theme_bw() +
    theme(legend.box = "horizontal")
  
  legend <- get_legend(p_est_avg)
  p_est_avg <- p_est_avg + theme(legend.position='none')
  
  p_bias <- ggplot(df, aes(x = a, y = bias)) +  
    geom_line(aes(color=method)) +
    geom_point(aes(color=method)) + 
    labs(title="|Bias|") +
    scale_x_continuous(limits = c(0, 5), breaks = 0:5) +
    scale_color_manual(name='Method',
                       breaks=c('CV', 'U_G', 'U_L'),
                       values=c('CV'=color_cv, 'U_G'=color_u_g, 'U_L'=color_u_l)) +
    theme_bw() +
    theme(legend.position='none') 
  
  p_se <- ggplot(df, aes(x = a)) +  
    geom_line(aes(y = SE, color=method, linetype='Empirical'),alpha=0.7) +
    geom_point(aes(y = SE, color=method, linetype='Empirical'),alpha=0.7) + 
    geom_line(aes(y = oracal_SE, color=method, linetype='Oracal'),alpha=0.7) +
    geom_point(aes(y = oracal_SE, color=method, linetype='Oracal'),alpha=0.7) + 
    labs(title="Standard Error") +
    scale_x_continuous(limits = c(0, 5), breaks = 0:5) +
    scale_color_manual(name='Method',
                       breaks=c('CV', 'U_G', 'U_L'),
                       values=c('CV'=color_cv, 'U_G'=color_u_g, 'U_L'=color_u_l)) +
    scale_linetype_manual(breaks=c('Empirical', 'Oracal'),
                          values=c('Empirical'= 1, 'Oracal'=2))+
    theme_bw() +
    theme(legend.position='none') 
  
  
  p_bias_d_df <- ggplot(df, aes(x = a)) +  
    geom_line(aes(y = bias_se_ratio, color=method, linetype='Empirical'), alpha=0.7) +
    geom_point(aes(y = bias_se_ratio, color=method, linetype='Empirical'), alpha=0.7) + 
    geom_line(aes(y = oracal_bias_se_ratio, color=method, linetype='Oracal'), alpha=0.7) +
    geom_point(aes(y = oracal_bias_se_ratio, color=method, linetype='Oracal'), alpha=0.7) + 
    labs(title="|Bias| / Standard Error") +
    scale_x_continuous(limits = c(0, 5), breaks = 0:5) +
    scale_color_manual(name='Method',
                       breaks=c('CV', 'U_G', 'U_L'),
                       values=c('CV'=color_cv, 'U_G'=color_u_g, 'U_L'=color_u_l)) +
    scale_linetype_manual(breaks=c('Empirical', 'Oracal'),
                          values=c('Empirical'= 1, 'Oracal'=2))+
    theme_bw() +
    theme(legend.position='none') 
  
  p_cr <- ggplot(df, aes(x = a)) +  
    geom_line(aes(y = cover_rate, color=method, linetype='Empirical'), alpha=0.7) +
    geom_point(aes(y = cover_rate, color=method, linetype='Empirical'), alpha=0.7) + 
    geom_line(aes(y = oracal_cover_rate, color=method, linetype='Oracal'), alpha=0.7) +
    geom_point(aes(y = oracal_cover_rate, color=method, linetype='Oracal'), alpha=0.7) + 
    labs(title="95% CI Coverage Rate")+
    scale_x_continuous(limits = c(0, 5), breaks = 0:5) +
    scale_color_manual(name='Method',
                       breaks=c('CV', 'U_G', 'U_L'),
                       values=c('CV'=color_cv, 'U_G'=color_u_g, 'U_L'=color_u_l)) +
    scale_linetype_manual(breaks=c('Empirical', 'Oracal'),
                          values=c('Empirical'= 1, 'Oracal'=2))+
    theme_bw() +
    theme(legend.position='none') 
  
  
  p <- grid.arrange(p_est_avg, p_bias, p_se, p_bias_d_df, p_cr, legend, p_lambda,
                    layout_matrix = rbind(c(1,1,1,NA,6,6),
                                          c(1,1,1,7,7,7),
                                          c(1,1,1,7,7,7),
                                          c(2,2,2,3,3,3),
                                          c(2,2,2,3,3,3),
                                          c(2,2,2,3,3,3),
                                          c(4,4,4,5,5,5),
                                          c(4,4,4,5,5,5),
                                          c(4,4,4,5,5,5)),
                    top = textGrob(paste0("HAL-based plug-in estimator performence for E[Y|a,W]"), 
                                   gp=gpar(fontsize=15, fontface = 'bold')))
  
  if(!is.na(save_plot)){
    ggsave(save_plot, plot=p, width = 8, height = 8)
  }
  return(p)
  
}


estimation_qqplot_cv_u_gl_alla <- function(results_list, save_plot=NA){
  df <- data.frame()
  for (method in c("CV", "U_G", "U_L")){
    result_all <-  do.call("rbind", results_list[[method]]$all_results) %>% as.data.frame()
    result_all$method = method
    
    df <- rbind(df, result_all)
  }
  
  
  p <- ggplot(df, aes(sample = y_hat)) + 
    stat_qq() + stat_qq_line() +
    labs(x = "Theoretical Quantiles",
         y = "Sample Quantiles",
         title = paste0("Q-Q Plot for Estimated E[Y|a,W]")) +
    facet_grid(method ~ a)
  
  if(!is.na(save_plot)){
    ggsave(save_plot, plot=p, width = 15, height = 5, dpi = 1200)
  }
  
  return(p)
}



plot_perforences_alllambda_1a <- function(df, a_para, z_para, add_oracal=F, u_g_scaler=NA, u_l_scaler=NA){
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
  
  if(!is.na(u_g_scaler)){
    p_est_avg <- p_est_avg +
      geom_vline(xintercept = u_g_scaler, lty=2, col = "#00BA38") +
      # geom_text(aes(x=u_g_scaler, label="Globally undersmoothed HAL\n", y=+Inf), colour="#00BA38", angle=90, text=element_text(size=5)) +
      geom_vline(xintercept = 1, lty=2, col = "#F8766D")
    p_bias <- p_bias +
      geom_vline(xintercept = u_g_scaler, lty=2, col = "#00BA38")+
      geom_vline(xintercept = 1, lty=2, col = "#F8766D")
    p_se <- p_se + 
      geom_vline(xintercept = u_g_scaler, lty=2, col = "#00BA38")+
      geom_vline(xintercept = 1, lty=2, col = "#F8766D")
    p_bias_d_df <- p_bias_d_df +  
      geom_vline(xintercept = u_g_scaler, lty=2, col = "#00BA38")+
      geom_vline(xintercept = 1, lty=2, col = "#F8766D")
    p_cr <- p_cr +
      geom_vline(xintercept = u_g_scaler, lty=2, col = "#00BA38")+
      geom_vline(xintercept = 1, lty=2, col = "#F8766D")
  }
  
  if(!is.na(u_l_scaler)){
    p_est_avg <- p_est_avg +
      geom_vline(xintercept = u_l_scaler, lty=2, col = "#619CFF") 
    p_bias <- p_bias +
      geom_vline(xintercept = u_l_scaler, lty=2, col = "#619CFF")
    p_se <- p_se + 
      geom_vline(xintercept = u_l_scaler, lty=2, col = "#619CFF")
    p_bias_d_df <- p_bias_d_df +  
      geom_vline(xintercept = u_l_scaler, lty=2, col = "#619CFF")
    p_cr <- p_cr +
      geom_vline(xintercept = u_l_scaler, lty=2, col = "#619CFF")
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

