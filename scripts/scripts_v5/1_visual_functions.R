
plot_performences_cv_u_gl_alla <- function(df, save_plot=NA){
  
  color_cv =  "#F8766D"
  color_u_g = "#00BA38"
  color_u_l = '#619CFF'
  
  a_max <- max(df$a)
  
  p_lambda <- ggplot(data=df, aes(x=a)) +
    geom_line(aes(y=lambda_scaler, color=method), alpha = 0.7) +
    geom_point(aes(y=lambda_scaler, color=method), shape=17, size=2, alpha= 0.7) +
    ylim(0, 1.2) +
    labs(x="a", y="lambda_scaler", 
         title = "Lambda scaler",
         subtitle = paste0("upon CV_lambda ", round(mean(df$lambda[df$method == 'CV']), 6)))+
    scale_x_continuous(limits = c(0, a_max), breaks = 0:a_max) +
    scale_color_manual(name='Method',
                       breaks=c('CV', 'U_G', 'U_L'),
                       values=c('CV'=color_cv, 'U_G'=color_u_g, 'U_L'=color_u_l)) +
    scale_linetype_manual(breaks=c('Oracal', 'Empirical', 'Bootstrap'),
                          values=c('Oracal'=1, 'Empirical'=5, 'Bootstrap'=3)) +
    theme_bw()+
    theme(legend.position='none') 
  
  p_est_avg <- ggplot(data=df, aes(x=a)) +
    geom_line(aes(y=psi0), alpha = 0.5, color="darkgrey") +
    geom_ribbon(aes(ymin=ci_lwr, ymax=ci_upr, color=method, fill=method, linetype='Empirical'), width=0.7, alpha=0.1) +
    geom_ribbon(aes(ymin=ci_lwr_bt, ymax=ci_upr_bt, color=method, fill=method, linetype = "Bootstrap"),  width=0.7, alpha=0.1) +
    geom_ribbon(aes(ymin=oracal_ci_lwr, ymax=oracal_ci_upr, color=method, fill=method, linetype = "Oracal"),  width=0.7, alpha=0.1) +
    geom_point(aes(y=psi0), color = "black") +
    geom_point(aes(y=y_hat, color=method), shape=17, size=2, alpha= 0.7) +
    labs(x="a", y="E[Y|a, W]", title = "Estimation") +
    scale_x_continuous(limits = c(0, a_max), breaks = 0:a_max) +
    scale_color_manual(name='Method',
                       breaks=c('CV', 'U_G', 'U_L'),
                       values=c('CV'=color_cv, 'U_G'=color_u_g, 'U_L'=color_u_l)) +
    scale_linetype_manual(breaks=c('Oracal', 'Empirical', 'Bootstrap'),
                          values=c('Oracal'=1, 'Empirical'=5, 'Bootstrap'=3)) +
    theme_bw() +
    theme(legend.box = "horizontal",
          legend.position='none')
  
  p_bias <- ggplot(df, aes(x = a, y = bias)) +  
    geom_line(aes(color=method)) +
    geom_point(aes(color=method)) + 
    labs(title="|Bias|") +
    scale_x_continuous(limits = c(0, a_max), breaks = 0:a_max) +
    scale_color_manual(name='Method',
                       breaks=c('CV', 'U_G', 'U_L'),
                       values=c('CV'=color_cv, 'U_G'=color_u_g, 'U_L'=color_u_l)) +
    theme_bw() +
    theme(legend.position='none') 
  
  p_se <- ggplot(df, aes(x = a)) +  
    geom_line(aes(y = SE, color=method, linetype='Empirical'),alpha=0.7) +
    geom_point(aes(y = SE, color=method, linetype='Empirical'),alpha=0.7) + 
    geom_line(aes(y = SE_bt, color=method, linetype='Bootstrap'),alpha=0.7) +
    geom_point(aes(y = SE_bt, color=method, linetype='Bootstrap'),alpha=0.7) + 
    geom_line(aes(y = oracal_SE, color=method, linetype='Oracal'),alpha=0.7) +
    geom_point(aes(y = oracal_SE, color=method, linetype='Oracal'),alpha=0.7) + 
    labs(title="Standard Error") +
    scale_x_continuous(limits = c(0, a_max), breaks = 0:a_max) +
    scale_color_manual(name='Method',
                       breaks=c('CV', 'U_G', 'U_L'),
                       values=c('CV'=color_cv, 'U_G'=color_u_g, 'U_L'=color_u_l)) +
    scale_linetype_manual(breaks=c('Oracal', 'Empirical', 'Bootstrap'),
                          values=c('Oracal'=1, 'Empirical'=5, 'Bootstrap'=3)) +
    theme_bw() + 
    theme(legend.box = "horizontal")
  
  legend <- get_legend(p_se)
  p_se <- p_se + theme(legend.position='none')
  
  
  p_bias_d_df <- ggplot(df, aes(x = a)) +  
    geom_line(aes(y = bias_se_ratio, color=method, linetype='Empirical'), alpha=0.7) +
    geom_point(aes(y = bias_se_ratio, color=method, linetype='Empirical'), alpha=0.7) + 
    geom_line(aes(y = bias_se_ratio_bt, color=method, linetype='Bootstrap'), alpha=0.7) +
    geom_point(aes(y = bias_se_ratio_bt, color=method, linetype='Bootstrap'), alpha=0.7) + 
    geom_line(aes(y = oracal_bias_se_ratio, color=method, linetype='Oracal'), alpha=0.7) +
    geom_point(aes(y = oracal_bias_se_ratio, color=method, linetype='Oracal'), alpha=0.7) + 
    labs(title="|Bias| / Standard Error") +
    scale_x_continuous(limits = c(0, a_max), breaks = 0:a_max) +
    scale_color_manual(name='Method',
                       breaks=c('CV', 'U_G', 'U_L'),
                       values=c('CV'=color_cv, 'U_G'=color_u_g, 'U_L'=color_u_l)) +
    scale_linetype_manual(breaks=c('Oracal', 'Empirical', 'Bootstrap'),
                          values=c('Oracal'=1, 'Empirical'=5, 'Bootstrap'=3)) +
    theme_bw() +
    theme(legend.position='none') 
  
  p_cr <- ggplot(df, aes(x = a)) +  
    geom_rect(data=NULL,aes(xmin=-Inf,xmax=Inf,ymin=0.95,ymax=Inf), fill="khaki1", alpha = 0.1)+ # fill="darkseagreen1"
    geom_line(aes(y = cover_rate, color=method, linetype='Empirical'), alpha=0.7) +
    geom_point(aes(y = cover_rate, color=method, linetype='Empirical'), alpha=0.7) + 
    geom_line(aes(y = cover_rate_bt, color=method, linetype='Bootstrap'), alpha=0.7) +
    geom_point(aes(y = cover_rate_bt, color=method, linetype='Bootstrap'), alpha=0.7) + 
    geom_line(aes(y = oracal_cover_rate, color=method, linetype='Oracal'), alpha=0.7) +
    geom_point(aes(y = oracal_cover_rate, color=method, linetype='Oracal'), alpha=0.7) + 
    labs(title="95% CI Coverage Rate")+
    scale_x_continuous(limits = c(0, a_max), breaks = 0:a_max) +
    scale_color_manual(name='Method',
                       breaks=c('CV', 'U_G', 'U_L'),
                       values=c('CV'=color_cv, 'U_G'=color_u_g, 'U_L'=color_u_l)) +
    scale_linetype_manual(breaks=c('Oracal', 'Empirical', 'Bootstrap'),
                          values=c('Oracal'=1, 'Empirical'=5, 'Bootstrap'=3)) +
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
                    top = textGrob(paste0("HAL-based plug-in estimator performences for E[Y|a,W] \n"), 
                                   gp=gpar(fontsize=17)))
  
  if(!is.na(save_plot)){
    ggsave(save_plot, plot=p, width = 8, height = 9, dpi = 800)
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
  
  df <- df[!is.na(df$a), ]
  p <- ggplot(df, aes(sample = y_hat)) + 
    stat_qq() + stat_qq_line() +
    labs(x = "Theoretical Quantiles",
         y = "Sample Quantiles",
         title = paste0("Q-Q Plot for Estimated E[Y|a,W]")) +
    facet_grid(method ~ a, scales = 'free') +
    theme(plot.title = element_text(hjust = 0.5, size = 18),
          axis.text = element_text(size=10),
          axis.title = element_text(size=14),
          strip.text = element_text(size = 12))
  
  if(!is.na(save_plot)){
    ggsave(save_plot, plot=p, width = 13, height = 4, dpi = 800)
  }
  
  return(p)
}


plot_perforences_alllambda <- function(df, u_g_scaler=NA, u_l_scalers=NA, save_plot=NA, max_bias_sd=NA){
  
  legend_undersmoothing = ggplot(df) +  
    geom_line(aes(x = lambda_scaler, y = y_hat, color = "Global")) + 
    geom_line(aes(x = lambda_scaler, y = y_hat, color = "Local")) + 
    geom_line(aes(x = lambda_scaler, y = y_hat, color = "None")) + 
    theme_bw() +
    scale_colour_manual(name="Undersmoothing",
                        values=c(Global="#00BA38", Local="#619CFF", None="#F8766D"))
  legend_undersmoothing <- get_legend(legend_undersmoothing)
  
  p_est_avg_list = list()
  p_bias_list = list()
  p_se_list = list()
  p_bias_se_list = list()
  p_cr_list = list()
  legend = NA
  
  for (i in 1:length(eval_points)) {
    df_a <- df %>% filter(a == eval_points[i])
    u_l_scaler = u_l_scalers[i]
    
    p_est_avg = ggplot(df_a) +  
      geom_ribbon(aes(x = lambda_scaler, ymin=ci_lwr, ymax=ci_upr, color='Empirical', fill = 'Empirical'), width=0.7, alpha=0.5) +
      geom_ribbon(aes(x = lambda_scaler, ymin=ci_lwr_bt, ymax=ci_upr_bt,  color='Bootstrap', fill = 'Bootstrap'), width=0.7, alpha=0.5) +
      geom_ribbon(aes(x = lambda_scaler, ymin=oracal_ci_lwr, ymax=oracal_ci_upr,  color='Oracal', fill = 'Oracal'), width=0.7, alpha=0.5) +
      geom_line(aes(x = lambda_scaler, y = y_hat), color = "grey") + 
      geom_point(aes(x = lambda_scaler, y = y_hat)) + 
      geom_hline(aes(yintercept=psi0)) + 
      scale_x_continuous(breaks=seq(0,1.2,by=0.25)) +
      scale_color_manual(breaks=c('Oracal', 'Empirical', 'Bootstrap'),
                         values=c('Oracal'='darkolivegreen3', 'Empirical'='lightsalmon', 'Bootstrap'='darkorchid')) +
      scale_fill_manual(breaks=c('Oracal', 'Empirical', 'Bootstrap'),
                         values=c('Oracal'='darkolivegreen3', 'Empirical'='lightsalmon', 'Bootstrap'='darkorchid')) +
      theme_bw() +
      labs(x = "", title = paste0('a = ', eval_points[i])) +
      theme(axis.title=element_blank(),
            plot.title = element_text(hjust = 0.5),
            legend.position='none')


    p_bias <- ggplot(df_a, aes(x = lambda_scaler, y = bias)) +  
      geom_line(color = "grey") + 
      geom_point() + 
      scale_x_continuous(breaks=seq(0,1.2,by=0.25)) +
      theme_bw()+
      theme(axis.title=element_blank())
    
    p_se <- ggplot(df_a, aes(x = lambda_scaler)) +  
      geom_line(aes(y = SE, color='Empirical')) + 
      geom_point(aes(y = SE, color='Empirical')) + 
      geom_line(aes(y = SE_bt, color='Bootstrap')) + 
      geom_point(aes(y = SE_bt, color='Bootstrap')) +
      geom_line(aes(y = oracal_SE, color='Oracal')) + 
      geom_point(aes(y = oracal_SE, color='Oracal')) +
      scale_x_continuous(breaks=seq(0,1.2,by=0.25)) +
      scale_color_manual(name='method',
                         breaks=c('Oracal', 'Empirical', 'Bootstrap'),
                         values=c('Oracal'='darkolivegreen3', 'Empirical'='lightsalmon', 'Bootstrap'='darkorchid')) +
      theme_bw()+
      theme(axis.title=element_blank())
    
    legend <- get_legend(p_se)
    p_se <- p_se + theme(legend.position='none')
    
    
    p_bias_se <- ggplot(df_a, aes(x = lambda_scaler)) + 
      geom_line(aes(y = bias_se_ratio, color = "Empirical")) + 
      geom_point(aes(y = bias_se_ratio, color = "Empirical")) +
      geom_line(aes(y = bias_se_ratio_bt, color = "Bootstrap")) + 
      geom_point(aes(y = bias_se_ratio_bt, color = "Bootstrap")) +
      geom_line(aes(y = oracal_bias_se_ratio, color = "Oracal")) + 
      geom_point(aes(y = oracal_bias_se_ratio, color = "Oracal")) +
      geom_hline(aes(yintercept=1/log(n))) +
      scale_color_manual(name='method',
                         breaks=c('Oracal', 'Empirical', 'Bootstrap'),
                         values=c('Oracal'='darkolivegreen3', 'Empirical'='lightsalmon', 'Bootstrap'='darkorchid')) +
      scale_x_continuous(breaks=seq(0,1.2,by=0.25)) +
      theme_bw() +
      theme(axis.title=element_blank(),
            legend.position='none')
    
    if(!is.na(max_bias_sd)){
      p_bias_se <- p_bias_se + ylim(0,max_bias_sd)
    }
    
    p_cr <- ggplot(df_a, aes(x = lambda_scaler)) +  
      geom_rect(data=NULL,aes(xmin=-Inf,xmax=Inf,ymin=0.95,ymax=Inf), fill="khaki1", alpha = 0.1)+ # fill="darkseagreen1"
      geom_line(aes(y = cover_rate, color = "Empirical")) + 
      geom_point(aes(y = cover_rate, color = "Empirical")) + 
      geom_line(aes(y = cover_rate_bt, color = "Bootstrap")) + 
      geom_point(aes(y = cover_rate_bt, color = "Bootstrap")) +
      geom_line(aes(y = oracal_cover_rate, color = "Oracal")) + 
      geom_point(aes(y = oracal_cover_rate, color = "Oracal")) +
      scale_color_manual(name='method',
                         breaks=c('Oracal', 'Empirical', 'Bootstrap'),
                         values=c('Oracal'='darkolivegreen3', 'Empirical'='lightsalmon', 'Bootstrap'='darkorchid')) +
      scale_x_continuous(breaks=seq(0,1.2,by=0.25)) +
      theme_bw()  + 
      theme(axis.title=element_blank(),
            legend.position='none')
    
    if(!any(is.na(u_g_scaler) & is.na(u_l_scalers)) ){
      p_est_avg <- p_est_avg +      
        geom_vline(xintercept = u_g_scaler, lty=2, col = "#00BA38") +
        geom_vline(xintercept = u_l_scaler, lty=2, col = "#619CFF") + 
        geom_vline(xintercept = 1, lty=2, col = "#F8766D") 
      
      p_bias <- p_bias +      
        geom_vline(xintercept = u_g_scaler, lty=2, col = "#00BA38") +
        geom_vline(xintercept = u_l_scaler, lty=2, col = "#619CFF") + 
        geom_vline(xintercept = 1, lty=2, col = "#F8766D") 
      
      p_se <- p_se +      
        geom_vline(xintercept = u_g_scaler, lty=2, col = "#00BA38") +
        geom_vline(xintercept = u_l_scaler, lty=2, col = "#619CFF") + 
        geom_vline(xintercept = 1, lty=2, col = "#F8766D") 
      
      p_bias_se <- p_bias_se +      
        geom_vline(xintercept = u_g_scaler, lty=2, col = "#00BA38") +
        geom_vline(xintercept = u_l_scaler, lty=2, col = "#619CFF") + 
        geom_vline(xintercept = 1, lty=2, col = "#F8766D") 
      
      p_cr <- p_cr +      
        geom_vline(xintercept = u_g_scaler, lty=2, col = "#00BA38") +
        geom_vline(xintercept = u_l_scaler, lty=2, col = "#619CFF") + 
        geom_vline(xintercept = 1, lty=2, col = "#F8766D") 
      
        
    }
    
    p_est_avg_list[[i]] = p_est_avg
    p_bias_list[[i]] = p_bias
    p_se_list[[i]] = p_se
    p_bias_se_list[[i]] = p_bias_se
    p_cr_list[[i]] = p_cr
  }

  
  g1 <- arrangeGrob(grobs = p_est_avg_list, nrow=1, left = grid::textGrob("Estimation, 95% CI", rot=90, gp=gpar(fontsize=12)))
  g2 <- arrangeGrob(grobs = p_bias_list, nrow=1, left = grid::textGrob("|Bias|", rot=90, gp=gpar(fontsize=12)))
  g3 <- arrangeGrob(grobs = p_se_list, nrow=1, left = grid::textGrob("Standard Error", rot=90, gp=gpar(fontsize=12)))
  g4 <- arrangeGrob(grobs = p_bias_se_list, nrow=1, left = grid::textGrob("|Bias| / Standard Error", rot=90, gp=gpar(fontsize=12)))
  g5 <- arrangeGrob(grobs = p_cr_list, nrow=1, left = grid::textGrob("Coverage rate", rot=90, gp=gpar(fontsize=12)),
                    bottom = grid::textGrob("lambda scalers", gp=gpar(fontsize=15)))
  
  p <- grid.arrange(g1, g2, g3, g4, g5, legend, legend_undersmoothing, 
                    layout_matrix = rbind(c(1,NA),
                                          c(2,7),
                                          c(3,6),
                                          c(4,NA),
                                          c(5,NA)),
                    widths=c(13, 1), 
                    top = textGrob("HAL-based plug-in estimator performances for E[Y|a,W] \n", 
                                   gp=gpar(fontsize=18)))  
  
  if(!is.na(save_plot)){
    ggsave(save_plot, plot=p, width = 25, height = 9, dpi = 500)
  }
  
  return(p)
}

