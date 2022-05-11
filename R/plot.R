#' @title 
#' \code{plot} the results obtained by EpiInvert
#'
#' @param x an object of class \code{estimate_EpiInvert}.
#'
#' @param what one of the following drawing options: 
#' 
#'  \itemize{
#'  \item{all}{: a plot combining the main EpiInvert results}
#'  \item{R}{: a plot of the reproduction number Rt estimation }
#'  \item{incid}{: a plot combining the obtained incidence curves }
#'  \item{SI}{: the serial interval used in the EpiInvert estimation }
#'  
#' }
#' 
#' @param date_start the start date to plot 
#' 
#' @param date_end the final date to plot 
#' 
#' @return a plot.
#' 
#' @seealso \code{\link{EpiInvert}} 


#' @export
plot <- function(x, what = "all",date_start="1000-01-01",date_end="3000-01-01"){
  
  if(what == "R"){

    d2 <- data.frame(
      date = as.Date(x$dates),
      Rt = x$Rt,
      Rt_CI95 = x$Rt_CI95
    )
   
    d <- dplyr::filter(d2,date >= as.Date(date_start,format="%Y-%m-%d") & 
                         date <= as.Date(date_end,format="%Y-%m-%d"))

    g <- ggplot2::ggplot(d, ggplot2::aes(x=date, y=Rt)) +
      ggplot2::geom_line() +
      ggplot2::ylim(0,1.02*max(d$Rt)) +
      ggplot2::scale_x_date(date_labels = '%Y-%m-%d') +
      ggplot2::theme(axis.title.x= 
      ggplot2::element_blank(),axis.title.y= ggplot2::element_blank())+
      ggplot2::geom_ribbon(ggplot2::aes(ymin = Rt - Rt_CI95, ymax = Rt + Rt_CI95), 
                           fill=I(rgb(0.7, 0.7, 0.7)))+
      ggplot2::geom_line(ggplot2::aes(y = Rt))+ 
      ggplot2::geom_hline(yintercept=1)+
      ggplot2::geom_text(x = -Inf, y = Inf, label='Rt',hjust=-0.5, vjust=1.5)

    p <- ggplot2::ggplotGrob(g)
    grid::grid.draw(rbind(p))

  }

  else if(what == "incid"){

    d2 <- data.frame(
      date = as.Date(x$dates),
      i_original = x$i_original,
      i_festive = x$i_festive,
      i_bias_free = x$i_bias_free,
      i_restored = x$i_restored,
      Rt = x$Rt,
      seasonality = x$seasonality,
      epsilon = x$epsilon,
      Rt_CI95 = x$Rt_CI95
    )

    d <-dplyr::filter(d2,date >= as.Date(date_start) & date <= as.Date(date_end))

    N <- length(d$date)
    m=data.frame( date = as.Date(d$date),
                  values=c( d$i_festive,
                            d$i_original,
                            d$i_bias_free,
                            d$i_restored
                  ),
                  incidence=c(rep("2. festive bias free",N),
                              rep("1. original",N),
                              rep("3. weekly + festive biases free",N),
                              rep("4. restored",N))
    )

    g <- ggplot2::ggplot(m,ggplot2::aes(date,values,col=incidence))+
      ggplot2::geom_line()+
      ggplot2::scale_color_manual(values=c("black","green","blue","red"))+
      ggplot2::scale_x_date(date_labels = '%Y-%m-%d')+
      ggplot2::theme(axis.title.y=ggplot2::element_blank(),
                     axis.title.x=ggplot2::element_blank())+
      ggplot2::geom_line(size=0.7)


    p <- ggplot2::ggplotGrob(g)
    grid::grid.draw(rbind(p))
  }
  else if(what == "SI"){
    days <- c(x$shift_si_distr)

    # Filling the vector using a for loop
    for(i in 2:length(x$si_distr)) {
      days <- append(days, i+x$shift_si_distr-1)
    }

    dsi <- data.frame(
      days = days,
      si_distr = x$si_distr
    )

    g<-ggplot2::ggplot(dsi, ggplot2::aes(x=days, y=si_distr)) + 
	     ggplot2::geom_line() +
       ggplot2::ylim(0,1.02*max(x$si_distr)) +
       ggplot2::theme(axis.title.y=ggplot2::element_blank())+
 	     ggplot2::labs(x='days', y='') +
       ggplot2::geom_text(x = -Inf, y = Inf, 
	              label='serial interval distribution',hjust = -0.1, vjust = 1.5)

        p <- ggplot2::ggplotGrob(g)
    grid:: grid.draw(rbind(p))
  }

  else if(what == "all"){
    d2 <- data.frame(
      date = as.Date(x$dates),
      i_original = x$i_original,
      i_festive = x$i_festive,
      i_bias_free = x$i_bias_free,
      i_restored = x$i_restored,
      Rt = x$Rt,
      seasonality = x$seasonality,
      epsilon = x$epsilon,
      Rt_CI95 = x$Rt_CI95
    )

    d <- dplyr::filter(d2,date >= as.Date(date_start) & date <= as.Date(date_end))
    
    maxI <- max(d$i_original)

    g1 <- ggplot2::ggplot(d, ggplot2::aes(x=date, y=i_original)) + 
	     ggplot2::geom_line() + 
       ggplot2::theme(axis.title.x=ggplot2::element_blank(),
	       axis.text.x=ggplot2::element_blank(),
	       axis.ticks.x=ggplot2::element_blank(),
	       axis.title.y=ggplot2::element_blank())+
       ggplot2::scale_y_continuous(labels = function(x) format(x, scientific = TRUE),
                                   limits = c(0, 1.2*maxI)) +
       ggplot2::geom_text(x = -Inf, y = Inf, 
	       label='original incidence',hjust = -0.11, vjust = 1.5)
    
  	g2 <- ggplot2::ggplot(d, ggplot2::aes(x=date, y=i_festive)) + 
  	  ggplot2::geom_line() + 
  	  ggplot2::theme(axis.title.x=ggplot2::element_blank(),
  	  axis.text.x=ggplot2::element_blank(),
  	  axis.ticks.x=ggplot2::element_blank(),
  	  axis.title.y=ggplot2::element_blank())+
  	  ggplot2::scale_y_continuous(labels = function(x) format(x, scientific = TRUE),
  	                              limits = c(0, 1.2*maxI)) + 
  	  ggplot2::geom_text(x = -Inf, y = Inf, label='festive bias free incidence',
  	  hjust = -0.11, vjust = 1.5)
  	  
    g3 <- ggplot2::ggplot(d, ggplot2::aes(x=date, y=i_bias_free)) + 
  	  ggplot2::geom_line() + 
  	  ggplot2::theme(axis.title.x=ggplot2::element_blank(),
  	  axis.text.x=ggplot2::element_blank(),axis.ticks.x=ggplot2::element_blank(),
  	  axis.title.y=ggplot2::element_blank())+ 
  	  ggplot2::scale_y_continuous(labels = function(x) format(x, scientific = TRUE),
  	                              limits = c(0, 1.2*maxI))+ 
  	  ggplot2::geom_text(x = -Inf, y = Inf, 
  	  label='weekly + festive biases free incidence',hjust = -0.07, vjust = 1.5)
  	  
    g4 <- ggplot2::ggplot(d, ggplot2::aes(x=date, y=i_restored)) + 
  	  ggplot2::geom_line() + ggplot2::theme(axis.title.x=ggplot2::element_blank(),
  	  axis.text.x=ggplot2::element_blank(),
  	  axis.ticks.x=ggplot2::element_blank(),axis.title.y=ggplot2::element_blank())+ 
  	  ggplot2::scale_y_continuous(labels = function(x) format(x, scientific = TRUE),
  	                              limits = c(0, 1.2*maxI)) + 
  	  ggplot2::geom_text(x = -Inf, y = Inf, 
  	  label='restored incidence',hjust = -0.12, vjust = 1.5)
  	  
    g5 <- ggplot2::ggplot(d, ggplot2::aes(x=date, y=Rt)) + 
      ggplot2::ylim(0,1.2*max(d$Rt)) +
  	  ggplot2::geom_line() + ggplot2::theme(axis.title.x=ggplot2::element_blank(),
  	  axis.text.x=ggplot2::element_blank(),axis.ticks.x=ggplot2::element_blank(),
  	  axis.title.y=ggplot2::element_blank())+
  	  ggplot2::geom_ribbon(ggplot2::aes(ymin = Rt - Rt_CI95, ymax = Rt + Rt_CI95), 
  	  fill=I(rgb(0.7, 0.7, 0.7)))+
  	  ggplot2::geom_line(ggplot2::aes(y = Rt))+
  	  ggplot2::geom_hline(yintercept=1)+ 
  	  ggplot2::geom_text(x = -Inf, y = Inf, label='Rt',hjust = -0.9, vjust = 1.5)
  	  
    g6 <- ggplot2::ggplot(d, ggplot2::aes(x=date, y=seasonality)) + 
      ggplot2::ylim(0,1.2*max(d$seasonality)) +
  	  ggplot2::geom_line() + ggplot2::theme(axis.title.x=ggplot2::element_blank(),
  	  axis.text.x=ggplot2::element_blank(),
  	  axis.ticks.x=ggplot2::element_blank(),
  	  axis.title.y=ggplot2::element_blank()) + 
  	  ggplot2::geom_text(x = -Inf, y = Inf, label='seasonality',hjust = -0.2, vjust = 1.5)
  	  
    g7 <- ggplot2::ggplot(d, ggplot2::aes(x=date, y=epsilon)) + 
      ggplot2::ylim(min(d$epsilon),1.2*max(d$epsilon)) +
  	  ggplot2::geom_line() + 
  	  ggplot2::scale_x_date(date_labels = '%Y-%m-%d')+
  	  ggplot2::theme(axis.title.y=ggplot2::element_blank())+
      ggplot2::labs(x='', y='normalized noise') + 
  	  ggplot2::geom_text(x = -Inf, y = Inf, label=paste('normalized noise ( a = ',trunc(100*x$power_a)/100,')'),hjust = -0.1, vjust = 1.5)


    p1 <- ggplot2::ggplotGrob(g1)
    p2 <- ggplot2::ggplotGrob(g2)
    p3 <- ggplot2::ggplotGrob(g3)
    p4 <- ggplot2::ggplotGrob(g4)
    p5 <- ggplot2::ggplotGrob(g5)
    p6 <- ggplot2::ggplotGrob(g6)
    p7 <- ggplot2::ggplotGrob(g7)
    grid::grid.draw(rbind(p1,p2,p3,p4,p5,p6,p7))
  }
  else{
    stop('The drawing option is not recognized')
  }


}
