plot.EpiInvert <- function(x, what = c("all", "incid", "R", "SI"),date_start="1000-01-01",date_end="3000-01-01"){

  library(ggplot2)
  library(dplyr)
  library(gridExtra)
  library(grid)

  if(what == "R"){

    d2 <- data.frame(
      date = as.Date(x$dates),
      Rt = x$Rt,
      Rt_CI95 = x$Rt_CI95
    )
    d <- filter(d2,date >= as.Date(date_start) & date <= as.Date(date_end))

    g <- ggplot(d, aes(x=date, y=Rt)) +
      geom_line() +
      scale_x_date(date_labels = '%Y-%m-%d') +
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
      geom_ribbon(aes(ymin = Rt - Rt_CI95, ymax = Rt + Rt_CI95), fill=I(rgb(0.7, 0.7, 0.7)))+
      geom_line(aes(y = Rt))+geom_hline(yintercept=1)+
      geom_text(x = -Inf, y = Inf, label='Rt',hjust = -0.1, vjust = 1.5)

    #ggplot(d, aes(x=date, y=Rt)) + geom_line() +
    #geom_ribbon(aes(ymin = Rt - Rt_CI95, ymax = Rt + Rt_CI95), fill=I(rgb(0.7, 0.7, 0.7)))+geom_line(aes(y = Rt))+geom_hline(yintercept=1)
    #scale_x_date(date_labels = '%Y-%m-%d')+theme(axis.title.y=element_blank())+
    #labs(x='', y='normalized noise') +
    # geom_text(x = -Inf, y = Inf, label=paste('Rt'),hjust = -0.1, vjust = 1.5)

    p <- ggplotGrob(g)
    grid.draw(rbind(p))

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

    d <- filter(d2,date >= as.Date(date_start) & date <= as.Date(date_end))

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

    g <- ggplot(m,aes(date,values,col=incidence))+geom_line()+
      scale_color_manual(values=c("black","green","blue","red"))+
      scale_x_date(date_labels = '%Y-%m-%d')+theme(axis.title.y=element_blank(),axis.title.x=element_blank())+geom_line(size=0.7)


    p <- ggplotGrob(g)
    grid.draw(rbind(p))
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


    g<-ggplot(dsi, aes(x=days, y=si_distr)) + geom_line() +theme(axis.title.y=element_blank())+labs(x='days', y='') +
      geom_text(x = -Inf, y = Inf, label='serial interval distribution',hjust = -0.1, vjust = 1.5)
    p <- ggplotGrob(g)
    grid.draw(rbind(p))

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

    d <- filter(d2,date >= as.Date(date_start) & date <= as.Date(date_end))

    g1 <- ggplot(d, aes(x=date, y=i_original)) + geom_line() + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.y=element_blank())+
      scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
      geom_text(x = -Inf, y = Inf, label='original incidence',hjust = -0.1, vjust = 1.5)
    g2 <- ggplot(d, aes(x=date, y=i_festive)) + geom_line() + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.y=element_blank())+
      scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) + geom_text(x = -Inf, y = Inf, label='festive bias free incidence',hjust = -0.1, vjust = 1.5)
    g3 <- ggplot(d, aes(x=date, y=i_bias_free)) + geom_line() + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.y=element_blank())+ scale_y_continuous(labels = function(x) format(x, scientific = TRUE))+ geom_text(x = -Inf, y = Inf, label='weekly + festive biases free incidence',hjust = -0.1, vjust = 1.5)
    g4 <- ggplot(d, aes(x=date, y=i_restored)) + geom_line() + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.y=element_blank())+ scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) + geom_text(x = -Inf, y = Inf, label='restored incidence',hjust = -0.1, vjust = 1.5)
    g5 <- ggplot(d, aes(x=date, y=Rt)) + geom_line() + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.y=element_blank())+geom_ribbon(aes(ymin = Rt - Rt_CI95, ymax = Rt + Rt_CI95), fill=I(rgb(0.7, 0.7, 0.7)))+geom_line(aes(y = Rt))+geom_hline(yintercept=1)+ geom_text(x = -Inf, y = Inf, label='Rt',hjust = -0.1, vjust = 1.5)
    g6 <- ggplot(d, aes(x=date, y=seasonality)) + geom_line() + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.y=element_blank()) + geom_text(x = -Inf, y = Inf, label='seasonality',hjust = -0.1, vjust = 1.5)
    g7 <- ggplot(d, aes(x=date, y=epsilon)) + geom_line() + scale_x_date(date_labels = '%Y-%m-%d')+theme(axis.title.y=element_blank())+labs(x='', y='normalized noise') + geom_text(x = -Inf, y = Inf, label=paste('normalized noise ( a = ',trunc(100*x$power_a)/100,')'),hjust = -0.1, vjust = 1.5)


    p1 <- ggplotGrob(g1)
    p2 <- ggplotGrob(g2)
    p3 <- ggplotGrob(g3)
    p4 <- ggplotGrob(g4)
    p5 <- ggplotGrob(g5)
    p6 <- ggplotGrob(g6)
    p7 <- ggplotGrob(g7)
    grid.draw(rbind(p1,p2,p3,p4,p5,p6,p7))
  }
  else{
    stop('The drawing option is not recognized')
  }


}
