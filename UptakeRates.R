rawdat = list.files(paste("./BTC_out/",sep=""), full.names = T)


# Lets start calculating uptake rates

# GB
GBL_06 <- read.csv("./BTC_out/GBL_BTC_20210603.csv")

model.drm1 <- drm (Uadd ~ NO3, data = GBL_06, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NO3 = seq(0, max(GBL_06$NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

summary(model.drm1)
(model.drm1)

Uadd_plotGB <- ggplot(GBL_06, aes(x=NO3, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#a67d17") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("Glenbrook L Creek") +
  theme_bw()



# GB
GBL_06 <- read.csv("./BTC_out/GBL_NO3_BTC_220623.csv")

model.drm1 <- drm (Uadd ~ NO3, data = GBL_06, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NO3 = seq(0, max(GBL_06$NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

summary(model.drm1)
(model.drm1)

Uadd_plotGB <- ggplot(GBL_06, aes(x=NO3, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#a67d17") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("Glenbrook L Creek") +
  theme_bw()





# GB
GBL_06 <- read.csv("./BTC_out/GBL_BTC_NH4_20220623.csv")

model.drm1 <- drm (Uadd ~ NH4, data = GBL_06, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NO3 = seq(0, max(GBL_06$NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

summary(model.drm1)
(model.drm1)

Uadd_plotGB <- ggplot(GBL_06, aes(x=NH4, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#a67d17") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("Glenbrook L Creek") +
  theme_bw()


GBL_07 <- read.csv("./BTC_out/GB_BTC_20210722.csv")
summary(GBL_07)

model.drm1 <- drm (Uadd ~ NO3, data = GBL_07a, fct = MM.2())
summary(model.drm1)



summary(model.drm1)
(model.drm1)

Uadd_plotGB <- ggplot(GBL_07, aes(x=NO3, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
 # geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#a67d17") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("Glenbrook L Creek") +
  theme_bw()




GBL_12 <- read.csv("./BTC_out/GBL_BTC_221212.csv")
summary(GBL_12)
GBL_12a<-GBL_12[c(-2),]

model.drm1 <- drm (Uadd ~ NH4, data = GBL_12a, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame( NH4= seq(0, max(GBL_12a$NH4), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)


summary(model.drm1)
(model.drm1)

Uadd_plotGB <- ggplot(GBL_12a, aes(x=NH4, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  #geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#a67d17") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("Glenbrook L Creek") +
  theme_bw()


GBL_10 <- read.csv("./BTC_out/GBL_NO3_BTC_221003.csv")
summary(GBL_10)
GBL_10a<-GBL_10[c(-1),]

model.drm1 <- drm (Uadd ~ NO3, data = GBL_10a, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame( NO3= seq(0, max(GBL_10a$NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

Uadd_plotGB <- ggplot(GBL_10a, aes(x=NO3, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#a67d17") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("Glenbrook L Creek") +
  theme_bw()


GBL_10 <- read.csv("./BTC_out/GBL_BTC_NH4_20221003.csv")
GBL_10a<-GBL_10[c(-1,-2),]

model.drm1 <- drm (Uadd ~ NH4, data = GBL_10a, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame( NH4= seq(0, max(GBL_12a$NH4), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

Uadd_plotGB <- ggplot(GBL_10a, aes(x=NH4, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  #geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#a67d17") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("Glenbrook L Creek") +
  theme_bw()

#############
GBL_11 <- read.csv("./BTC_out/GBL_BTC_NH4_20221104.csv")
GBL_11a<-GBL_11[c(-1,-2),]

model.drm1 <- drm (Uadd ~ NH4, data = GBL_11a, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame( NH4= seq(0, max(GBL_12a$NH4), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

Uadd_plotGB <- ggplot(GBL_11a, aes(x=NH4, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  #geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#a67d17") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("Glenbrook L Creek") +
  theme_bw()

################
GBL_12 <- read.csv("./BTC_out/GBL_NO3_BTC_221212.csv")
GBL_10a<-GBL_10[c(-1,-2),]

model.drm1 <- drm (Uadd ~ NO3, data = GBL_12, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame( NO3= seq(0, max(GBL_12$NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

Uadd_plotGB <- ggplot(GBL_12, aes(x=NO3, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#a67d17") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("Glenbrook L Creek") +
  theme_bw()


GBL_04 <- read.csv("./BTC_out/GBL_NO3_BTC_220407.csv")
summary(GBL_04)

model.drm1 <- drm (Uadd ~ NO3, data = GBL_04, fct = MM.2())
summary(model.drm1)

summary(model.drm1)
(model.drm1)

Uadd_plotGB <- ggplot(GBL_04, aes(x=NO3, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  #geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#a67d17") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  #labs(x = expression(paste("Nitrate (mgL)")),
  #     y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("Glenbrook L Creek") +
  theme_bw()


GBL_04 <- read.csv("./BTC_out/GBL_BTC_NH4_20220407.csv")
summary(GBL_04)

model.drm1 <- drm (Uadd ~ NH4, data = GBL_04, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame( NH4= seq(0, max(GBL_12a$NH4), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)


summary(model.drm1)
(model.drm1)

Uadd_plotGB <- ggplot(GBL_04, aes(x=NH4, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  #geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#a67d17") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  #labs(x = expression(paste("Nitrate (mgL)")),
  #     y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("Glenbrook L Creek") +
  theme_bw()



GBL_12 <- read.csv("./BTC_out/GBL_NO3_BTC_221212.csv")
summary(GBL_12)

GBL_12a<-GBL_12[c(-2),]

model.drm1 <- drm (Uadd ~ NO3, data = GBL_12a, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame( NO3= seq(0, max(GBL_12a$NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)


summary(model.drm1)
(model.drm1)

Uadd_plotGB <- ggplot(GBL_12a, aes(x=NO3, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#a67d17") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  #labs(x = expression(paste("Nitrate (mgL)")),
  #     y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("Glenbrook L Creek") +
  theme_bw()


GBL_10 <- read.csv("./BTC_out/GBL_NO3_BTC_BWL221003.csv")
summary(GBL_10)
GBL_10a<-na.omit(GBL_10[c(-1,-2),])

model.drm1 <- drm (Uadd ~ NO3, data = GBL_10a, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame( NO3= seq(0, max(GBL_10a$NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)


summary(model.drm1)
(model.drm1)

Uadd_plotGB <- ggplot(GBL_10a, aes(x=NO3, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#a67d17") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  #labs(x = expression(paste("Nitrate (mgL)")),
  #     y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("Glenbrook L Creek") +
  theme_bw()




GBU_10 <- read.csv("./BTC_out/GBU_BTC_20221003_MW.csv")
summary(GBU_10)

GBU_10a<-na.omit(GBU_10[c(-11,-13),])

model.drm1 <- drm (Uadd ~ NH4, data = GBU_10a, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame( NH4= seq(0, max(GBU_10a$NH4), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)


summary(model.drm1)
(model.drm1)

Uadd_plotGB <- ggplot(GBU_10a, aes(x=NH4, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NH4, y = Uadd), colour = "#a67d17") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  #labs(x = expression(paste("Nitrate (mgL)")),
  #     y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("Glenbrook U Creek") +
  theme_bw()



GBU_11 <- read.csv("./BTC_out/GBU_BTC_20221104_MW.csv")
summary(GBU_11)

GBU_10a<-na.omit(GBU_10[c(-11,-13),])

model.drm1 <- drm (Uadd ~ NH4, data = GBU_10a, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame( NH4= seq(0, max(GBU_10a$NH4), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)


summary(model.drm1)
(model.drm1)

Uadd_plotGB <- ggplot(GBU_11, aes(x=NH4, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  #geom_line(data = mm2, aes(x = NH4, y = Uadd), colour = "#a67d17") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  #labs(x = expression(paste("Nitrate (mgL)")),
  #     y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("Glenbrook U Creek") +
  theme_bw()





GBU_06 <- read.csv("./BTC_out/GBU_BTC_NH4_20220623.csv")
summary(GBU_06)
GBU_06a<-na.omit(GBU_06[c(-2),])

model.drm1 <- drm (Uadd ~ NH4, data = GBU_06, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame( NH4= seq(0, max(GBU_06$NH4), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)


Uadd_plotGB <- ggplot(GBU_06, aes(x=log(NH4+1), y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  #geom_line(data = mm2, aes(x = NH4, y = Uadd), colour = "#a67d17") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  #labs(x = expression(paste("Nitrate (mgL)")),
  #     y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("Glenbrook U Creek") +
  theme_bw()





GBU_04 <- read.csv("./BTC_out/GBU_NO3_BTC_220407.csv")
summary(GBU_04)

model.drm1 <- drm (Uadd ~ NH4, data = GBU_10a, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame( NH4= seq(0, max(GBU_10a$NH4), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)


summary(model.drm1)
(model.drm1)

Uadd_plotGB <- ggplot(GBU_04, aes(x=NO3, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  #geom_line(data = mm2, aes(x = NH4, y = Uadd), colour = "#a67d17") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  #labs(x = expression(paste("Nitrate (mgL)")),
  #     y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("Glenbrook U Creek") +
  theme_bw()



GBU_06 <- read.csv("./BTC_out/GBU_NO3_BTC_220623.csv")
summary(GBU_06)

model.drm1 <- drm (Uadd ~ NO3, data = GBU_06, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame( NH4= seq(0, max(GBU_10a$NH4), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)


summary(model.drm1)
(model.drm1)

Uadd_plotGB <- ggplot(GBU_06, aes(x=NO3, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  #geom_line(data = mm2, aes(x = NH4, y = Uadd), colour = "#a67d17") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  #labs(x = expression(paste("Nitrate (mgL)")),
  #     y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("Glenbrook U Creek") +
  theme_bw()


# GB
GBL_03 <- read.csv("./BTC_out/GBL_BTC_NH4_20230327.csv")

model.drm1 <- drm (Uadd ~ NH4, data = GBL_03, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NO3 = seq(0, max(GBL_06$NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

summary(model.drm1)
(model.drm1)

Uadd_plotGB <- ggplot(GBL_03, aes(x=NH4, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#a67d17") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("Glenbrook L Creek") +
  theme_bw()


GBL_03 <- read.csv("./BTC_out/GBL_NO3_BTC_230327.csv")
summary(GBL_03)
GBL_10a<-na.omit(GBL_10[c(-1,-2),])

model.drm1 <- drm (Uadd ~ NO3, data = GBL_03, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame( NO3= seq(0, max(GBL_10a$NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

Uadd_plotGB <- ggplot(GBL_03, aes(x=NO3, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#a67d17") +
  geom_point(size = 3, shape= 17, col = "#a67d17") +
  ggtitle("Glenbrook L Creek") +
  theme_bw()












###
###
###

# BWL
BWL_05 <- read.csv("./BTC_out/BWL_NO3_BTC_BWL220526.csv")

model.drm1 <- drm (Uadd ~ NO3, data = BWL_05, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NO3 = seq(0, max(BWL_05$NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

Uadd_plotBW <- ggplot(BWL_05, aes(x=NO3, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  #geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#3283a8") +
  geom_point(size = 3, shape= 17, col = "#3283a8",) +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("BWL") +
  theme_bw()

BWL_05 <- read.csv("./BTC_out/BWL_NH4_BTC_BWL220526.csv")
BWL_05a<-BWL_05[c(-8),]

model.drm1 <- drm (Uadd ~ NH4, data = BWL_05a, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NH4 = seq(0, max(BWL_05a$NH4), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

Uadd_plotBW <- ggplot(BWL_05a, aes(x=NH4, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NH4, y = Uadd), colour = "#3283a8") +
  geom_point(size = 3, shape= 17, col = "#3283a8",) +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("BWL") +
  theme_bw()


BWL_08 <- read.csv("./BTC_out/BWL_NH4_BTC_BWL220824.csv")
#BWL_05a<-BWL_05[c(-8),]
model.drm1 <- drm (Uadd ~ NH4, data = BWL_08, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NH4 = seq(0, max(BWL_08$NH4), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

Uadd_plotBW <- ggplot(BWL_08, aes(x=NH4, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NH4, y = Uadd), colour = "#3283a8") +
  geom_point(size = 3, shape= 17, col = "#3283a8",) +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("BWL") +
  theme_bw()


BWL_10 <- read.csv("./BTC_out/BWL_NH4_BTC_BWL221012.csv")
#BWL_05a<-BWL_05[c(-8),]
model.drm1 <- drm (Uadd ~ NH4, data = BWL_10, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NH4 = seq(0, max(BWL_10$NH4), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

Uadd_plotBW <- ggplot(BWL_10, aes(x=NH4, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NH4, y = Uadd), colour = "#3283a8") +
  geom_point(size = 3, shape= 17, col = "#3283a8",) +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("BWL") +
  theme_bw()


BWL_10 <- read.csv("./BTC_out/BWL_NO3_BTC_BWL221012.csv")
BWL_10a<-BWL_10[c(-1,-2),]

model.drm1 <- drm (Uadd ~ NO3, data = BWL_10, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NO3 = seq(0, max(BWL_10$NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

Uadd_plotBW <- ggplot(BWL_10, aes(x=NO3, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#3283a8") +
  geom_point(size = 3, shape= 17, col = "#3283a8",) +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("BWL") +
  theme_bw()

###

BWL_11 <- read.csv("./BTC_out/BWL_NH4_BTC_BWL221121v2.csv")
BWL_11a<-BWL_11[c(-1, -2, -3),]
model.drm1 <- drm (Uadd ~ NH4, data = BWL_11a, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NH4 = seq(0, max(BWL_11a$NH4), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

Uadd_plotBW <- ggplot(BWL_11a, aes(x=NH4, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NH4, y = Uadd), colour = "#3283a8") +
  geom_point(size = 3, shape= 17, col = "#3283a8",) +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("BWL") +
  theme_bw()

BWL_11 <- read.csv("./BTC_out/BWL_NO3_BTC_BWL221121.csv")
BWL_11a<-BWL_11[c(-1,-2,-15, -14),]

model.drm1 <- drm (Uadd ~ NO3, data = BWL_11a, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NO3 = seq(0, max(BWL_11a$NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

Uadd_plotBW <- ggplot(BWL_11a, aes(x=NO3, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#3283a8") +
  geom_point(size = 3, shape= 17, col = "#3283a8",) +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("BWL") +
  theme_bw()

###
BWL_12 <- read.csv("./BTC_out/BWL_NH4_BTC_BWL221219.csv")
BWL_12a<-BWL_12[c(-1,-2,-3,-4,-18),]
model.drm1 <- drm (Uadd ~ NH4, data = BWL_12a, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NH4 = seq(0, max(BWL_12a$NH4), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

Uadd_plotBW <- ggplot(BWL_12a, aes(x=NH4, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NH4, y = Uadd), colour = "#3283a8") +
  geom_point(size = 3, shape= 17, col = "#3283a8",) +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("BWL") +
  theme_bw()


BWL_12 <- read.csv("./BTC_out/BWL_NO3_BTC_BWL221219.csv")
BWL_12a<-BWL_12[c(-1,-2,-14,-15),]

model.drm1 <- drm (Uadd ~ NO3, data = BWL_12a, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NO3 = seq(0, max(BWL_12a$NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

Uadd_plotBW <- ggplot(BWL_12a, aes(x=NO3, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#3283a8") +
  geom_point(size = 3, shape= 17, col = "#3283a8",) +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("BWL") +
  theme_bw()

###

BWL_02 <- read.csv("./BTC_out/BWL_NH4_BTC_BWL230215.csv")
BWL_12a<-BWL_12[c(-1,-2,-3,-4,-18),]
model.drm1 <- drm (Uadd ~ NH4, data = BWL_02, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NH4 = seq(0, max(BWL_12a$NH4), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

Uadd_plotBW <- ggplot(BWL_02, aes(x=NH4, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NH4, y = Uadd), colour = "#3283a8") +
  geom_point(size = 3, shape= 17, col = "#3283a8",) +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("BWL") +
  theme_bw()


BWL_02 <- read.csv("./BTC_out/BWL_NO3_BTC_BWL230215.csv")
BWL_12a<-BWL_12[c(-1,-2,-14,-15),]

model.drm1 <- drm (Uadd ~ NO3, data = BWL_02, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NO3 = seq(0, max(BWL_02$NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

Uadd_plotBW <- ggplot(BWL_02, aes(x=NO3, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#3283a8") +
  geom_point(size = 3, shape= 17, col = "#3283a8",) +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("BWL") +
  theme_bw()

####
BWL_04 <- read.csv("./BTC_out/BWL_NH4_BTC_BWL230405.csv")
BWL_12a<-BWL_12[c(-1,-2,-3,-4,-18),]
model.drm1 <- drm (Uadd ~ NH4, data = BWL_04, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NH4 = seq(0, max(BWL_04$NH4), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

Uadd_plotBW <- ggplot(BWL_04, aes(x=NH4, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NH4, y = Uadd), colour = "#3283a8") +
  geom_point(size = 3, shape= 17, col = "#3283a8",) +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("BWL") +
  theme_bw()


BWL_04 <- read.csv("./BTC_out/BWL_NO3_BTC_BWL230405.csv")

model.drm1 <- drm (Uadd ~ NO3, data = BWL_04, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NO3 = seq(0, max(BWL_02$NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

Uadd_plotBW <- ggplot(BWL_04, aes(x=NO3, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#3283a8") +
  geom_point(size = 3, shape= 17, col = "#3283a8",) +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("BWL") +
  theme_bw()

######
#####
BWL_07 <- read.csv("./BTC_out/BWL_NH4_BTC_BWL230718.csv")
model.drm1 <- drm (Uadd ~ NH4, data = BWL_07, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NH4 = seq(0, max(BWL_04$NH4), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

Uadd_plotBW <- ggplot(BWL_07, aes(x=NH4, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  #geom_line(data = mm2, aes(x = NH4, y = Uadd), colour = "#3283a8") +
  geom_point(size = 3, shape= 17, col = "#3283a8",) +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("BWL") +
  theme_bw()

BWL_07 <- read.csv("./BTC_out/BWL_NO3_BTC_BWL230718.csv")
BWL_07a<-BWL_07[c(-12),]


model.drm1 <- drm (Uadd ~ NO3, data = BWL_07a, fct = MM.2())
summary(model.drm1)

mm2 <- data.frame(NO3 = seq(0, max(BWL_07a$NO3), length.out = 100))
mm2$Uadd <- predict(model.drm1, newdata = mm2)

Uadd_plotBW <- ggplot(BWL_07a, aes(x=NO3, y=Uadd)) +
  #ylim(0,20) + xlim(0,1500) +
  geom_line(data = mm2, aes(x = NO3, y = Uadd), colour = "#3283a8") +
  geom_point(size = 3, shape= 17, col = "#3283a8",) +
  labs(x = expression(paste("Nitrate (mgL)")),
       y= expression(paste("Uadd (mg m^-2 min^-1)"))) +
  ggtitle("BWL") +
  theme_bw()

