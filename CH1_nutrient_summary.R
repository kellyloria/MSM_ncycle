bg_nuts <- readRDS("/Users/kellyloria/Documents/LittoralMetabModeling/RawData/WaterChem/NS_chem_dat_nh4_24.rds") %>%
  #mutate(date = as.Date(date, format="%m/%d/%y")) %>%
  filter(location=="stream")

bg_nuts_bwL <- bg_nuts%>%filter(site=="BWL")
range(bg_nuts_bwL$date)

bg_nuts_bwU <- bg_nuts%>%filter(site=="BWU")
range(bg_nuts_bwU$date)

bg_nuts_gbL <- bg_nuts%>%filter(site=="GBL")
range(bg_nuts_gbL$date)