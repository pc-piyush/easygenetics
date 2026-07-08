setwd("S:/Epi HST Staff/Data Management and Analysis/Martin-Spring 2020/Files")
data1 <- read.csv("owid-covid-data_5_countries.csv", header = TRUE)

data1$date1 <- format(as.Date(data1$date, format = "%m/%d/%Y"),"%Y-%m")

#Summarize d/m/yr to m/yr
data_grouped <- data1 %>%
  group_by(date1,location) %>%
  summarize(new_cases = sum(new_cases))



#Sum total_cases and new_cases
data_sum <- data1 %>% group_by(location) %>%
  summarize(total_cases = sum(total_cases),
            new_cases = sum(new_cases))

data_sum <- data1 %>% group_by(location)

#remove columns
data1 <- data1[, -c(5:17)]

?source

