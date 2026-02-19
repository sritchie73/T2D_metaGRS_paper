library(data.table)
library(lubridate)

# Function to compute years between two dates - note this preserves
# human notions of whole years
years_between <- function(d1, d2) {
  as.period(interval(as.Date(d1), as.Date(d2)), unit="years") / years(1)
}

# Likewise, add_years is consistent with the above, i.e.
add_years <- function(d1, follow) {
  as.IDate(as.Date(d1) + years(follow))
}

# Get the midpoint between two dates
midpoint <- function(d1, d2) {
  m <- as.IDate(length(d1))
  not_na <- !is.na(d1) & !is.na(d2)
  m[not_na] <- as.IDate(date_decimal((decimal_date(d2[not_na]) - decimal_date(d1[not_na]))/2 + decimal_date(d1[not_na])))
  return(m)
}
