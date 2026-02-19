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
