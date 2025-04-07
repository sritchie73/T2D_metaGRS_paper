maf <- function(eaf) {
  0.5 - abs(eaf - 0.5)  # i.e. 0.99 and 0.01 both become 0.01 
}
