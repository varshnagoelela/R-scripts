#/**
#* ##############################################################
#* get the date and time in a specific format (20120810_110426).
#* Created by Varshna Goelela
#* Updated: 2012 August 10
#* ##############################################################
#**/

#System Date for filenaming
dateTime <- function (x) {
  dateTIME= format(Sys.time(), "%Y%m%d_%H%M%S")
}