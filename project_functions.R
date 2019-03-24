## Scrape dividend yield for the S&P 500 (Daily data in %)
div_yield_extraction <- function (div_yield_URL) {
  
  
  ## Read HTML DOM
  div_yield_page  <- read_html(div_yield_URL)
  
  ## Remove inner labels compromising table structure
  span_nodes      <- div_yield_page %>% html_nodes(xpath = '//*[contains(@class,"value") or contains(@class,"estimate")]')
  xml_remove(span_nodes)
  
  ## Fetch S&P 500 dividend yield table
  div_yield_table <- div_yield_page %>% html_node("#datatable") %>% html_table()
  
  ## Dividend Yield
  div_yield <- as.numeric(sub("%", "", div_yield_table$Yield[1], fixed = T)) / 100
  
  return(div_yield)
  
  
}


## Fetch closest expiration date from quantmod
fbl_next_expdate <- function(Symbols) {
  
  checkExp <- !hasArg(".expiry.known") || !match.call(expand.dots = TRUE)$.expiry.known
  urlExp <- paste0("https://query2.finance.yahoo.com/v7/finance/options/", 
                   Symbols[1])
  if (!checkExp) 
    urlExp <- paste0(urlExp, "?&date=", Exp)
  tbl <- jsonlite::fromJSON(urlExp)
  all.expiries <- tbl$optionChain$result$expirationDates[[1]]
  all.expiries.posix <- .POSIXct(as.numeric(all.expiries), 
                                 tz = "UTC")
  
  return(as.Date(all.expiries.posix[1]))
  
}
