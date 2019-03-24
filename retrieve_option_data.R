# ================================ Packages ================================


library(quantmod)
library(rvest)
library(xml2)
library(lubridate)
library(data.table)
library(purrr)

## Working Directory Setting
wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)


## Source Project Functions
source(paste0(wd, '/project_functions.R'))



# ================================ Fetch Option Prices ================================


## Date in which data is retrieved
pull_date <- as.Date("24/03/2019", format = "%d/%m/%Y")


## Turning off Yahoo Warning
options("getSymbols.yahoo.warning"=FALSE)
options("getSymbols.warning4.0"=FALSE)

## Download S&P500 prices from Quantmod
getSymbols("^GSPC", warnings = FALSE, verbose = FALSE)
stock_price <- as.numeric(tail(GSPC[, 6], 1))


## Compute next option expiration date
exp_date <- fbl_next_expdate("^GSPC") 


## Download option prices from Quantmod for the next expiration date
chain   <- getOptionChain("^GSPC", Exp = exp_date)

## Strike Prices
chain   <- rbindlist(chain)
K       <- data.frame(K = sort(unique(chain$Strike)))
optionPrices <- data.frame(K = unique(chain$Bid))


## Time to maturity
Tm <- as.numeric(exp_date - pull_date) / 360

# ================================ Fetch Div Yield and Risk Free-Rate ================================


## Download risk free rate from FRED (10-Year Treasury Constant Maturity Rate - Data in bps)
getSymbols("TB3MS", src = "FRED", warnings = FALSE, verbose = FALSE)
rf <- as.numeric(tail(TB3MS, 1)) / (100)

## Dividend Yield
div_yield <- div_yield_extraction("http://www.multpl.com/s-p-500-dividend-yield/table")

## Relevant parameters for option pricing
params <- data.frame(S0 = stock_price, r = rf, q = div_yield, "T" = Tm)

## Write to csv
fwrite(K,      file = "strike_prices.csv")
fwrite(params, file = "params.csv"       )
fwrite(optionPrices, file = "option_prices.csv")  
