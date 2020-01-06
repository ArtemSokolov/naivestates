## Top-scoring cells in naivestates
##
## by Artem Sokolov

library( tidyverse )
library( naivestates )

fnData <- system.file( "examples/example1_data.csv.gz", package="naivestates" )
X <- read_csv( fnData, col_types=cols() ) %>%
    rename_at( vars(starts_with("Cell_25546ON")), str_sub, 13 )

fnMap <- system.file( "examples/example1_chnlmap.csv", package="naivestates" )
M <- read_csv( fnMap, col_types=cols() ) %>% mutate_at( "Channel", str_sub, 13 )

Fits <- GMMfit( X, CellId, M$Channel, qq=0.005 )
exprPostProb <- GMMreshape( Fits )

task1 <- c("Immune", "Tumor", "Stroma")
R <- taskPostProb( exprPostProb, CellId, M, task1 ) %>% na.omit()
XR <- inner_join( X, R, by="CellId" ) %>% slice( seq(1, n(), by=5) )

## Isolate S100 fit
S100mn <- 7; S100mx <- 11
FS100 <- filter( Fits, Marker == "S100" ) %>% pluck( "Values", 1 ) %>%
    filter( Value >= S100mn, Value <= S100mx ) %>% arrange( Value )

## Isolate CATENIN fit
CATmn <- 7.25; CATmx <- 10
FCAT <- filter( Fits, Marker == "CATENIN" ) %>% pluck( "Values", 1 ) %>%
    filter( Value >= CATmn, Value <= CATmx ) %>% arrange( Value )

## Short-hand elements
etxt <- function(s, ...) {element_text( size = s, face = "bold", ... )}
eblk <- function() element_blank()

## Themes
theme_nox <- function() {
    theme_bw() +
        theme(axis.title.x=eblk(), axis.text.x=eblk(), axis.ticks.x=eblk(),
              axis.title.y=etxt(14), axis.text.y=etxt(12) ) }

theme_noy <- function() {
    theme_bw() +
        theme(axis.title.y=eblk(), axis.text.y=eblk(), axis.ticks.y=eblk(),
              axis.title.x=etxt(14), axis.text.x=etxt(12)) }

## Main plot
gg0 <- ggplot( X1, aes(x=S100, y=CATENIN, color=Call) ) +
    geom_point() + theme_bw() + ggthemes::scale_color_few() +
    xlim( c(S100mn, S100mx) ) + ylim( c(CATmn, CATmx) ) +
    guides( color=guide_legend(override.aes = list(size=5)) ) +
    theme( axis.title=etxt(14), axis.text=etxt(12),
          legend.title=etxt(16), legend.text=etxt(14) )

## Subpanels
ggl <- cowplot::get_legend(gg0)
gg1 <- gg0 + guides(color=FALSE)
gg2 <- ggplot( FS100, aes(x=Value, y=Prob) ) + geom_line(lwd=1.1) + theme_nox()
gg3 <- ggplot( FCAT, aes(x=Value, y=Prob) ) + geom_line(lwd=1.1) + theme_noy() + coord_flip()

## Overall plot
gg <- egg::ggarrange( gg2, cowplot::ggdraw(ggl), gg1, gg3, ncol=2,
                     heights=c(1,2), widths=c(2,1) )
ggsave( "tumor_call.png", gg, width=9, height=7 )
