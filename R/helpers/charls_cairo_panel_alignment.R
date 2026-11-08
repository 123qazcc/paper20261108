# CHARLS-only adaptation: measure actual grid viewports on the export Cairo device.
# Existing HRS/ELSA geometry helpers and output files are not modified.
source("private/render_tools/nature-figure/scripts/panel_alignment.R")
require_charls_panel_alignment <- function(plot,stem,width_in,height_in,family="Arial") {
  probe<-tempfile(fileext=".pdf")
  grDevices::cairo_pdf(probe,width=width_in,height=height_in,family=family)
  open<-TRUE
  on.exit({if(open)grDevices::dev.off();unlink(probe)},add=TRUE)
  grob<-if(inherits(plot,"patchwork"))patchwork::patchworkGrob(plot)else ggplot2::ggplotGrob(plot)
  cells<-grob$layout[grepl("^panel(?:-[0-9]+)*$",grob$layout$name,perl=TRUE),,drop=FALSE]
  cells<-cells[order(cells$t,cells$l),,drop=FALSE]
  stopifnot(nrow(cells)>0L)
  grid::grid.newpage();grid::grid.draw(grob);grid::grid.force()
  panels<-lapply(seq_len(nrow(cells)),function(i){
    z<-cells[i,];vp<-paste0(z$name,".",z$t,"-",z$l,"-",z$b,"-",z$r)
    grid::seekViewport(vp)
    lo<-grid::deviceLoc(grid::unit(0,"npc"),grid::unit(0,"npc"),valueOnly=TRUE)
    hi<-grid::deviceLoc(grid::unit(1,"npc"),grid::unit(1,"npc"),valueOnly=TRUE)
    grid::upViewport(0);rect<-c(lo$x,lo$y,hi$x,hi$y)*72
    stopifnot(all(is.finite(rect)),rect[3]>rect[1],rect[4]>rect[2])
    list(id=letters[i],bbox_pt=rect,grid_id="charls-rendered-grid",
      row_start=as.integer(z$t-1L),row_stop=as.integer(z$b),
      col_start=as.integer(z$l-1L),col_stop=as.integer(z$r))
  })
  grDevices::dev.off();open<-FALSE
  manifest<-list(schema_version=1L,backend="r-cairo-rendered-viewports",
    figure=list(width_pt=width_in*72,height_pt=height_in*72),panels=panels,exemptions=list())
  jsonlite::write_json(manifest,paste0(stem,".alignment-layout.json"),auto_unbox=TRUE,pretty=TRUE,digits=NA)
  python<-"python"
  auditor<-"private/render_tools/nature-figure/scripts/audit_panel_alignment.py"
  status<-system2(python,c(shQuote(auditor),shQuote(paste0(stem,".alignment-layout.json")),
    "--json-out",shQuote(paste0(stem,".alignment.json")),"--strict"))
  if(status!=0L)stop("Rendered panel alignment did not pass: ",stem)
  invisible(manifest)
}
