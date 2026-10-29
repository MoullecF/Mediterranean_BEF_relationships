# ==============================================================================
# Figure S3: A priori conceptual path structure of the piecewise SEM
#
#   Draws the hypothesised DAG, not the fitted model
#
# ==============================================================================

# ==============================================================================
# 1. Load libraries and prepare the environment
# ==============================================================================

source(here::here("R/0-Load libraries.R"))

rm(list = ls())

# ==============================================================================
# 2. Define palette and graph, render and save
# ==============================================================================

# palette
col_driver_env   <- "#BFD8E8"   # environmental drivers
col_driver_anthro<- "#F2C4A0"   # anthropogenic drivers
col_taxo         <- "#FC8D62"   # taxonomic diversity
col_func         <- "#8DA0CB"   # functional diversity
col_trait        <- "#66C2A5"   # trait composition (CWM)
col_function     <- "#3B9AB2"   # ecosystem function

graph_code <- "
digraph apriori_psem {

  graph [compound = true, rankdir = TB, splines = spline,
         nodesep = 0.30, ranksep = '0.30:1.00:0.90:1.00', fontname = 'Helvetica']
  node  [shape = box, style = 'filled,rounded', fontname = 'Helvetica',
         fontsize = 11, penwidth = 0.8, color = 'grey35', margin = '0.14,0.07']
  edge  [fontname = 'Helvetica', fontsize = 10, color = 'grey30', penwidth = 1.1]

  # TIER 1: drivers
  subgraph cluster_drivers {
    label = 'Environmental and anthropogenic predictors'
    labelloc = 't'; fontsize = 12; fontname = 'Helvetica-Bold'
    style = 'rounded'; color = 'grey60'; penwidth = 1

    # row 1 - environmental
    subgraph {
      rank = same
      depth    [label = 'Depth',              fillcolor = '<ENV>']
      sst_m    [label = 'SST\\nmean',           fillcolor = '<ENV>']
      sst_sd   [label = 'SST\\nvariability',    fillcolor = '<ENV>']
      sbt_m    [label = 'SBT\\nmean',           fillcolor = '<ENV>']
      so_m     [label = 'Salinity\\nmean',      fillcolor = '<ENV>']
      so_sd    [label = 'Salinity\\nvariability', fillcolor = '<ENV>']
      chla     [label = 'Chlorophyll-a\\nmean',      fillcolor = '<ENV>']
    }

    # row 2 - anthropogenic
    subgraph {
      rank = same
      otb      [label = 'Fishing pressure (OTB)', fillcolor = '<ANT>']
      ssf      [label = 'Fishing pressure (SSF)',  fillcolor = '<ANT>']
      grav     [label = 'Human gravity index', fillcolor = '<ANT>']
    }

    depth -> sst_m -> sst_sd -> sbt_m -> so_m -> so_sd -> chla [style = invis]
    otb -> ssf -> grav [style = invis]
    sst_sd -> otb  [style = invis]
    so_sd  -> grav [style = invis]
  }

  # TIER 2: biodiversity
  subgraph cluster_biodiv {
    label = 'Biodiversity metrics'
    labelloc = 't'; fontsize = 12; fontname = 'Helvetica-Bold'
    style = 'rounded'; color = 'grey60'; penwidth = 1

    subgraph {
      rank = same
      tric [label = 'Taxonomic\\nrichness',      fillcolor = '<TAX>']
      teve [label = 'Taxonomic\\nevenness',      fillcolor = '<TAX>']
      cwml [label = 'CWM maximum\\nbody length (log)', fillcolor = '<TRA>']
      tric -> teve -> cwml [style = invis]
    }

    subgraph {
      rank = same
      fric [label = 'Functional\\nrichness', fillcolor = '<FUN>']
      feve [label = 'Functional\\nevenness', fillcolor = '<FUN>']
      fric -> feve [style = invis]
    }
  }

  # TIER 3: functioning
  biomass [label = 'Total fish biomass (log)',
           fillcolor = '<FUNC>', fontsize = 12, penwidth = 1.4]

  # legend 
  legend [shape = box, style = 'filled,rounded', fillcolor = 'white',
          color = 'grey70', penwidth = 1, fontsize = 10, margin = '0.20,0.14',
          label = 'Expected sign\\l \\l+      positive path predicted\\l-      negative path predicted\\l+/-   competing mechanisms, no directional prediction\\l- - -  residual covariance (no direction)\\l']

  # driver-level pathways 
  ssf -> teve [ltail = cluster_drivers, lhead = cluster_biodiv,
               penwidth = 2.2, color = 'grey45',
               label = '  all predictors -->\\n  all biodiversity metrics  ']

  # one arrow drivers -> biomass (direct paths)
  grav -> biomass [ltail = cluster_drivers,
                   penwidth = 2.2, color = 'grey45', style = 'solid',
                   label = '  direct paths  ']

  # within-biodiversity paths
  tric -> fric [label = ' + ', penwidth = 1.6]
  teve -> feve [label = ' + ', penwidth = 1.6]

  # residual covariances, no directional hypothesis
  edge [dir = both, style = dashed, color = 'grey55', penwidth = 0.9,
        label = '', constraint = false]
  tric -> teve
  fric -> feve
  tric -> feve
  # teve -> fric
  cwml -> fric
  cwml -> feve
  cwml -> teve

  # biodiversity -> function
  edge [dir = forward, style = solid, color = 'grey20', penwidth = 1.4,
        constraint = true]
  tric -> biomass [label = ' + ']
  teve -> biomass [label = ' +/- ']
  fric -> biomass [label = ' + ']
  feve -> biomass [label = ' +/- ']
  cwml -> biomass [label = ' +/- ']

  # legend anchored in layout
  { rank = same; biomass; legend }
  biomass -> legend [style = invis]
}
"

graph_code <- gsub("<ENV>",  col_driver_env,    graph_code, fixed = TRUE)
graph_code <- gsub("<ANT>",  col_driver_anthro, graph_code, fixed = TRUE)
graph_code <- gsub("<TAX>",  col_taxo,          graph_code, fixed = TRUE)
graph_code <- gsub("<FUN>",  col_func,          graph_code, fixed = TRUE)
graph_code <- gsub("<TRA>",  col_trait,         graph_code, fixed = TRUE)
graph_code <- gsub("<FUNC>", col_function,      graph_code, fixed = TRUE)

g <- DiagrammeR::grViz(graph_code)

# save
svg <- DiagrammeRsvg::export_svg(g)
# rsvg::rsvg_png(charToRaw(svg), file = here::here("Figures/Figure_S3.png"), width = 2400)
# rsvg::rsvg_pdf(charToRaw(svg), file = here::here("Figures/Figure_S3_apriori_pSEM_structure.pdf"))