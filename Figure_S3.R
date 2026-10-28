# ==============================================================================
# Figure S_3 — A priori conceptual path structure of the piecewise SEM
#
# Draws the hypothesised (pre-fitting) DAG, not the fitted model.
# Three tiers:
#   (1) exogenous environmental + anthropogenic drivers
#   (2) endogenous biodiversity variables (taxonomic, functional, trait)
#   (3) ecosystem function (log total fish biomass)
#
# Driver -> biodiversity and driver -> biomass paths are drawn as bundled
# cluster-level arrows (compound = true), because every driver was allowed to
# act on every biodiversity variable and on biomass; drawing the 60 individual
# arrows would be unreadable. Within-tier and biodiversity -> biomass paths,
# which carry the hypothesis-critical predictions, are drawn explicitly and
# labelled with their expected sign (see Tables S_ and S_).
#
# Layout: the exogenous tier is split into two rows (environmental drivers on
# the first row, anthropogenic drivers on the second) via rank = same groups;
# the legend is a plain node pinned to the biomass rank so it sits in the
# bottom-right whitespace of the figure.
#
# Author: Fabien Moullec
# ==============================================================================

library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

# ------------------------------------------------------------------ palette --
col_driver_env   <- "#BFD8E8"   # environmental drivers
col_driver_anthro<- "#F2C4A0"   # anthropogenic drivers
col_taxo         <- "#FC8D62"   # taxonomic diversity
col_func         <- "#8DA0CB"   # functional diversity
col_trait        <- "#66C2A5"   # trait composition (CWM)
col_function     <- "#3B9AB2"   # ecosystem function

graph_code <- "
digraph apriori_psem {

  # ranksep is a per-rank list: tight gap between the two driver rows,
  # wider gaps between the three conceptual tiers
  graph [compound = true, rankdir = TB, splines = spline,
         nodesep = 0.30, ranksep = '0.30:1.00:0.90:1.00', fontname = 'Helvetica']
  node  [shape = box, style = 'filled,rounded', fontname = 'Helvetica',
         fontsize = 11, penwidth = 0.8, color = 'grey35', margin = '0.14,0.07']
  edge  [fontname = 'Helvetica', fontsize = 10, color = 'grey30', penwidth = 1.1]

  # ---------------------------------------------------------- TIER 1: drivers
  subgraph cluster_drivers {
    label = 'Environmental and anthropogenic predictors'
    labelloc = 't'; fontsize = 12; fontname = 'Helvetica-Bold'
    style = 'rounded'; color = 'grey60'; penwidth = 1

    # row 1 — environmental
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

    # row 2 — anthropogenic
    subgraph {
      rank = same
      otb      [label = 'Fishing pressure (OTB)', fillcolor = '<ANT>']
      ssf      [label = 'Fishing pressure (SSF)',  fillcolor = '<ANT>']
      grav     [label = 'Human gravity index', fillcolor = '<ANT>']
    }

    # invisible edges: flat chains fix the left-to-right order within each row,
    # vertical ones stack row 2 underneath row 1
    depth -> sst_m -> sst_sd -> sbt_m -> so_m -> so_sd -> chla [style = invis]
    otb -> ssf -> grav [style = invis]
    sst_sd -> otb  [style = invis]
    so_sd  -> grav [style = invis]
  }

  # ----------------------------------------------------- TIER 2: biodiversity
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

  # ------------------------------------------------------- TIER 3: functioning
  biomass [label = 'Total fish biomass (log)',
           fillcolor = '<FUNC>', fontsize = 12, penwidth = 1.4]

  # ------------------------------------------------------------------ legend --
  # plain node (not a cluster) so it can share the biomass rank
  legend [shape = box, style = 'filled,rounded', fillcolor = 'white',
          color = 'grey70', penwidth = 1, fontsize = 10, margin = '0.20,0.14',
          label = 'Expected sign\\l \\l+      positive path predicted\\l-      negative path predicted\\l+/-   competing mechanisms, no directional prediction\\l- - -  residual covariance (no direction)\\l']

  # ------------------------------------------- bundled driver-level pathways --
  # one arrow drivers -> biodiversity block (tail on the lower driver row so the
  # biodiversity cluster is ranked below the whole driver cluster)
  ssf -> teve [ltail = cluster_drivers, lhead = cluster_biodiv,
               penwidth = 2.2, color = 'grey45',
               label = '  all predictors -->\\n  all biodiversity metrics  ']

  # one arrow drivers -> biomass (direct paths)
  grav -> biomass [ltail = cluster_drivers,
                   penwidth = 2.2, color = 'grey45', style = 'solid',
                   label = '  direct paths  ']

  # ------------------------------------------------ within-biodiversity paths --
  tric -> fric [label = ' + ', penwidth = 1.6]
  teve -> feve [label = ' + ', penwidth = 1.6]

  # residual (unresolved) covariances, no directional hypothesis
  # constraint = false: they must not influence rank assignment
  edge [dir = both, style = dashed, color = 'grey55', penwidth = 0.9,
        label = '', constraint = false]
  tric -> teve
  fric -> feve
  tric -> feve
  # teve -> fric
  cwml -> fric
  cwml -> feve
  cwml -> teve

  # ------------------------------------------------ biodiversity -> function --
  edge [dir = forward, style = solid, color = 'grey20', penwidth = 1.4,
        constraint = true]
  tric -> biomass [label = ' + ']
  teve -> biomass [label = ' +/- ']
  fric -> biomass [label = ' + ']
  feve -> biomass [label = ' +/- ']
  cwml -> biomass [label = ' +/- ']

  # ------------------------------------------------- legend anchored in layout --
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
g

# ------------------------------------------------------------------- export --
svg <- DiagrammeRsvg::export_svg(g)
rsvg::rsvg_png(charToRaw(svg), file = "Figures/Figure_S3.png", width = 2400)
rsvg::rsvg_pdf(charToRaw(svg), file = "Figures/Figure_S3_apriori_pSEM_structure.pdf")

# ==============================================================================
# NOTE
# Every driver -> biodiversity and driver -> biomass path shown here as a bundled
# arrow was offered to AIC-based selection within each sub-model; the retained
# subsets are reported in Table S_. Geographical Sub-Area (GSA) entered every
# sub-model as a main effect and in interaction with the biodiversity metrics
# (and with taxonomic richness in the functional richness sub-model, and with
# taxonomic evenness in the functional evenness sub-model); these moderating
# terms are omitted from the figure for legibility and are described in the
# Methods.
# ==============================================================================
