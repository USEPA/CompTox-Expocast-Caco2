# Supplemental Material

## Impact of Gut Permeability on Estimation of Oral Bioavailability for Chemicals in Commerce and the Environment

### Abstract

Administered equivalent dose (AED) estimation using *in vitro* hazard and high 
throughput toxicokinetics (HTTK) can be improved by refining assumptions 
regarding fraction absorbed ($F_{abs}$) through the intestine, a component of 
oral bioavailability ($F_{bio}$). Although *in vivo* data to inform $F_{abs}$ 
are often unavailable for non-pharmaceuticals, *in vitro* measures of apparent 
permeability ($P_{app}$) using the Caco-2 cell line have been highly correlated 
with $F_{abs}$ when used in *in vitro-in vivo* extrapolation (IVIVE) modeling. 
To address data gaps for non-drug chemicals, bidirectional $P_{app}$ was 
evaluated using the Caco-2 assay for over 400 chemicals. A random forest 
quantitative structure-property relationship (QSPR) model was developed using 
these and peer-reviewed pharmaceutical data. Both Caco-2 data ($R^2$=0.37) and 
the QSPR model ($R^2$=0.3) were better at predicting human bioavailability 
compared to *in vivo* rat data ($R^2$=0.2). The httk-predicted plasma steady 
state concentrations ($C_{ss}$) for IVIVE were updated, leading to modest 
changes for poorly absorbed chemicals. Experimental data were evaluated for 
sources of measurement uncertainty, which was then accounted for using the Monte 
Carlo method. Revised AEDs were subsequently compared with exposure estimates to
evaluate influence on bioactivity:exposure ratios as a surrogate for risk. 
Ultimately, $P_{app}$ incorporation to adjust for $F_{bio}$ in httk modeling 
improves AED estimations used in HT risk prioritization. 

### Data Analysis Vignettes

The following [R](https://cran.r-project.org/ "R") [vignettes](https://r-pkgs.org/vignettes.html "Vignettes") were written in [R Markdown](https://rmarkdown.rstudio.com/ "R Markdown") using [RStudio](https://posit.co/downloads/ "Download RStudio"):

1. [Caco-2 Data Organization and Analysis](https://github.com/USEPA/comptox-expocast-caco2/blob/main/Vignette1-MakeDataFigures.Rmd "Caco-2 Data").
2. [Machine Learning QSPR for Caco-2 Permeability](https://github.com/USEPA/comptox-expocast-caco2/blob/main/Vignette3-CreateQSPR.Rmd "Machine Learning QSPR").
3. [Using QSPR for New Chemicals](https://github.com/USEPA/comptox-expocast-caco2/blob/main/Vignette3-UseCaco2QSPRforNewChemicals.Rmd "Using QSPR for New Chemicals").
4. [Comparison of Oral Absorption Predictions to In Vivo Data](https://github.com/USEPA/comptox-expocast-caco2/blob/main/Vignette3-MakeEvaluationFigures.Rmd "Comparison to In Vivo Data").
5. [Impact of Bioavailability on Risk Prioritization](https://github.com/USEPA/comptox-expocast-caco2/blob/main/Vignette5-BioactivityExposureRatio.Rmd "Risk Prioritization").

### Authors

#### Principal Investigators 
John Wambaugh [wambaugh.john@epa.gov]

Elaina Kenyon [kenyon.elaina@epa.gov]

Barbara Wetmore [wetmore.barbara@epa.gov]

#### Full Author List

Gregory Honda<sup>1</sup>

Elaina M. Kenyon<sup>1</sup>

Sarah Davidson-Fritz<sup>1</sup>

Roger Dinallo<sup>2</sup>

Hisham El-Masri<sup>1</sup>

Evgenia Korel-Bexell<sup>1</sup>

Li Li<sup>2</sup>

Robert Pearce<sup>1</sup>

Risa Sayre<sup>1</sup>

Christopher Strock<sup>2</sup>

Russell Thomas<sup>1</sup>

John F. Wambaugh<sup>1</sup>

Barbara A. Wetmore<sup>1</sup>

1.	U.S. Environmental Protection Agency, Office of Research and Development, Center for Computational Toxicology and Exposure, Research Triangle Park, NC 27711, USA
2.	Cyprotex, Framingham, MA, USA

### License

License: GPL-3 <https://www.gnu.org/licenses/gpl-3.0.en.html>