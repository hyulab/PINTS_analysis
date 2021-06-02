This code repo hosts scripts and analysis protocols used to generate results that are reported in our study: *Systematic comparison of experimental assays and analytical pipelines for identification of active ­­­­enhancers genome-wide*

## Structure of this repo
### Protocols
This folder stores major protocols/pipelines that we used to analyze sequencing data. Files are in json format; they are structured as follows:
* `name`
* `description`
* `step`: a list of steps in this protocol/pipeline
    * `step_order`: the order of current step in this protocol
    * `software`
    * `parameter`
* `reference`: a dictionary of predefined variables that are referred to in this protocol
    * key: name of the reference
    * value: description about this reference

### Analysis
Scripts that we used to analyze data generated from protocols and produce aggregated files for downstream plotting.
* `assay_sensitivity.py`: Analysis about assays' sensitivity in detecting eRNAs.
* `assay_specificity.py`: 
    * Analysis about factors which affect assays' sensitivity.
    * Signal-to-noise ratio.
    * False discovery rates.
* `cross_analysis_via_read_distribution.py`: Analysis about the effects of rRNA.
* `tools_resolution_robustness.py`: Analysis about tools' applicability to different assays, resolution, and robustness.
* `tools_sensitivity_specificity.py`: Generate ROCs for different tools
* `compendium.py`: Gather summaries about the compendium
* `utils.py`: misc functions that are used in multiple modules
* `utils_roc.py`: python classes for ROC profilers

### Software
To incorporate previously published tools into our pipeline, sometimes we have to modify these tools a little bit. We put all our modifications in this folder.
