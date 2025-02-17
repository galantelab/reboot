# WELCOME TO REBOOT (regression and survival tool with a multivariate bootstrap approach)
## A flexible and easy-to-use algorithm to identify and validate genes / transcripts expression signatures and associate it with patient survival. Reboot innovates by using a multivariate strategy with penalized Cox regression (Lasso method) combined with a bootstrap approach. Several statistical tests and images for signature score are generated.
### documentation [^1]: https://galantelab.github.io/reboot
### web tool: https://www.bioinfo.mochsl.org.br/reboot

<div align="center">
  <img src="docs/New_Reboot.inkscape.v2.300dpi-min.jpeg" alt="TOOL LOGO" width="500" height="500"/>
</div>

### Run Reboot in 4 easy STEPS:
1. **Download** this repository either directly or via our Docker container
2. **Install** Reboot and its dependencies
3. **Prepare** your input files
4. **Run** Reboot modules!

### Quick usage (CMDs):
- `git clone https://github.com/galantelab/reboot.git` (direct) | `docker pull galantelab/reboot` ([Docker](https://docs.docker.com/get-started/get-docker/) must be installed)
- `sudo sh reboot/install.sh` (if 'direct' installation was chosen, dependencies must be manually installed)
- `reboot.R complete <options>` (direct) | `docker run --rm -v $(pwd):$(pwd) galantelab/reboot reboot.R complete <options>`

**PS:** for a complete step by step walkthrough, please refer to our detailed guide at our [documentation](https://galantelab.github.io/reboot/).

### Citation:
If you use either the GUI or the CLI of Reboot, please cite us: **Reboot: a straightforward approach to identify genes and splicing isoforms associated with cancer patient prognosis.** Felipe R. C. dos Santos, Gabriela D. A. Guardia, Filipe F. dos Santos, Daniel T. Ohara, and Pedro A. F. Galante. *NAR Cancer*. DOI: [10.1093/narcan/zcab024](https://doi.org/10.1093/narcan/zcab024). PMID: [34316711](https://pubmed.ncbi.nlm.nih.gov/34316711/).

## CONTACT
For more information, you are welcome to visit us at our [Lab's website](https://www.bioinfo.mochsl.org.br/)!
[^1]: You can contact any of the authors by email.
