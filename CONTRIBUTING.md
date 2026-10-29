# Contributing to cellweave

The coding rules are in `CLAUDE.md`: English everywhere, the crate graph, the error handling, the tests, and the course in `docs/`, which every change to the model has to update. This page adds what is particular to data and to calibration.

## No data in the repository

The repository holds no data of anyone else's, whatever its licence. It holds code, configurations, and recipes.

- **A dataset enters as a recipe**, `datasets/<name>.toml`: the paper it comes with, the licence its authors give it, where each file is published, the sha256 fingerprint each file must have, and how `scripts/fetch-dataset.py` turns them into the tables the comparison scripts read. Whoever runs the script downloads the data for their own use, from the source, into `data/`, which git ignores. The CI fails if a data file is committed or a recipe leaves something out.
- **Fetch from where the authors published**: the journal's supplementary files, Europe PMC, Zenodo, the authors' own repository at a fixed commit. A direct link to the file, not scraping a publisher's pages.
- **Numbers quoted from a paper**, with the citation, are fine in the course and in the comments of a configuration: that is citing, not copying.
- **Your own data** can be compared with the model by pointing the scripts at any directory. What you have the right to use is yours to judge.
- **Never personal data.** Nothing about an identifiable person, and no individual-level human data, in a recipe, a configuration, a test, an issue or a pull request.

## Contributing a calibration

A calibration is a set of parameter values and the evidence for them. A pull request that brings one includes:

1. **The values**, in a configuration or network file under `examples/`, each changed value with a comment saying where it comes from.
2. **The data**: a recipe in `datasets/` when the data are openly published; otherwise the reference, and how access to them is obtained.
3. **The method**: which parameters were free, over which ranges, what was compared with what, and how the best values were chosen.
4. **The fit**: the whole comparison, what does not fit as well as what does, and the predictions written before fitting, if there were any.
5. **The chapter** of the course that tells it, under the rules of `CLAUDE.md`.

**Data under restricted access**, such as human data held by EGA or dbGaP, are welcome as the basis of a calibration, labelled as one the maintainers cannot reproduce. The contributor confirms that their data access agreement allows sharing what the pull request contains, and that it holds aggregate results only. The defaults the repository ships stay calibrated on openly published data, so that anyone can check them.
